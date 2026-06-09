import os
import sys
import platform
import shutil

import numpy as np
from concurrent.futures import ProcessPoolExecutor


# ============================================================================
# Backend selection
# ----------------------------------------------------------------------------
# If a Java runtime, TauP.jar and jpype are all available, use the (faster)
# Java TauP backend. Otherwise fall back to obspy.taup, which is slower.
# ============================================================================
def _detect_java_backend():
    """
    Return (use_java, jar_path).
    use_java is True only when `java` is on PATH, TauP.jar exists in the
    environment's bin/Scripts directory and jpype can be imported.
    """
    if shutil.which("java") is None:
        return False, None

    if platform.system() == "Windows":
        jar_path = os.path.join(sys.exec_prefix, "Scripts", "TauP.jar")
    else:
        jar_path = os.path.join(sys.exec_prefix, "bin", "TauP.jar")

    if not os.path.exists(jar_path):
        return False, None

    try:
        import jpype  # noqa: F401
    except ImportError:
        return False, None

    return True, jar_path


_USE_JAVA, _JAR_PATH = _detect_java_backend()

if _USE_JAVA:
    import jpype
    import jpype.imports  # noqa: F401  use jpype to call java class

    if not jpype.isJVMStarted():
        jpype.startJVM("--enable-native-access=ALL-UNNAMED", classpath=[_JAR_PATH])
    from edu.sc.seis.TauP import TauP_Time  # type: ignore
else:
    from obspy.taup import TauPyModel
    from obspy.taup.taup_create import TauPCreate


_DEG_PER_KM = 1.0 / 111.19492664455874

# Phase lists shared by both backends.
_PHASES_P = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
_PHASES_S = ["s", "S", "sS", "pS", "Sg", "Sn", "Sdiff", "SKS"]


# ============================================================================
# obspy backend
# ============================================================================
# --- Global Model Cache ---
_MODEL_CACHE = {}


def _get_model(model_name, rebuild_npz=False):
    """
    Retrieve a cached TauPyModel instance (obspy backend only).
    """
    if model_name in _MODEL_CACHE:
        return _MODEL_CACHE[model_name]

    real_model_path = model_name

    # Auto-build logic: if input is an .nd file
    if model_name.endswith(".nd"):
        nd_file = model_name
        npz_file = os.path.splitext(nd_file)[0] + ".npz"

        # Check/Build .npz
        if not os.path.exists(npz_file) or rebuild_npz:
            # Only allow building in the main process logic or ensure file lock (simplified here)
            # The calling function ensures this is done before forking in most cases.
            _taup_create_npz_file_obspy(nd_file)

        real_model_path = npz_file

    try:
        model_instance = TauPyModel(model=real_model_path)
        _MODEL_CACHE[model_name] = model_instance
        return model_instance
    except Exception as e:
        print(f"Error loading model '{real_model_path}': {e}")
        sys.exit(1)


def _taup_create_npz_file_obspy(nd_file):
    npz_file = os.path.splitext(nd_file)[0] + ".npz"
    try:
        taup_creator = TauPCreate(
            input_filename=nd_file,
            output_filename=npz_file,
            verbose=False,
            # Using conservative parameters to avoid crashes in low-velocity zones
            # min_delta_p=0.05,
            # max_depth_interval=1.0,
            # max_interp_error=0.03
        )
        taup_creator.load_velocity_model()
        taup_creator.run()
    except Exception as e:
        print(f"Model build failed: {e}")
        sys.exit(1)
    return npz_file


def _cal_first_p_obspy(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    # Force deeper point as source (reciprocity)
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km

    # This will use the per-process cache
    model = _get_model(model_name)
    dist_deg = dist_km * _DEG_PER_KM

    try:
        arrivals_p = model.get_travel_times(
            source_depth_in_km=event_depth_km,
            receiver_depth_in_km=receiver_depth_km,
            distance_in_degree=dist_deg,
            phase_list=_PHASES_P,
        )
        first_p = arrivals_p[0].time if arrivals_p else np.nan
    except Exception:
        first_p = np.nan

    return first_p


def _cal_first_p_s_obspy(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    # Force deeper point as source (reciprocity)
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km

    # This will use the per-process cache
    model = _get_model(model_name)
    dist_deg = dist_km * _DEG_PER_KM

    # 1. First P
    try:
        arrivals_p = model.get_travel_times(
            source_depth_in_km=event_depth_km,
            receiver_depth_in_km=receiver_depth_km,
            distance_in_degree=dist_deg,
            phase_list=_PHASES_P,
        )
        first_p = arrivals_p[0].time if arrivals_p else np.nan
    except Exception:
        first_p = np.nan

    # 2. First S
    try:
        arrivals_s = model.get_travel_times(
            source_depth_in_km=event_depth_km,
            receiver_depth_in_km=receiver_depth_km,
            distance_in_degree=dist_deg,
            phase_list=_PHASES_S,
        )
        first_s = arrivals_s[0].time if arrivals_s else np.nan
    except Exception:
        first_s = np.nan

    return first_p, first_s


# ============================================================================
# Java (TauP) backend
# ============================================================================
def taup_time_java(
    event_depth_km, dist_km, phases_list, receiver_depth_km=0, model_name="ak135"
):
    """
    Full travel-time query using the Java TauP backend.
    Only available when the Java backend has been selected.
    """
    if not _USE_JAVA:
        raise RuntimeError(
            "taup_time_java requires a Java runtime with TauP.jar and jpype."
        )
    ttobj = TauP_Time(model_name)
    ttobj.setSourceDepth(event_depth_km)
    ttobj.setReceiverDepth(receiver_depth_km)
    ttobj.setPhaseNames(phases_list)
    ttobj.calculate(dist_km * _DEG_PER_KM)

    N_arr = ttobj.getNumArrivals()
    results = {"phase": [], "puristphase": [], "time": [], "rayparameter": []}
    for i in range(N_arr):
        arr = ttobj.getArrival(i)
        results["phase"].append(str(arr.getName()))
        results["puristphase"].append(str(arr.getPuristName()))
        results["time"].append(float(arr.getTime()))
        results["rayparameter"].append(float(arr.getRayParam()))
    return results


def _first_arrival_java(event_depth_km, dist_km, receiver_depth_km, model_name, phases):
    ttobj = TauP_Time(model_name)
    ttobj.setSourceDepth(event_depth_km)
    ttobj.setReceiverDepth(receiver_depth_km)
    ttobj.setPhaseNames(phases)
    ttobj.calculate(dist_km * _DEG_PER_KM)
    if ttobj.getNumArrivals() == 0:
        return np.nan
    return float(ttobj.getArrival(0).getTime())


def _cal_first_p_java(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    return _first_arrival_java(
        event_depth_km, dist_km, receiver_depth_km, model_name, _PHASES_P
    )


def _cal_first_p_s_java(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    first_p = _first_arrival_java(
        event_depth_km, dist_km, receiver_depth_km, model_name, _PHASES_P
    )
    first_s = _first_arrival_java(
        event_depth_km, dist_km, receiver_depth_km, model_name, _PHASES_S
    )
    return first_p, first_s


# ============================================================================
# Public API (dispatches to the selected backend)
# ============================================================================
def taup_create_npz_file(nd_file):
    """
    Prepare a velocity model from an .nd file.

    obspy backend: builds and returns the corresponding .npz file.
    Java backend : TauP reads .nd files directly, so the .nd path is returned
                   unchanged (no build step needed).
    """
    if _USE_JAVA:
        return nd_file
    return _taup_create_npz_file_obspy(nd_file)


def cal_first_p(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """Calculate the first P arrival time for a single distance."""
    if _USE_JAVA:
        return _cal_first_p_java(event_depth_km, dist_km, receiver_depth_km, model_name)
    else:
        return _cal_first_p_obspy(
            event_depth_km, dist_km, receiver_depth_km, model_name
        )


def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """Calculate the first P and S arrival times for a single distance."""
    if _USE_JAVA:
        return _cal_first_p_s_java(
            event_depth_km, dist_km, receiver_depth_km, model_name
        )
    else:
        return _cal_first_p_s_obspy(
            event_depth_km, dist_km, receiver_depth_km, model_name
        )


# --- Worker Function for Parallelization (obspy backend only) ---
def _calculate_chunk(dist_chunk, event_depth_km, receiver_depth_km, model_name):
    """
    Worker function to process a chunk of distances.
    Must be at the top level to be picklable by multiprocessing.
    """
    chunk_tp = []
    chunk_ts = []

    # Iterate through the subset of distances
    for dist in dist_chunk:
        fp, fs = cal_first_p_s(event_depth_km, dist, receiver_depth_km, model_name)
        chunk_tp.append(fp)
        chunk_ts.append(fs)

    return np.array(chunk_tp, dtype=np.float32), np.array(chunk_ts, dtype=np.float32)


def create_tpts_table(
    path_green,
    event_depth_km,
    receiver_depth_km,
    dist_km_list,
    model_name="ak135",
    check_finished=False,
    max_workers=None,  # Added parameter to control parallelism (obspy backend)
):
    # Ensure directory exists
    dir_path = os.path.join(
        path_green, "%.2f" % event_depth_km, "%.2f" % receiver_depth_km
    )
    if not os.path.exists(dir_path):
        os.makedirs(dir_path, exist_ok=True)

    path_tp_table = os.path.join(str(dir_path), "tp_table.bin")
    path_ts_table = os.path.join(str(dir_path), "ts_table.bin")

    if (
        check_finished
        and os.path.exists(path_tp_table)
        and os.path.exists(path_ts_table)
    ):
        return

    if _USE_JAVA:
        # The JVM does not play well with forked worker processes, so run
        # serially in the current (JVM-hosting) process.
        tp_table = np.zeros(len(dist_km_list), dtype=np.float32)
        ts_table = np.zeros(len(dist_km_list), dtype=np.float32)
        for i in range(len(dist_km_list)):
            first_p, first_s = _cal_first_p_s_java(
                event_depth_km, dist_km_list[i], receiver_depth_km, model_name
            )
            tp_table[i] = first_p
            ts_table[i] = first_s
        tp_table.tofile(path_tp_table)
        ts_table.tofile(path_ts_table)
        return

    # ---- obspy backend (parallel) ----
    # [CRITICAL] Pre-load/Build model in the MAIN process first.
    # This prevents a race condition where multiple workers try to build
    # the .npz file simultaneously if it doesn't exist.
    _get_model(model_name)

    # Determine number of workers (default to CPU count)
    if max_workers is None:
        max_workers = os.cpu_count() or 1

    # If data is small, don't use multiprocessing overhead
    if len(dist_km_list) < 50:
        max_workers = 1

    tp_table_parts = []
    ts_table_parts = []

    # Using ProcessPoolExecutor for parallel execution
    if max_workers > 1:
        # Split the distance list into chunks for each worker
        chunks = np.array_split(dist_km_list, max_workers)

        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit tasks
            futures = [
                executor.submit(
                    _calculate_chunk,
                    chunk,
                    event_depth_km,
                    receiver_depth_km,
                    model_name,
                )
                for chunk in chunks
            ]

            # Loop through futures in order of submission to keep order.
            for future in futures:
                res_tp, res_ts = future.result()
                tp_table_parts.append(res_tp)
                ts_table_parts.append(res_ts)
    else:
        # Serial execution fallback
        res_tp, res_ts = _calculate_chunk(
            dist_km_list, event_depth_km, receiver_depth_km, model_name
        )
        tp_table_parts.append(res_tp)
        ts_table_parts.append(res_ts)

    # Concatenate all parts
    tp_table = np.concatenate(tp_table_parts)
    ts_table = np.concatenate(ts_table_parts)

    tp_table.tofile(path_tp_table)
    ts_table.tofile(path_ts_table)


if __name__ == "__main__":
    pass