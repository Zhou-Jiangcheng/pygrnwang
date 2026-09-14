import os
import sys
import platform
import shutil
import subprocess
import tempfile
from functools import lru_cache

import numpy as np
from concurrent.futures import ProcessPoolExecutor


# ============================================================================
# Backend selection
# ----------------------------------------------------------------------------
# Java is invoked in a subprocess only when a travel-time query is requested.
# Without a JDK (java + javac), use the ObsPy backend.
# ============================================================================
def _detect_java_backend():
    """Locate the bundled JAR, with the environment's Scripts/bin as fallback."""
    if shutil.which("java") is None or shutil.which("javac") is None:
        return False, None
    candidates = [
        os.path.join(os.path.dirname(__file__), "exec", "TauP.jar"),
        os.path.join(
            sys.exec_prefix,
            "Scripts" if platform.system() == "Windows" else "bin",
            "TauP.jar",
        ),
    ]
    for jar_path in candidates:
        if os.path.isfile(jar_path):
            return True, os.path.abspath(jar_path)
    return False, None


_USE_JAVA, _JAR_PATH = _detect_java_backend()


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
        from obspy.taup import TauPyModel

        model_instance = TauPyModel(model=real_model_path)
        _MODEL_CACHE[model_name] = model_instance
        return model_instance
    except Exception as e:
        raise RuntimeError(f"Error loading model '{real_model_path}': {e}") from e


def _taup_create_npz_file_obspy(nd_file):
    npz_file = os.path.splitext(nd_file)[0] + ".npz"
    try:
        from obspy.taup.taup_create import TauPCreate

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
        raise RuntimeError(f"Model build failed for '{nd_file}': {e}") from e
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
@lru_cache(maxsize=4)
def _compile_java_bridge(jar_path):
    """Compile once per process into a private directory, removed on exit."""
    source = os.path.join(os.path.dirname(__file__), "java", "PygrnwangTauPBatch.java")
    build = tempfile.TemporaryDirectory(prefix="pygrnwang-taup-")
    try:
        completed = subprocess.run(
            ["javac", "-encoding", "UTF-8", "-cp", jar_path, "-d", build.name, source],
            capture_output=True,
            text=True,
            encoding="utf-8",
            check=False,
        )
        if completed.returncode:
            raise RuntimeError(f"Could not compile TauP bridge:\n{completed.stderr}")
    except Exception:
        build.cleanup()
        raise
    return build


def _query_java_batch(
    event_depth_km, distances_km, phase_groups, receiver_depth_km, model_name,
    first_only=False,
):
    """Query all distances and phase groups in one Java subprocess."""
    if not _USE_JAVA:
        raise RuntimeError("Java TauP requires java and javac on PATH and TauP.jar.")
    distances = np.asarray(distances_km, dtype=float).reshape(-1)
    results = [
        [dict(phase=[], puristphase=[], time=[], rayparameter=[]) for _ in phase_groups]
        for _ in distances
    ]
    if not len(distances):
        return results
    build = _compile_java_bridge(_JAR_PATH)
    completed = subprocess.run(
        [
            "java", "-cp", os.pathsep.join([build.name, _JAR_PATH]),
            "PygrnwangTauPBatch", os.fspath(model_name),
            str(float(event_depth_km)), str(float(receiver_depth_km)),
            "first" if first_only else "all",
            *[",".join(phases) for phases in phase_groups],
        ],
        input="".join(f"{dist * _DEG_PER_KM:.17g}\n" for dist in distances),
        capture_output=True,
        text=True,
        encoding="utf-8",
        check=False,
    )
    if completed.returncode:
        raise RuntimeError(f"TauP query failed:\n{completed.stderr}")
    seen = set()
    try:
        for line in completed.stdout.splitlines():
            index, group, phase, purist, time, rayparam = line.split("\t")
            index, group = int(index), int(group)
            if not (0 <= index < len(distances) and 0 <= group < len(phase_groups)):
                raise ValueError("Unexpected result index")
            result = results[index][group]
            seen.add((index, group))
            if phase:
                result["phase"].append(phase)
                result["puristphase"].append(purist)
                result["time"].append(float(time))
                result["rayparameter"].append(float(rayparam))
        if len(seen) != len(distances) * len(phase_groups):
            raise ValueError("Missing travel-time results")
    except (ValueError, IndexError) as exc:
        raise RuntimeError("TauP returned malformed or incomplete results") from exc
    return results


def taup_time_java(
    event_depth_km, dist_km, phases_list, receiver_depth_km=0, model_name="ak135"
):
    """Query all requested phase arrivals directly through Java TauP.

    Parameters
    ----------
    event_depth_km : float
        Requested source depth in km, positive down.
    dist_km : float
        Epicentral distance in km; query within the stored distance grid.
    phases_list : list of str
        TauP phase names to request, for example ['P', 'p', 'S'].
    receiver_depth_km : float, optional
        Requested receiver depth in km, positive down. Default: 0.
    model_name : str, optional
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library. Default: 'ak135'.

    Returns
    -------
    arrivals : dict
        Lists under phase, puristphase, time and rayparameter. Times are seconds
        and ray parameters are seconds/radian, not seconds/degree. Missing
        phases produce empty lists.

    Raises
    ------
    RuntimeError
        The model cannot be built or loaded, or the Java bridge fails.
    OSError
        A required file or executable cannot be accessed.

    Notes
    -----
    Requires java and javac plus the bundled JAR. The bridge is compiled lazily. Built-in TauP names and .nd model paths are supported. This explicitly Java-only entry does not fall back to ObsPy.
    """
    return _query_java_batch(
        event_depth_km, [dist_km], [phases_list], receiver_depth_km, model_name
    )[0][0]


def _first_arrival_java(event_depth_km, dist_km, receiver_depth_km, model_name, phases):
    result = _query_java_batch(
        event_depth_km, [dist_km], [phases], receiver_depth_km, model_name,
        first_only=True,
    )[0][0]
    return result["time"][0] if result["time"] else np.nan


def _cal_first_p_java(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    return _first_arrival_java(
        event_depth_km, dist_km, receiver_depth_km, model_name, _PHASES_P
    )


def _cal_first_s_java(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    return _first_arrival_java(
        event_depth_km, dist_km, receiver_depth_km, model_name, _PHASES_S
    )


def _first_p_s_java_batch(event_depth_km, distances_km, receiver_depth_km, model_name):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    results = _query_java_batch(
        event_depth_km, distances_km, [_PHASES_P, _PHASES_S], receiver_depth_km,
        model_name, first_only=True,
    )
    first_p = np.array([r[0]["time"][0] if r[0]["time"] else np.nan for r in results])
    first_s = np.array([r[1]["time"][0] if r[1]["time"] else np.nan for r in results])
    return first_p, first_s


def _cal_first_p_s_java(
    event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"
):
    first_p, first_s = _first_p_s_java_batch(
        event_depth_km, [dist_km], receiver_depth_km, model_name
    )
    return float(first_p[0]), float(first_s[0])


# ============================================================================
# Public API (dispatches to the selected backend)
# ============================================================================
def taup_create_npz_file(nd_file):
    """Prepare a custom velocity-model path for the selected TauP backend.

    Parameters
    ----------
    nd_file : str
        Path to a named-discontinuity velocity model; the selected TauP backend reads or converts this file.

    Returns
    -------
    model_path : str
        Java returns nd_file unchanged; ObsPy builds and returns a sibling .npz path.

    Raises
    ------
    RuntimeError
        The model cannot be built or loaded, or the Java bridge fails.
    OSError
        A required file or executable cannot be accessed.

    Notes
    -----
    Uses Java subprocesses when the bundled JAR, java and javac are available; otherwise uses ObsPy. Importing this module does not start a JVM. Runtime Java failures are reported rather than silently switching backends. Distances use 111.19492664455874 km per degree. The historical function name does not imply that every backend creates an NPZ file.
    """
    if _USE_JAVA:
        return nd_file
    return _taup_create_npz_file_obspy(nd_file)


def cal_first_p(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """Calculate the earliest arrival from the configured P phase set.

    Parameters
    ----------
    event_depth_km : float
        Requested source depth in km, positive down.
    dist_km : float
        Epicentral distance in km; query within the stored distance grid.
    receiver_depth_km : float, optional
        Requested receiver depth in km, positive down. Default: 0.0.
    model_name : str, optional
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library. Default: 'ak135'.

    Returns
    -------
    arrival : float
        Arrival time in seconds after origin; NaN when the selected phase set has no arrival.

    Raises
    ------
    RuntimeError
        The model cannot be built or loaded, or the Java bridge fails.
    OSError
        A required file or executable cannot be accessed.

    Notes
    -----
    Uses Java subprocesses when the bundled JAR, java and javac are available; otherwise uses ObsPy. Importing this module does not start a JVM. Runtime Java failures are reported rather than silently switching backends. Distances use 111.19492664455874 km per degree. Phases: p, P, pP, Pg, Pn, Pdiff, PKP. The deeper endpoint is treated as source using reciprocity.
    """
    if _USE_JAVA:
        return _cal_first_p_java(event_depth_km, dist_km, receiver_depth_km, model_name)
    else:
        return _cal_first_p_obspy(
            event_depth_km, dist_km, receiver_depth_km, model_name
        )


def cal_first_s(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """Calculate the earliest arrival from the configured S phase set.

    Parameters
    ----------
    event_depth_km : float
        Requested source depth in km, positive down.
    dist_km : float
        Epicentral distance in km; query within the stored distance grid.
    receiver_depth_km : float, optional
        Requested receiver depth in km, positive down. Default: 0.0.
    model_name : str, optional
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library. Default: 'ak135'.

    Returns
    -------
    arrival : float
        Arrival time in seconds after origin; NaN when the selected phase set has no arrival.

    Raises
    ------
    RuntimeError
        The model cannot be built or loaded, or the Java bridge fails.
    OSError
        A required file or executable cannot be accessed.

    Notes
    -----
    Uses Java subprocesses when the bundled JAR, java and javac are available; otherwise uses ObsPy. Importing this module does not start a JVM. Runtime Java failures are reported rather than silently switching backends. Distances use 111.19492664455874 km per degree. Phases: s, S, sS, pS, Sg, Sn, Sdiff, SKS. The deeper endpoint is treated as source using reciprocity.
    """
    if _USE_JAVA:
        return _cal_first_s_java(event_depth_km, dist_km, receiver_depth_km, model_name)
    else:
        return _cal_first_p_s_obspy(
            event_depth_km, dist_km, receiver_depth_km, model_name
        )[1]


def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """Calculate first P and S arrivals for one geometry.

    Parameters
    ----------
    event_depth_km : float
        Requested source depth in km, positive down.
    dist_km : float
        Epicentral distance in km; query within the stored distance grid.
    receiver_depth_km : float, optional
        Requested receiver depth in km, positive down. Default: 0.0.
    model_name : str, optional
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library. Default: 'ak135'.

    Returns
    -------
    first_p, first_s : float
        P and S arrival seconds after origin. Each can independently be NaN.

    Raises
    ------
    RuntimeError
        The model cannot be built or loaded, or the Java bridge fails.
    OSError
        A required file or executable cannot be accessed.

    Notes
    -----
    Uses Java subprocesses when the bundled JAR, java and javac are available; otherwise uses ObsPy. Importing this module does not start a JVM. Runtime Java failures are reported rather than silently switching backends. Distances use 111.19492664455874 km per degree. Uses the same phase sets and endpoint reciprocity as cal_first_p and cal_first_s.
    """
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
    """Write first P/S arrival tables for a depth pair and distance sequence.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    event_depth_km : float
        Requested source depth in km, positive down.
    receiver_depth_km : float
        Requested receiver depth in km, positive down.
    dist_km_list : sequence of float
        Epicentral distances in km; output order follows this sequence.
    model_name : str, optional
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library. Default: 'ak135'.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.
    max_workers : int or None, optional
        ObsPy process limit; None uses CPU count. Fewer than 50 distances use one worker; Java always batches in one process. Default: None.

    Returns
    -------
    None
        Writes tp_table.bin and ts_table.bin under source/receiver depth folders
        named with two decimal places. Each file is one float32 second value
        per input distance; NaN records a missing arrival.

    Raises
    ------
    RuntimeError
        The model cannot be built or loaded, or the Java bridge fails.
    OSError
        A required file or executable cannot be accessed.

    Notes
    -----
    Uses Java subprocesses when the bundled JAR, java and javac are available; otherwise uses ObsPy. Importing this module does not start a JVM. Runtime Java failures are reported rather than silently switching backends. Distances use 111.19492664455874 km per degree. Java evaluates all distances in one subprocess; ObsPy may split them among workers. Model-cache files are built before workers start. Guard multiprocessing calls with if __name__ == "__main__".
    """
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
        # Reuse one Java process and model for the entire distance table.
        tp_table, ts_table = _first_p_s_java_batch(
            event_depth_km, dist_km_list, receiver_depth_km, model_name
        )
        tp_table.astype(np.float32).tofile(path_tp_table)
        ts_table.astype(np.float32).tofile(path_ts_table)
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
