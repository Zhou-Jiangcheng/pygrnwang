import os
import sys
import numpy as np
from obspy.taup import TauPyModel
from obspy.taup.taup_create import TauPCreate
from concurrent.futures import ProcessPoolExecutor  # 引入并行处理模块

# --- Global Model Cache ---
_MODEL_CACHE = {}


def _get_model(model_name, rebuild_npz=False):
    """
    Retrieve a cached TauPyModel instance.
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
            taup_create_npz_file(nd_file)

        real_model_path = npz_file

    try:
        model_instance = TauPyModel(model=real_model_path)
        _MODEL_CACHE[model_name] = model_instance
        return model_instance
    except Exception as e:
        print(f"Error loading model '{real_model_path}': {e}")
        sys.exit(1)


def taup_create_npz_file(nd_file):
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

def cal_first_p(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """
    Calculates P and S times for a single distance.
    This function remains unchanged and is called by workers.
    """
    # Force deeper point as source (reciprocity)
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km

    # This will use the per-process cache
    model = _get_model(model_name)
    dist_deg = dist_km / 111.19492664455874

    # 1. First P
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    try:
        arrivals_p = model.get_travel_times(
            source_depth_in_km=event_depth_km,
            receiver_depth_in_km=receiver_depth_km,
            distance_in_degree=dist_deg,
            phase_list=phases_list_p
        )
        first_p = arrivals_p[0].time if arrivals_p else np.nan
    except Exception:
        first_p = np.nan

    return first_p

def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0.0, model_name="ak135"):
    """
    Calculates P and S times for a single distance.
    This function remains unchanged and is called by workers.
    """
    # Force deeper point as source (reciprocity)
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km

    # This will use the per-process cache
    model = _get_model(model_name)
    dist_deg = dist_km / 111.19492664455874

    # 1. First P
    phases_list_p = ["p", "P", "pP", "Pg", "Pn", "Pdiff", "PKP"]
    try:
        arrivals_p = model.get_travel_times(
            source_depth_in_km=event_depth_km,
            receiver_depth_in_km=receiver_depth_km,
            distance_in_degree=dist_deg,
            phase_list=phases_list_p
        )
        first_p = arrivals_p[0].time if arrivals_p else np.nan
    except Exception:
        first_p = np.nan

    # 2. First S
    phases_list_s = ["s", "S", "sS", "pS", "Sg", "Sn", "Sdiff", "SKS"]
    try:
        arrivals_s = model.get_travel_times(
            source_depth_in_km=event_depth_km,
            receiver_depth_in_km=receiver_depth_km,
            distance_in_degree=dist_deg,
            phase_list=phases_list_s
        )
        first_s = arrivals_s[0].time if arrivals_s else np.nan
    except Exception:
        first_s = np.nan

    return first_p, first_s


# --- Worker Function for Parallelization ---
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
        max_workers=None  # Added parameter to control parallelism
):
    # Ensure directory exists
    dir_path = os.path.join(path_green, "%.2f" % event_depth_km, "%.2f" % receiver_depth_km)
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

    # print(f"Calculating table (Depth: {event_depth_km}km) with {max_workers} workers...")

    tp_table_parts = []
    ts_table_parts = []

    # Using ProcessPoolExecutor for parallel execution
    if max_workers > 1:
        # Split the distance list into chunks for each worker
        chunks = np.array_split(dist_km_list, max_workers)
        
        with ProcessPoolExecutor(max_workers=max_workers) as executor:
            # Submit tasks
            futures = [
                executor.submit(_calculate_chunk, chunk, event_depth_km, receiver_depth_km, model_name)
                for chunk in chunks
            ]
            
            # Collect results as they complete (order doesn't matter for collection, 
            # but we need to reassemble in order. Map or simple loop works best here.)
            # Here we just loop through futures in order of submission to keep order.
            for future in futures:
                res_tp, res_ts = future.result()
                tp_table_parts.append(res_tp)
                ts_table_parts.append(res_ts)
    else:
        # Serial execution fallback
        res_tp, res_ts = _calculate_chunk(dist_km_list, event_depth_km, receiver_depth_km, model_name)
        tp_table_parts.append(res_tp)
        ts_table_parts.append(res_ts)

    # Concatenate all parts
    tp_table = np.concatenate(tp_table_parts)
    ts_table = np.concatenate(ts_table_parts)

    tp_table.tofile(path_tp_table)
    ts_table.tofile(path_ts_table)


if __name__ == "__main__":
    # Example usage
    my_nd_file = r"C:\Users\zjc\my_data\wenchuan.nd"
    
    # Generate a dummy distance list
    dists = np.linspace(10, 2000, 500) # 500 points
    
    try:
        import time
        s = time.time()
        
        create_tpts_table(
            path_green="./output_tables",
            event_depth_km=10.0,
            receiver_depth_km=0.0,
            dist_km_list=dists,
            model_name=my_nd_file,
            max_workers=4  # Use 4 cores
        )
        
        print(f"Calculation finished in {time.time() - s:.2f}s")
        
    except Exception as e:
        import traceback
        traceback.print_exc()
        print(f"Error in main execution: {e}")