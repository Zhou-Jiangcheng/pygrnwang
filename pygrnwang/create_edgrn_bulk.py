import os
import pickle
import json
import datetime

from tqdm import tqdm
from multiprocessing import Pool

from .create_edgrn import create_inp_edgrn2, call_edgrn2
from .utils import group, convert_earth_model_nd2nd_without_Q


def _call_edgrn2_star(args):
    return call_edgrn2(*args)


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_edgrn2(
    processes_num,
    path_green,
    grn_source_depth_range,
    grn_source_delta_depth,
    grn_dist_range,
    grn_delta_dist,
    obs_depth_list,
    wavenumber_sampling_rate=12,
    path_nd=None,
    earth_model_layer_num=None,
):
    # print("preprocessing edgrn2")

    """Prepare the edgrn library grid, input files and job groups.

    Parameters
    ----------
    processes_num : int
        Positive worker count used to group jobs; MPI rank count must match the prepared group width.
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    grn_source_depth_range : list of float
        Minimum and maximum source depths in km.
    grn_source_delta_depth : float
        Positive source-depth grid increment in km.
    grn_dist_range : list of float
        Minimum and maximum epicentral distances in km.
    grn_delta_dist : float
        Positive regular epicentral-distance increment in km.
    obs_depth_list : list of float
        Nonempty receiver depth list in km, positive down.
    wavenumber_sampling_rate : float, optional
        Dimensionless spatial Nyquist oversampling factor for wavenumber integration. Default: 12.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.

    Returns
    -------
    group_list : list
        Jobs grouped by processes_num; the same groups are saved as a pickle file.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions.
    """
    for obs_depth in obs_depth_list:
        sub_sub_dir = str(os.path.join(path_green, "edgrn2", "%.2f" % obs_depth))
        os.makedirs(sub_sub_dir, exist_ok=True)
        create_inp_edgrn2(
            path_green,
            obs_depth,
            grn_dist_range,
            grn_delta_dist,
            grn_source_depth_range,
            grn_source_delta_depth,
            wavenumber_sampling_rate,
            path_nd,
            earth_model_layer_num,
        )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    green_info = {
        "processes_num": processes_num,
        "grn_source_depth_range": grn_source_depth_range,
        "grn_source_delta_depth": grn_source_delta_depth,
        "grn_dist_range": grn_dist_range,
        "grn_delta_dist": grn_delta_dist,
        "obs_depth_list": obs_depth_list,
        "wavenumber_sampling_rate": wavenumber_sampling_rate,
        "path_nd": path_nd,
        "path_nd_without_Q": path_nd_without_Q,
        "earth_model_layer_num": earth_model_layer_num,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    group_list_edgrn = group(obs_depth_list, processes_num)
    with open(os.path.join(path_green, "group_list_edgrn.pkl"), "wb") as fw:
        pickle.dump(group_list_edgrn, fw)  # type: ignore
    return group_list_edgrn


def create_grnlib_edgrn2_sequential(path_green, check_finished=False):
    """Compute the prepared edgrn library sequentially.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.

    Returns
    -------
    elapsed : datetime.timedelta
        Wall-clock duration of the computation loop.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edgrn.pkl"), "rb") as fr:
        group_list_edgrn = pickle.load(fr)
    for item in tqdm(group_list_edgrn, desc="Computing Green's function library"):
        for i in range(len(item)):
            # print("computing " + str(item[i]) + " km")
            call_edgrn2(item[i], path_green, check_finished)
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))
    return e - s


def create_grnlib_edgrn2_parallel(path_green, check_finished=False):
    """Compute the prepared edgrn library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.

    Returns
    -------
    elapsed : datetime.timedelta
        Wall-clock duration of the computation loop.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution. On Windows call under an if __name__ == "__main__" guard.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edgrn.pkl"), "rb") as fr:
        group_list_edgrn = pickle.load(fr)
    tasks = []
    for grp in group_list_edgrn:
        for d in grp:
            tasks.append((d, path_green, check_finished))

    processes = None
    try:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            processes = json.load(fr).get("processes_num", None)
    except Exception:
        pass

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_edgrn2_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Computing Green's function library",
        ):
            pass
    e = datetime.datetime.now()
    return e - s


def create_grnlib_edgrn2_parallel_multi_nodes(path_green, check_finished=False):
    """Compute the prepared edgrn library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.

    Returns
    -------
    elapsed : datetime.timedelta
        Wall-clock duration of the computation loop.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    RuntimeError
        mpi4py is unavailable.
    ValueError
        MPI rank count does not match the prepared group width.

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    with open(os.path.join(path_green, "group_list_edgrn.pkl"), "rb") as fr:
        group_list_edgrn = pickle.load(fr)
    for ind_group in range(len(group_list_edgrn)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num != len(group_list_edgrn[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!"
                % (processes_num, len(group_list_edgrn[0]))
            )
        print("ind_group:%d rank:%d" % (ind_group, rank))
        # the last group holds the remainder and may be shorter than processes_num
        if rank >= len(group_list_edgrn[ind_group]):
            continue
        call_edgrn2(
            obs_depth=group_list_edgrn[ind_group][rank],
            path_green=path_green,
            check_finished=check_finished,
        )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))
    return e - s
