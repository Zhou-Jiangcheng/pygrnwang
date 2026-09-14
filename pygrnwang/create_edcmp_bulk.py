import os
import pickle
import json
import datetime
import warnings
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm

from .create_edcmp import create_inp_edcmp2, call_edcmp2, convert_edcmp2
from .utils import group, cal_grid


def _call_edcmp2_star(args):
    return call_edcmp2(*args)


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_edcmp2(
    processes_num: int,
    path_green: str,
    grn_source_depth_range,
    grn_source_delta_depth,
    grn_dist_range,
    grn_delta_dist,
    obs_depth_list,
    output_observables=(1, 0, 0, 0),
    layered=True,
    lam=30516224000,
    mu=33701888000,
):
    """Prepare the edcmp library grid, input files and job groups.

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
    output_observables : sequence of int, optional
        Four 0/1 flags in displacement, strain, stress, tilt order. Default: (1, 0, 0, 0).
    layered : bool, optional
        True uses the EDGRN layered-medium library; False selects the homogeneous half-space formula. Default: True.
    lam : float, optional
        First Lame parameter in Pa for the homogeneous half-space. Default: 30516224000.
    mu : float, optional
        Shear modulus in Pa for the homogeneous half-space. Default: 33701888000.

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
    See the edcmp tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. Run EDGRN preparation first: this function reads and updates its metadata, even for homogeneous-half-space mode.
    """
    # Load the current Green's function library information.
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)

    item_list = []
    event_depth_list = cal_grid(
        grn_source_depth_range[0],
        grn_source_depth_range[1],
        grn_source_delta_depth,
    )
    for event_depth in event_depth_list:
        for obs_depth in obs_depth_list:
            for mt_ind in range(5):
                create_inp_edcmp2(
                    path_green=path_green,
                    event_depth=event_depth,
                    obs_depth=obs_depth,
                    dist_range=grn_dist_range,
                    delta_dist=grn_delta_dist,
                    mt_ind=mt_ind,
                    output_observables=output_observables,
                    layered=layered,
                    lam=lam,
                    mu=mu,
                )
                item_list.append([event_depth, obs_depth, mt_ind])

    # Update the green_info dictionary with new observation and model parameters.
    # seek_edcmp2 indexes the edcmp2 output on the grid used here, so record it;
    # pre_process_edgrn2 wrote the edgrn2 grid into the same keys, and the two
    # disagreeing means the caller passed different parameters to the two steps.
    grid_keys = {
        "grn_source_depth_range": list(grn_source_depth_range),
        "grn_source_delta_depth": grn_source_delta_depth,
        "grn_dist_range": list(grn_dist_range),
        "grn_delta_dist": grn_delta_dist,
        "obs_depth_list": list(obs_depth_list),
    }
    for key, value in grid_keys.items():
        old = green_info.get(key, None)
        if old is not None and old != value:
            warnings.warn(
                "pre_process_edcmp2 got %s=%r but pre_process_edgrn2 recorded %r; "
                "the edcmp2 grid is the one seek_edcmp2 will use." % (key, value, old)
            )
        green_info[key] = value
    green_info["layered"] = layered
    if layered:
        green_info["lam"] = None
        green_info["mu"] = None
    else:
        green_info["lam"] = lam
        green_info["mu"] = mu
    green_info["output_observables"] = output_observables
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    # Group the observation depths into subgroups for parallel processing.
    group_list_edcmp = group(item_list, processes_num)
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "wb") as fw:
        pickle.dump(group_list_edcmp, fw)  # type: ignore

    return group_list_edcmp


def create_grnlib_edcmp2_sequential(path_green, check_finished=False):
    # s = datetime.datetime.now()
    """Compute the prepared edcmp library sequentially.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    for item in tqdm(group_list_edcmp, desc="Computing static stress"):
        for i in range(len(item)):
            # print("computing " + str(item[i]) + " km")
            item[i] = item[i] + [path_green, check_finished]
            call_edcmp2(*item[i])
    # e = datetime.datetime.now()
    # print("run time:%s" % str(e - s))


def create_grnlib_edcmp2_parallel(
    path_green, check_finished=False, convert_bulk=True, remove=False
):
    """Compute the prepared edcmp library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.
    convert_bulk : bool, optional
        Generate combined EDCMP float32 libraries after jobs finish. Default: True.
    remove : bool, optional
        Delete source ASCII files after converting them; keep False when inspecting backend output. Default: False.

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
    See the edcmp tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution. On Windows call under an if __name__ == "__main__" guard.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    tasks = []
    for grp in group_list_edcmp:
        for d in grp:
            tasks.append(tuple(d + [path_green, check_finished]))

    processes = None
    try:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            processes = json.load(fr).get("processes_num", None)
    except Exception:
        pass

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_edcmp2_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Computing static lib",
        ):
            pass
    if convert_bulk:
        convert_pd2bin_edcmp2_all(path_green, remove=remove)
    e = datetime.datetime.now()
    return e - s


def create_grnlib_edcmp2_parallel_multi_nodes(path_green, check_finished=False):
    """Compute the prepared edcmp library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

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
    See the edcmp tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    for ind_group in range(len(group_list_edcmp)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num != len(group_list_edcmp[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!"
                % (processes_num, len(group_list_edcmp[0]))
            )
        print("ind_group:%d rank:%d" % (ind_group, rank))
        # the last group holds the remainder and may be shorter than processes_num
        if rank >= len(group_list_edcmp[ind_group]):
            continue
        call_edcmp2(
            event_depth=group_list_edcmp[ind_group][rank][0],
            obs_depth=group_list_edcmp[ind_group][rank][1],
            mt_ind=group_list_edcmp[ind_group][rank][2],
            path_green=path_green,
            check_finished=check_finished,
        )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


_EDCMP_OUTPUT_TYPE_NAMES = ["disp", "strain", "stress", "tilt"]
_EDCMP_CHA_NUM = {"disp": 3, "strain": 6, "stress": 6, "tilt": 2}


def convert_pd2bin_edcmp2_all(path_green, remove=False):
    """Convert all completed edcmp outputs to float32 binary libraries.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    remove : bool, optional
        Delete source ASCII files after converting them; keep False when inspecting backend output. Default: False.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units.
    """
    print("converting ascii files to binary float32 files")
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    grn_source_depth_range = green_info["grn_source_depth_range"]
    grn_source_delta_depth = green_info["grn_source_delta_depth"]
    event_depth_list = cal_grid(
        grn_source_depth_range[0],
        grn_source_depth_range[1],
        grn_source_delta_depth,
    )
    obs_depth_list = green_info["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]

    grn_dist_range = green_info["grn_dist_range"]
    grn_dist_delta = green_info["grn_delta_dist"]
    n_dist = len(
        cal_grid(grn_dist_range[0], grn_dist_range[1], grn_dist_delta)
    )

    output_observables = np.nonzero(np.array(green_info["output_observables"]))[0]
    n_dep = len(event_depth_list)
    n_obs = len(obs_depth_list)

    # Pre-allocate one bulk array per active output_type:
    # shape (n_dep, n_obs, 5, n_dist, cha_num)
    bulk_data = {}
    for o in output_observables:
        ot = _EDCMP_OUTPUT_TYPE_NAMES[int(o)]
        bulk_data[ot] = np.zeros(
            (n_dep, n_obs, 5, n_dist, _EDCMP_CHA_NUM[ot]), dtype=np.float32
        )

    for i, event_depth in enumerate(event_depth_list):
        for j, obs_depth in enumerate(obs_depth_list):
            for k in range(5):
                path_sub_dir = os.path.join(
                    path_green,
                    "edcmp2",
                    "%.2f" % event_depth,
                    "%.2f" % obs_depth,
                    "%d" % k,
                )
                for o in output_observables:
                    v_ijko = convert_edcmp2(
                        path_sub_dir=path_sub_dir,
                        output_type_ind=int(o),
                        remove=remove,
                    )
                    bulk_data[_EDCMP_OUTPUT_TYPE_NAMES[int(o)]][i, j, k] = v_ijko

    for ot, arr in bulk_data.items():
        out_path = os.path.join(path_green, "edcmp2_%s.bin" % ot)
        arr.tofile(out_path)
        print("Saved %s, shape: %s" % (os.path.basename(out_path), str(arr.shape)))
