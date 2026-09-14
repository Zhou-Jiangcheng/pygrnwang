import os
import json
import pickle
import datetime
from multiprocessing import Pool

from tqdm import tqdm

from .create_spgrn2020 import create_dir_spgrn, create_inp_spgrn2020, call_spgrn2020
from .read_green_info_spgrn import read_green_info_spgrn
from .utils import group, convert_earth_model_nd2nd_without_Q


def _call_spgrn2020_star(args):
    return call_spgrn2020(*args)


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_spgrn2020(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    source_radius,
    cal_gf,
    time_window,
    green_before_p,
    source_duration,
    dist_range,
    delta_dist_range,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    """Prepare the spgrn2020 library grid, input files and job groups.

    Parameters
    ----------
    processes_num : int
        Positive worker count used to group jobs; MPI rank count must match the prepared group width.
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    event_depth_list : list of float
        Source depth nodes in km, positive down; supply a nonempty sorted list.
    receiver_depth_list : list of float
        Receiver depth nodes in km, positive down; supply a nonempty sorted list.
    spec_time_window : float
        Spectral calculation duration in seconds; it must cover the requested output window.
    sampling_interval : float
        Time step in seconds; choose it consistently with the highest modeled frequency.
    max_frequency : float
        Highest modeled frequency in Hz; must be compatible with the time step.
    max_slowness : float
        Maximum modeled slowness in s/km; in SPGRN, nonpositive requests the complete wavefield.
    anti_alias : float
        Dimensionless time-domain alias suppression factor; use a small positive value below 1.
    gravity_fc : float
        Critical frequency in Hz below which self-gravity is included together with the harmonic cutoff.
    gravity_harmonic : int
        Critical spherical harmonic degree for self-gravity.
    cal_sph : int
        1 enables spheroidal (P-SV) modes; 0 disables them.
    cal_tor : int
        1 enables toroidal (SH) modes; 0 disables them.
    source_radius : float
        Source patch radius in km.
    cal_gf : int
        1 updates source spectra; 0 reuses spectra with identical model and spectral parameters.
    time_window : float
        Output time-window duration in seconds.
    green_before_p : float
        Positive seconds before direct P at the start of an SPGRN2020 trace.
    source_duration : float
        Squared half-sinusoid source-time-function duration in seconds; zero requests the backend minimum.
    dist_range : list of float
        Minimum and maximum epicentral distances in km.
    delta_dist_range : list of float
        Smallest and largest distance increments in km at the near and far limits; equal values produce a regular grid.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.
    physical_dispersion : int, optional
        0 disables, 1 enables the backend physical-dispersion correction associated with attenuation. Default: 0.

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
    See the spgrn2020 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions.
    """
    item_list = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            create_dir_spgrn(event_depth, receiver_depth, path_green)
            create_inp_spgrn2020(
                path_green,
                event_depth,
                receiver_depth,
                spec_time_window,
                sampling_interval,
                max_frequency,
                max_slowness,
                anti_alias,
                gravity_fc,
                gravity_harmonic,
                cal_sph,
                cal_tor,
                source_radius,
                cal_gf,
                time_window,
                green_before_p,
                source_duration,
                dist_range,
                delta_dist_range,
                path_nd,
                earth_model_layer_num,
                physical_dispersion,
            )
            item_list.append([event_depth, receiver_depth])
    group_list = group(item_list, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    params = {
        "processes_num": processes_num,
        "path_green": path_green,
        "event_depth_list": event_depth_list,
        "receiver_depth_list": receiver_depth_list,
        "spec_time_window": spec_time_window,
        "sampling_interval": sampling_interval,
        "max_frequency": max_frequency,
        "max_slowness": max_slowness,
        "anti_alias": anti_alias,
        "gravity_fc": gravity_fc,
        "gravity_harmonic": gravity_harmonic,
        "cal_sph": cal_sph,
        "cal_tor": cal_tor,
        "source_radius": source_radius,
        "cal_gf": cal_gf,
        "time_window": time_window,
        # "sampling_num": round(time_window / sampling_interval + 1),
        "green_before_p": green_before_p,
        "source_duration": source_duration,
        "dist_range": dist_range,
        "delta_dist_range": delta_dist_range,
        "path_nd": path_nd,
        "earth_model_layer_num": earth_model_layer_num,
        "physical_dispersion": physical_dispersion,
        "path_nd_without_Q": path_nd_without_Q,
    }
    json_str = json.dumps(params, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    return group_list


def update_green_info_lib_json(path_green, event_depth, receiver_depth):
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    green_info_dep = read_green_info_spgrn(
        str(
            os.path.join(
                path_green, "GreenFunc", "%.2f" % event_depth, "%.2f" % receiver_depth
            )
        ),
        event_depth,
    )
    green_info["dist_list"] = green_info_dep["dist_list"]
    green_info["samples_num"] = green_info_dep["samples_num"]
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)
    return green_info


def create_grnlib_spgrn2020_sequential(path_green, check_finished=False):
    """Compute the prepared spgrn2020 library sequentially.

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
    See the spgrn2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in group_list:
        for i in range(len(item)):
            print("computing " + str(item[i]))
            call_spgrn2020(item[i][0], item[i][1], path_green, check_finished)
    update_green_info_lib_json(path_green, group_list[0][0][0], group_list[0][0][1])
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))


def create_grnlib_spgrn2020_parallel(path_green, check_finished=False):
    """Compute the prepared spgrn2020 library with local worker processes.

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
    See the spgrn2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution. On Windows call under an if __name__ == "__main__" guard.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)

    tasks = []
    for grp in group_list:
        for item in grp:
            tasks.append(tuple(item + [path_green, check_finished]))

    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    processes = green_info.get("processes_num", None)

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_spgrn2020_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Computing SPGRN2020 library",
        ):
            pass

    update_green_info_lib_json(
        path_green,
        float(green_info["event_depth_list"][0]),
        float(green_info["receiver_depth_list"][0]),
    )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_spgrn2020_parallel_multi_nodes(path_green, check_finished=False):
    """Compute the prepared spgrn2020 library with MPI.

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
    See the spgrn2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for ind_group in range(len(group_list)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num != len(group_list[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!" % (processes_num, len(group_list[0]))
            )
        print("ind_group:%d rank:%d" % (ind_group, rank))
        # the last group holds the remainder and may be shorter than processes_num
        if rank >= len(group_list[ind_group]):
            continue
        call_spgrn2020(
            event_depth=group_list[ind_group][rank][0],
            receiver_depth=group_list[ind_group][rank][1],
            path_green=path_green,
            check_finished=check_finished,
        )
    update_green_info_lib_json(path_green, group_list[0][0][0], group_list[0][0][1])
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


if __name__ == "__main__":
    pass
