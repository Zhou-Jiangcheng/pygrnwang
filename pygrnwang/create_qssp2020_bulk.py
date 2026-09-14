import os
import glob
import pickle
import json
import datetime
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm

from .create_qssp2020 import (
    mt_com_list,
    create_inp_qssp2020,
    create_dir_qssp2020,
    call_qssp2020,
    convert_pd2bin_qssp2020,
)
from .utils import group, convert_earth_model_nd2nd_without_Q, cal_grid
from .pytaup import taup_create_npz_file, create_tpts_table


def _call_qssp2020_star(args):
    return call_qssp2020(*args)


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_spec(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    turning_point_filter,
    turning_point_d1,
    turning_point_d2,
    free_surface_filter,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    min_harmonic,
    max_harmonic,
    source_radius,
    source_duration,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    item_list_spec = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            create_dir_qssp2020(event_depth, receiver_depth, path_green)
            create_inp_qssp2020(
                "spec",
                path_green,
                event_depth,
                receiver_depth,
                spec_time_window,
                sampling_interval,
                max_frequency,
                max_slowness,
                anti_alias,
                turning_point_filter,
                turning_point_d1,
                turning_point_d2,
                free_surface_filter,
                gravity_fc,
                gravity_harmonic,
                cal_sph,
                cal_tor,
                min_harmonic,
                max_harmonic,
                source_radius,
                1,
                source_duration,
                [0 for _ in range(11)],
                time_window,
                time_reduction,
                dist_range,
                delta_dist,
                path_nd,
                earth_model_layer_num,
                physical_dispersion,
            )
            item_list_spec.append([event_depth, receiver_depth, "spec"])

    group_list_spec = group(item_list_spec, processes_num)
    with open(os.path.join(path_green, "group_list_spec.pkl"), "wb") as fw:
        pickle.dump(group_list_spec, fw)  # type: ignore
    return group_list_spec


def pre_process_func(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    turning_point_filter,
    turning_point_d1,
    turning_point_d2,
    free_surface_filter,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    min_harmonic,
    max_harmonic,
    source_radius,
    source_duration,
    output_observables: list,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    item_list_func = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            for mt_com in mt_com_list:
                create_inp_qssp2020(
                    mt_com,
                    path_green,
                    event_depth,
                    receiver_depth,
                    spec_time_window,
                    sampling_interval,
                    max_frequency,
                    max_slowness,
                    anti_alias,
                    turning_point_filter,
                    turning_point_d1,
                    turning_point_d2,
                    free_surface_filter,
                    gravity_fc,
                    gravity_harmonic,
                    cal_sph,
                    cal_tor,
                    min_harmonic,
                    max_harmonic,
                    source_radius,
                    0,
                    source_duration,
                    output_observables,
                    time_window,
                    time_reduction,
                    dist_range,
                    delta_dist,
                    path_nd,
                    earth_model_layer_num,
                    physical_dispersion,
                )
                item_list_func.append([event_depth, receiver_depth, mt_com])

    group_list_func = group(item_list_func, processes_num)
    with open(os.path.join(path_green, "group_list_func.pkl"), "wb") as fw:
        pickle.dump(group_list_func, fw)  # type: ignore
    return group_list_func


def pre_process_qssp2020(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    turning_point_filter,
    turning_point_d1,
    turning_point_d2,
    free_surface_filter,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    min_harmonic,
    max_harmonic,
    source_radius,
    source_duration,
    output_observables: list,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
    check_finished_tpts_table=False,
):
    """Prepare the qssp2020 library grid, input files and job groups.

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
    turning_point_filter : int
        1 selects the QSSP turning-depth filter; 0 disables it.
    turning_point_d1 : float
        Minimum allowed turning depth in km when the turning-point filter is enabled.
    turning_point_d2 : float
        Maximum allowed turning depth in km when the turning-point filter is enabled.
    free_surface_filter : int
        QSSP switch: 1 includes free-surface reflection, 0 removes it.
    gravity_fc : float
        Critical frequency in Hz below which self-gravity is included together with the harmonic cutoff.
    gravity_harmonic : int
        Critical spherical harmonic degree for self-gravity.
    cal_sph : int
        1 enables spheroidal (P-SV) modes; 0 disables them.
    cal_tor : int
        1 enables toroidal (SH) modes; 0 disables them.
    min_harmonic : int
        Control for estimating the low-frequency baseline of the
        frequency-dependent upper harmonic cutoff. This is not the lowest
        retained degree: the spectral sum still starts at degree zero.
        Converge it together with max_harmonic for the requested observables.
    max_harmonic : int
        Upper harmonic cutoff cap. It also influences the spatial
        differential-transform order through the internal maximum degree,
        so changing it can affect synthesis even when the stored spectral
        cutoff is unchanged. Converge it together with min_harmonic.
    source_radius : float
        Source patch radius in km.
    source_duration : float
        Squared half-sinusoid source-time-function duration in seconds; zero requests the backend minimum.
    output_observables : list of int
        Eleven 0/1 flags: displacement, velocity, acceleration, strain, strain rate, stress, stress rate, rotation, rotation rate, gravitation, gravimeter.
    time_window : float
        Output time-window duration in seconds.
    time_reduction : float
        QSSP trace start time in seconds relative to source origin, not a velocity.
    dist_range : list of float
        Minimum and maximum epicentral distances in km.
    delta_dist : float
        Positive regular distance increment in km. The last grid point can exceed the requested maximum by less than one increment.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.
    physical_dispersion : int, optional
        0 disables, 1 enables the backend physical-dispersion correction associated with attenuation. Default: 0.
    check_finished_tpts_table : bool, optional
        Reuse existing P/S table files without validating their model or grid provenance. Default: False.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. A first computation must include both spectral and time-domain stages. See the :doc:`QSSP2020 tutorial </backends/qssp2020>` for harmonic-cutoff convergence and comparison settings.
    """
    print("Preprocessing")
    pre_process_spec(
        processes_num,
        path_green,
        event_depth_list,
        receiver_depth_list,
        spec_time_window,
        sampling_interval,
        max_frequency,
        max_slowness,
        anti_alias,
        turning_point_filter,
        turning_point_d1,
        turning_point_d2,
        free_surface_filter,
        gravity_fc,
        gravity_harmonic,
        cal_sph,
        cal_tor,
        min_harmonic,
        max_harmonic,
        source_radius,
        source_duration,
        time_window,
        time_reduction,
        dist_range,
        delta_dist,
        path_nd,
        earth_model_layer_num,
        physical_dispersion,
    )
    pre_process_func(
        processes_num,
        path_green,
        event_depth_list,
        receiver_depth_list,
        spec_time_window,
        sampling_interval,
        max_frequency,
        max_slowness,
        anti_alias,
        turning_point_filter,
        turning_point_d1,
        turning_point_d2,
        free_surface_filter,
        gravity_fc,
        gravity_harmonic,
        cal_sph,
        cal_tor,
        min_harmonic,
        max_harmonic,
        source_radius,
        source_duration,
        output_observables,
        time_window,
        time_reduction,
        dist_range,
        delta_dist,
        path_nd,
        earth_model_layer_num,
        physical_dispersion,
    )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    # creating tp and ts tables
    npz_file = taup_create_npz_file(nd_file=path_nd_without_Q)
    dist_kms = cal_grid(dist_range[0], dist_range[1], delta_dist)
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                os.path.join(path_green, "GreenFunc"),
                event_depth,
                receiver_depth,
                dist_kms,
                npz_file,
                check_finished_tpts_table,
            )

    green_info = {
        "processes_num": processes_num,
        "event_depth_list": event_depth_list,
        "receiver_depth_list": receiver_depth_list,
        "spec_time_window": spec_time_window,
        "sampling_interval": sampling_interval,
        "max_frequency": max_frequency,
        "max_slowness": max_slowness,
        "anti_alias": anti_alias,
        "turning_point_filter": turning_point_filter,
        "turning_point_d1": turning_point_d1,
        "turning_point_d2": turning_point_d2,
        "free_surface_filter": free_surface_filter,
        "gravity_fc": gravity_fc,
        "gravity_harmonic": gravity_harmonic,
        "cal_sph": cal_sph,
        "cal_tor": cal_tor,
        "min_harmonic": min_harmonic,
        "max_harmonic": max_harmonic,
        "source_radius": source_radius,
        "source_duration": source_duration,
        "output_observables": output_observables,
        "time_window": time_window,
        "sampling_num": round(time_window / sampling_interval) + 1,
        "time_reduction": time_reduction,
        "grn_dist_range": dist_range,
        "grn_delta_dist": delta_dist,
        "path_nd": path_nd,
        "path_nd_without_Q": path_nd_without_Q,
        "earth_model_layer_num": earth_model_layer_num,
        "physical_dispersion": physical_dispersion,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)


def create_grnlib_qssp2020_sequential(
    path_green, cal_spec=True, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qssp2020 library sequentially.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    cal_spec : bool, optional
        Compute spectra before time-domain synthesis. Keep True for a new QSSP library; False requires compatible existing spectra. Default: True.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.
    convert_pd2bin : bool, optional
        Convert completed ASCII waveforms to the compact float32 reader format. Default: True.
    remove_pd : bool, optional
        Delete original ASCII output; retain it while validating a new calculation. Default: True.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    if cal_spec:
        with open(os.path.join(path_green, "group_list_spec.pkl"), "rb") as fr:
            group_list_spec = pickle.load(fr)
        for item in tqdm(
            group_list_spec,
            desc="Compute the Green's function library in the transformed domain.",
        ):
            for i in range(len(item)):
                item[i] = item[i] + [path_green, check_finished]
                call_qssp2020(*item[i])

    with open(os.path.join(path_green, "group_list_func.pkl"), "rb") as fr:
        group_list_func = pickle.load(fr)
    for item in tqdm(
        group_list_func, desc="Compute the Green's function library in the time domain."
    ):
        for i in range(len(item)):
            item[i] = item[i] + [path_green, check_finished]
            call_qssp2020(*item[i])

    if convert_pd2bin:
        convert_pd2bin_qssp2020_all(path_green)
    if remove_pd:
        remove_dat_files(path_green)


def create_grnlib_qssp2020_parallel(
    path_green, cal_spec=True, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qssp2020 library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    cal_spec : bool, optional
        Compute spectra before time-domain synthesis. Keep True for a new QSSP library; False requires compatible existing spectra. Default: True.
    check_finished : bool, optional
        Reuse outputs marked finished. Markers do not verify that inputs are unchanged. Default: False.
    convert_pd2bin : bool, optional
        Convert completed ASCII waveforms to the compact float32 reader format. Default: True.
    remove_pd : bool, optional
        Delete original ASCII output; retain it while validating a new calculation. Default: True.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution. On Windows call under an if __name__ == "__main__" guard.
    """
    tasks = []

    if cal_spec:
        with open(os.path.join(path_green, "group_list_spec.pkl"), "rb") as fr:
            group_list_spec = pickle.load(fr)
        for grp in group_list_spec:
            for item in grp:
                tasks.append(tuple(item + [path_green, check_finished]))

    processes = None
    try:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            processes = json.load(fr).get("processes_num", None)
    except Exception:
        pass

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_qssp2020_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Compute QSSP2020 Green's library in the transformed domain.",
        ):
            pass

    tasks = []
    with open(os.path.join(path_green, "group_list_func.pkl"), "rb") as fr:
        group_list_func = pickle.load(fr)
    for grp in group_list_func:
        for item in grp:
            tasks.append(tuple(item + [path_green, check_finished]))

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_qssp2020_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Compute QSSP2020 Green's function library in the time domain.",
        ):
            pass

    if convert_pd2bin:
        convert_pd2bin_qssp2020_all(path_green)
    if remove_pd:
        remove_dat_files(path_green)


def create_grnlib_qssp2020_spec_parallel_multi_nodes(path_green, check_finished=False):
    """Compute the prepared qssp2020 library with MPI.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    with open(os.path.join(path_green, "group_list_spec.pkl"), "rb") as fr:
        group_list_spec = pickle.load(fr)
    N_all = 0
    for ind_group in range(len(group_list_spec)):
        N_all = N_all + len(group_list_spec[ind_group])
    for ind_group in range(len(group_list_spec)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num < len(group_list_spec[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!"
                % (processes_num, len(group_list_spec[0]))
            )
        print(
            "computing spec lib ind_group:%d rank:%d event_depth:%.2f receiver_depth:%.2f"
            % (
                ind_group,
                rank,
                group_list_spec[ind_group][rank][0],
                group_list_spec[ind_group][rank][1],
            )
        )
        if ind_group * len(group_list_spec[0]) + rank < N_all:
            call_qssp2020(
                event_depth=group_list_spec[ind_group][rank][0],
                receiver_depth=group_list_spec[ind_group][rank][1],
                mt_com=group_list_spec[ind_group][rank][2],
                path_green=path_green,
                check_finished=check_finished,
            )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_qssp2020_func_parallel_multi_nodes(path_green, check_finished=False):
    """Compute the prepared qssp2020 library with MPI.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    with open(os.path.join(path_green, "group_list_func.pkl"), "rb") as fr:
        group_list_spec = pickle.load(fr)
    N_all = 0
    for ind_group in range(len(group_list_spec)):
        N_all = N_all + len(group_list_spec[ind_group])
    for ind_group in range(len(group_list_spec)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num < len(group_list_spec[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!"
                % (processes_num, len(group_list_spec[0]))
            )
        print(
            "computing time lib ind_group:%d rank:%d event_depth:%.2f receiver_depth:%.2f mt_com:%s"
            % (
                ind_group,
                rank,
                group_list_spec[ind_group][rank][0],
                group_list_spec[ind_group][rank][1],
                group_list_spec[ind_group][rank][2],
            )
        )
        if ind_group * len(group_list_spec[0]) + rank < N_all:
            call_qssp2020(
                event_depth=group_list_spec[ind_group][rank][0],
                receiver_depth=group_list_spec[ind_group][rank][1],
                mt_com=group_list_spec[ind_group][rank][2],
                path_green=path_green,
                check_finished=check_finished,
            )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def convert_pd2bin_qssp2020_all(path_green):
    """Convert all completed qssp2020 outputs to float32 binary libraries.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units.
    """
    print("converting ascii files to byte files")
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    event_dep_list = green_info["event_depth_list"]
    receiver_dep_list = green_info["receiver_depth_list"]
    output_observables = np.nonzero(np.array(green_info["output_observables"]))[0]
    for i in range(len(event_dep_list)):
        for j in range(len(receiver_dep_list)):
            for output_type_ind in output_observables:
                convert_pd2bin_qssp2020(
                    path_green,
                    event_dep_list[i],
                    receiver_dep_list[j],
                    int(output_type_ind),
                )


def remove_dat_files(path_green):
    print("removing dat files")
    path_func = os.path.join(path_green, "GreenFunc")
    for root, dirs, files in os.walk(path_func):
        for file in glob.glob(os.path.join(root, "*.dat")):
            os.remove(file)


if __name__ == "__main__":
    pass
