import os
import pickle
import json
import datetime
from multiprocessing import Pool

from tqdm import tqdm

from .create_qseis2025 import (
    create_dir_qseis2025,
    create_inp_qseis2025,
    call_qseis2025,
    convert_pd2bin_qseis2025,
)
from .pytaup import create_tpts_table
from .utils import group, convert_earth_model_nd2nd_without_Q, cal_grid


def _call_qseis2025_star(args):
    return call_qseis2025(*args)


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_qseis2025(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    dist_range,
    delta_dist,
    N_each_group,
    time_window,
    sampling_interval,
    output_observables,
    slowness_int_algorithm=0,
    eps_estimate_wavenumber=1e-6,
    source_radius_ratio=0.05,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=12,
    anti_alias=0.01,
    free_surface=0,
    wavelet_duration=0,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
    check_finished_tpts_table=False,
):
    """Prepare the qseis2025 library grid, input files and job groups.

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
    dist_range : list of float
        Minimum and maximum epicentral distances in km.
    delta_dist : float
        Positive regular distance increment in km. The last grid point can exceed the requested maximum by less than one increment.
    N_each_group : int
        Positive maximum number of distances in each backend input file.
    time_window : float
        Output time-window duration in seconds.
    sampling_interval : float
        Time step in seconds; choose it consistently with the highest modeled frequency.
    output_observables : list of int
        Five 0/1 flags in displacement (or velocity), volume, strain, stress, rotation order; rate versus integrated quantities depend on wavelet_type.
    slowness_int_algorithm : int, optional
        QSEIS integration selector: 0 for the full wavefield; 1 or 2 for narrow tapered slowness windows. Default: 0.
    eps_estimate_wavenumber : float, optional
        Full-wavefield wavenumber truncation tolerance; smaller values increase accuracy and computation. Default: 1e-06.
    source_radius_ratio : float, optional
        Dimensionless Gaussian spatial-smoothing ratio. At each frequency and
        receiver, the native radius is ``source_radius_ratio *
        min(sqrt(r**2 + (zs-zr)**2), Vp_source/(f+df))``, using epicentral
        distance r, source/receiver depths zs/zr, source-layer P-wave speed,
        frequency f and FFT increment df in consistent units. Coordinates and
        speed include any selected Earth flattening. The kernel is multiplied
        by ``exp(-(k*radius)**2/2)`` at wavenumber k; this is not a fixed-radius
        source disk. Zero disables smoothing. Larger values generally reduce
        the automatically estimated wavenumber cutoff and computation time.
        Default: 0.05.
    slowness_window : list of float or None, optional
        Four ordered slowness taper corners in s/km; None writes zeros for backend automatic limits. Default: None.
    time_reduction_velo : float, optional
        Reduction velocity in km/s; nonzero starts each trace at distance/velocity seconds, while zero disables reduction. Default: 0.
    wavenumber_sampling_rate : float, optional
        Dimensionless spatial Nyquist oversampling factor for wavenumber integration. Default: 12.
    anti_alias : float, optional
        Dimensionless time-domain alias suppression factor; use a small positive value below 1. Default: 0.01.
    free_surface : bool or int, optional
        Backend free-surface selection; see Notes for the backend-specific encoding. Default: 0.
    wavelet_duration : int, optional
        Wavelet duration in native time samples, not seconds. Nonpositive values request the backend default of two samples. Default: 0.
    wavelet_type : int, optional
        1 selects the normalized squared half-sinusoid; 2 selects its tapered Heaviside integral. Custom type 0 requires manually supplying wavelet samples. Default: 1.
    flat_earth_transform : bool, optional
        Apply the backend flat-Earth transformation and receiver-radius distance correction. Default: True.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.
    check_finished_tpts_table : bool, optional
        Reuse existing P/S table files without validating their model or grid provenance. Default: False.

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
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. free_surface=0 includes the free surface; 1 removes it; 2 removes it with amplitude correction for surface receivers.
    """
    print("Preprocessing")
    os.makedirs(path_green, exist_ok=True)

    N_dist, N_dist_group = None, None
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            # print("creating green func dir, info, inp for event_depth=%.2f receiver_depth=%.2f" %
            #       (event_depth, receiver_depth))
            N_dist, N_dist_group = create_dir_qseis2025(
                path_green=path_green,
                event_depth=event_depth,
                receiver_depth=receiver_depth,
                dist_range=dist_range,
                delta_dist=delta_dist,
                N_each_group=N_each_group,
            )
            create_inp_qseis2025(
                path_green=path_green,
                event_depth=event_depth,
                receiver_depth=receiver_depth,
                dist_range=dist_range,
                delta_dist=delta_dist,
                N_dist=N_dist,
                N_dist_group=N_dist_group,
                N_each_group=N_each_group,
                time_window=time_window,
                sampling_interval=sampling_interval,
                output_observables=output_observables,
                slowness_int_algorithm=slowness_int_algorithm,
                eps_estimate_wavenumber=eps_estimate_wavenumber,
                source_radius_ratio=source_radius_ratio,
                slowness_window=slowness_window,
                time_reduction_velo=time_reduction_velo,
                wavenumber_sampling_rate=wavenumber_sampling_rate,
                anti_alias=anti_alias,
                free_surface=free_surface,
                wavelet_duration=wavelet_duration,
                wavelet_type=wavelet_type,
                flat_earth_transform=flat_earth_transform,
                path_nd=path_nd,
                earth_model_layer_num=earth_model_layer_num,
            )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    # creating tp and ts tables
    dist_kms = cal_grid(dist_range[0], dist_range[1], delta_dist)
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                path_green,
                event_depth,
                receiver_depth,
                dist_kms,
                path_nd_without_Q,
                check_finished_tpts_table,
            )

    green_info = {
        "processes_num": processes_num,
        "event_depth_list": event_depth_list,
        "receiver_depth_list": receiver_depth_list,
        "grn_dist_range": dist_range,
        "grn_delta_dist": delta_dist,
        "N_dist": N_dist,
        "N_dist_group": N_dist_group,
        "N_each_group": N_each_group,
        "time_window": time_window,
        "sampling_interval": sampling_interval,
        "sampling_num": round(time_window / sampling_interval + 1),
        "slowness_int_algorithm": slowness_int_algorithm,
        "slowness_window": slowness_window,
        "time_reduction_velo": time_reduction_velo,
        "wavenumber_sampling_rate": wavenumber_sampling_rate,
        "anti_alias": anti_alias,
        "free_surface": free_surface,
        "wavelet_duration": wavelet_duration,
        "wavelet_type": wavelet_type,
        "flat_earth_transform": flat_earth_transform,
        "path_nd": path_nd,
        "path_nd_without_Q": path_nd_without_Q,
        "earth_model_layer_num": earth_model_layer_num,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    inp_list = []
    for event_dep in event_depth_list:
        for receiver_dep in receiver_depth_list:
            for nn in range(N_dist_group):
                inp_list.append([event_dep, receiver_dep, nn])
    # inp_list_sorted = sorted(inp_list, key=lambda x: abs(x[0] - x[1]))
    group_list = group(inp_list, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore
    return group_list


def create_grnlib_qseis2025_sequential(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qseis2025 library sequentially.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
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
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in tqdm(group_list, desc="Computing dynamic stress"):
        for i in range(len(item)):
            # print("computing " + str(item[i]))
            item[i] = item[i] + [path_green, check_finished]
            call_qseis2025(*item[i])
    if convert_pd2bin:
        convert_pd2bin_qseis2025_all(path_green, remove_pd)


def create_grnlib_qseis2025_parallel(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qseis2025 library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
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
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution. On Windows call under an if __name__ == "__main__" guard.
    """
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    tasks = []
    for grp in group_list:
        for item in grp:
            tasks.append(tuple(item + [path_green, check_finished]))

    processes = None
    try:
        with open(
            os.path.join(path_green, "green_lib_info.json"), "r", encoding="utf-8"
        ) as fr:
            processes = json.load(fr).get("processes_num", None)
    except Exception:
        pass

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_qseis2025_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Compute QSEIS2025 Green's library",
        ):
            pass

    if convert_pd2bin:
        convert_pd2bin_qseis2025_all(path_green, remove_pd)


def convert_pd2bin_qseis2025_all(path_green, remove=False):
    """Convert all completed qseis2025 outputs to float32 binary libraries.

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
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units.
    """
    print("Converting ascii files to bytes files")
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    event_depth_list = green_info["event_depth_list"]
    receiver_depth_list = green_info["receiver_depth_list"]
    for event_dep in event_depth_list:
        for receiver_dep in receiver_depth_list:
            sub_dir = str(
                os.path.join(path_green, "%.2f" % event_dep, "%.2f" % receiver_dep)
            )
            sub_sub_dirs = os.listdir(sub_dir)
            for sub_sub_dir in sub_sub_dirs:
                if "_table.bin" not in sub_sub_dir:
                    convert_pd2bin_qseis2025(
                        os.path.join(sub_dir, sub_sub_dir), remove=remove
                    )


def create_grnlib_qseis2025_parallel_multi_nodes(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qseis2025 library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
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
    RuntimeError
        mpi4py is unavailable.
    ValueError
        MPI rank count does not match the prepared group width.

    Notes
    -----
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    comm = MPI.COMM_WORLD
    processes_num = comm.Get_size()
    if processes_num != len(group_list[0]):
        raise ValueError(
            "processes_num is %d, item num in group is %d. \n"
            "Pleasse check the process num!" % (processes_num, len(group_list[0]))
        )
    rank = comm.Get_rank()
    for ind_group in range(len(group_list)):
        # the last group holds the remainder and may be shorter than processes_num
        if rank >= len(group_list[ind_group]):
            continue
        print("ind_group:%d rank:%d" % (ind_group, rank))
        call_qseis2025(
            event_depth=group_list[ind_group][rank][0],
            receiver_depth=group_list[ind_group][rank][1],
            n_group=group_list[ind_group][rank][2],
            path_green=path_green,
            check_finished=check_finished,
        )
    if convert_pd2bin:
        # every rank writes the same .bin files, let one rank do it after all are done
        comm.Barrier()
        if rank == 0:
            convert_pd2bin_qseis2025_all(path_green, remove_pd)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


if __name__ == "__main__":
    pass
