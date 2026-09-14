import os
import pickle
import json
import datetime
from multiprocessing import Pool

from tqdm import tqdm

from .create_qseis06 import (
    create_dir_qseis06,
    create_inp_qseis06,
    create_inp_qseis06_points,
    call_qseis06,
    convert_pd2bin_qseis06,
)
from .pytaup import taup_create_npz_file, create_tpts_table
from .utils import group, convert_earth_model_nd2nd_without_Q, cal_grid


# 新增：imap_unordered 的打包调用助手
def _call_qseis06_star(args):
    return call_qseis06(*args)


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def create_order_ind(order, diff_accu_order):
    if order == diff_accu_order // 2:
        order_ind = 0
    elif order < diff_accu_order // 2:
        order_ind = order + 1
    else:
        order_ind = order
    return order_ind


def pre_process_qseis06(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    dist_range,
    delta_dist,
    N_each_group,
    time_window,
    sampling_interval,
    slowness_int_algorithm=0,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=12,
    anti_alias=0.01,
    free_surface=True,
    wavelet_duration=0,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
    check_finished_tpts_table=False,
):
    """Prepare the qseis06 library grid, input files and job groups.

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
    slowness_int_algorithm : int, optional
        QSEIS integration selector: 0 for the full wavefield; 1 or 2 for narrow tapered slowness windows. Default: 0.
    slowness_window : list of float or None, optional
        Four ordered slowness taper corners in s/km; None writes zeros for backend automatic limits. Default: None.
    time_reduction_velo : float, optional
        Reduction velocity in km/s; nonzero starts each trace at distance/velocity seconds, while zero disables reduction. Default: 0.
    wavenumber_sampling_rate : float, optional
        Dimensionless spatial Nyquist oversampling factor for wavenumber integration. Default: 12.
    anti_alias : float, optional
        Dimensionless time-domain alias suppression factor; use a small positive value below 1. Default: 0.01.
    free_surface : bool or int, optional
        Backend free-surface selection; see Notes for the backend-specific encoding. Default: True.
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
    See the qseis06 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. free_surface=True retains the free surface; False filters its effects.
    """
    print("preprocessing")
    os.makedirs(path_green, exist_ok=True)

    N_dist, N_dist_group = None, None
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            # print("creating green func dir, info, inp for event_depth=%.2f receiver_depth=%.2f" %
            #       (event_depth, receiver_depth))
            N_dist, N_dist_group = create_dir_qseis06(
                path_green,
                event_depth,
                receiver_depth,
                dist_range,
                delta_dist,
                N_each_group,
                0,
            )
            path_sub_dir = str(
                os.path.join(path_green, "%.2f" % event_depth, "%.2f" % receiver_depth)
            )
            create_inp_qseis06(
                path_sub_dir=path_sub_dir,
                event_depth=event_depth,
                receiver_depth=receiver_depth,
                dist_range=dist_range,
                delta_dist=delta_dist,
                N_dist=N_dist,
                N_dist_group=N_dist_group,
                N_each_group=N_each_group,
                time_window=time_window,
                sampling_interval=sampling_interval,
                slowness_int_algorithm=slowness_int_algorithm,
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
                order=0,
            )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    # creating tp and ts tables
    npz_file = taup_create_npz_file(nd_file=path_nd_without_Q)
    dist_kms = cal_grid(dist_range[0], dist_range[1], delta_dist)
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                path_green,
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
                inp_list.append([event_dep, receiver_dep, nn, 0])
    group_list = group(inp_list, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore
    return group_list


def pre_process_qseis06_strain_rate(
    processes_num,
    path_green,
    path_bin,
    event_depth_list,
    receiver_depth_list,
    dist_range,
    delta_dist,
    N_each_group,
    time_window,
    sampling_interval,
    slowness_int_algorithm=0,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=12,
    anti_alias=0.01,
    free_surface=True,
    wavelet_duration=0,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
    k_dr=0.001,
    dz=0.1,  # km
    diff_accu_order=4,
    check_finished_tpts_table=False,
):
    """Prepare the qseis06 library grid, input files and job groups.

    Parameters
    ----------
    processes_num : int
        Positive worker count used to group jobs; MPI rank count must match the prepared group width.
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    path_bin : str
        Legacy executable-path metadata for the finite-difference library; execution is handled by the installed QSEIS wrapper.
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
    slowness_int_algorithm : int, optional
        QSEIS integration selector: 0 for the full wavefield; 1 or 2 for narrow tapered slowness windows. Default: 0.
    slowness_window : list of float or None, optional
        Four ordered slowness taper corners in s/km; None writes zeros for backend automatic limits. Default: None.
    time_reduction_velo : float, optional
        Reduction velocity in km/s; nonzero starts each trace at distance/velocity seconds, while zero disables reduction. Default: 0.
    wavenumber_sampling_rate : float, optional
        Dimensionless spatial Nyquist oversampling factor for wavenumber integration. Default: 12.
    anti_alias : float, optional
        Dimensionless time-domain alias suppression factor; use a small positive value below 1. Default: 0.01.
    free_surface : bool or int, optional
        Backend free-surface selection; see Notes for the backend-specific encoding. Default: True.
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
    k_dr : float, optional
        Dimensionless relative radial step: dr = distance * k_dr. Avoid zero epicentral distance. Default: 0.001.
    dz : float, optional
        Receiver-depth finite-difference increment in km. Default: 0.1.
    diff_accu_order : int, optional
        Central-difference accuracy order; one of 2, 4, 6 or 8. Default: 4.
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
    See the qseis06 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. free_surface=True retains the free surface; False filters its effects.
    """
    os.makedirs(path_green, exist_ok=True)
    if diff_accu_order not in [2, 4, 6, 8]:
        raise ValueError("diff_accu_order must be in [2,4,6,8]")

    N_dist, N_dist_group = None, None
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            print(
                "creating green func dir, info, inp for "
                "event_depth=%.2f receiver_depth=%.2f" % (event_depth, receiver_depth)
            )
            N_dist, N_dist_group = create_dir_qseis06(
                path_green,
                event_depth,
                receiver_depth,
                dist_range,
                delta_dist,
                N_each_group,
                diff_accu_order,
            )
            points = cal_grid(dist_range[0], dist_range[1], delta_dist)
            for n_group in range(N_dist_group):
                for order in range(diff_accu_order + 1):
                    points_n_o = points[
                        n_group * N_each_group : (n_group + 1) * N_each_group
                    ]
                    points_n_o = (
                        points_n_o + (order - diff_accu_order // 2) * points_n_o * k_dr
                    )
                    order_ind = create_order_ind(order, diff_accu_order)
                    # print(n_group, order, order_ind, dist_range)
                    create_inp_qseis06_points(
                        path_green=path_green,
                        event_depth=event_depth,
                        receiver_depth=receiver_depth,
                        n_group=n_group,
                        points=points_n_o,
                        time_window=time_window,
                        sampling_interval=sampling_interval,
                        slowness_int_algorithm=slowness_int_algorithm,
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
                        order=order_ind,
                    )
                if receiver_depth > 0:
                    for order in range(diff_accu_order + 1):
                        if order == diff_accu_order // 2:
                            continue
                        receiver_depth_inp = (
                            receiver_depth + (order - diff_accu_order // 2) * dz
                        )
                        order_ind = (
                            create_order_ind(order, diff_accu_order) + diff_accu_order
                        )
                        path_sub_dir = str(
                            os.path.join(
                                path_green,
                                "%.2f" % event_depth,
                                "%.2f" % receiver_depth,
                            )
                        )
                        create_inp_qseis06(
                            path_sub_dir=path_sub_dir,
                            event_depth=event_depth,
                            receiver_depth=receiver_depth_inp,
                            dist_range=dist_range,
                            delta_dist=delta_dist,
                            N_dist=N_dist,
                            N_dist_group=N_dist_group,
                            N_each_group=N_each_group,
                            time_window=time_window,
                            sampling_interval=sampling_interval,
                            slowness_int_algorithm=slowness_int_algorithm,
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
                            order=order_ind,
                        )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    if path_nd is not None:
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
        "dist_range": dist_range,
        "delta_dist": delta_dist,
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
        "k_dr": k_dr,
        "dz": dz,
        "diff_accu_order": diff_accu_order,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    inp_list = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            for n_group in range(N_dist_group):
                if receiver_depth > 0:
                    for order in range(2 * diff_accu_order + 1):
                        inp_list.append([event_depth, receiver_depth, n_group, order])
                else:
                    for order in range(diff_accu_order + 1):
                        inp_list.append([event_depth, receiver_depth, n_group, order])
    inp_list_sorted = sorted(inp_list, key=lambda x: abs(x[0] - x[1]))
    group_list = group(inp_list_sorted, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore
    return group_list


def create_grnlib_qseis06_sequential(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qseis06 library sequentially.

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
    See the qseis06 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
    """
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in tqdm(group_list, desc="Computing Green's Func Lib"):
        for i in range(len(item)):
            # print("computing " + str(item[i]))
            item[i] = item[i] + [path_green, check_finished]
            call_qseis06(*item[i])
    if convert_pd2bin:
        convert_pd2bin_qseis06_all(path_green, remove_pd)


def create_grnlib_qseis06_parallel(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    """Compute the prepared qseis06 library with local worker processes.

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
    See the qseis06 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution. On Windows call under an if __name__ == "__main__" guard.
    """
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    # 展平任务
    tasks = []
    for grp in group_list:
        for item in grp:
            tasks.append(tuple(item + [path_green, check_finished]))

    # 尝试读取进程数
    processes = None
    try:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            processes = json.load(fr).get("processes_num", None)
    except Exception:
        pass

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_qseis06_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Computing Green's Func Lib",
        ):
            pass

    if convert_pd2bin:
        convert_pd2bin_qseis06_all(path_green, remove_pd)


def create_grnlib_qseis06_parallel_multi_nodes(path_green, check_finished=False):
    """Compute the prepared qseis06 library with MPI.

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
    See the qseis06 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Check output files and logs after execution.
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
        call_qseis06(
            event_depth=group_list[ind_group][rank][0],
            receiver_depth=group_list[ind_group][rank][1],
            n_group=group_list[ind_group][rank][2],
            order=group_list[ind_group][rank][3],
            path_green=path_green,
            check_finished=check_finished,
        )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def convert_pd2bin_qseis06_all(path_green, remove=False):
    """Convert all completed qseis06 outputs to float32 binary libraries.

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
    See the qseis06 tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units.
    """
    print("converting ascii files to bytes files")
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
                    convert_pd2bin_qseis06(
                        os.path.join(sub_dir, sub_sub_dir), remove=remove
                    )


if __name__ == "__main__":
    pass
