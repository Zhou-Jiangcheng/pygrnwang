import os
import shutil
import platform
import pickle
import json
import datetime
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm
from mpi4py import MPI
import jpype

from .create_qseis import (
    create_dir_qseis06,
    create_inp_qseis06,
    create_inp_qseis06_points,
    call_qseis06,
    convert_pd2bin_qseis06,
)
from .pytaup import create_tpts_table
from .utils import group, convert_earth_model_nd2nd_without_Q


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
    check_finished_tpts_table=False,
):
    print("preprocessing")
    os.makedirs(path_green, exist_ok=True)
    if platform.system() == "Windows":
        path_bin_call = os.path.join(path_green, "qseis06.exe")
    else:
        path_bin_call = os.path.join(path_green, "qseis06.bin")
    shutil.copy(path_bin, path_bin_call)

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
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                path_green,
                event_depth,
                receiver_depth,
                dist_range,
                delta_dist,
                path_nd_without_Q,
                check_finished_tpts_table,
            )
    if jpype.isJVMStarted():
        jpype.shutdownJVM()

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
    os.makedirs(path_green, exist_ok=True)
    if diff_accu_order not in [2, 4, 6, 8]:
        raise ValueError("diff_accu_order must be in [2,4,6,8]")
    if platform.system() == "Windows":
        path_bin_call = os.path.join(path_green, "qseis06.exe")
    else:
        path_bin_call = os.path.join(path_green, "qseis06.bin")
    shutil.copy(path_bin, path_bin_call)

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
            points = np.linspace(dist_range[0], dist_range[1], N_dist)
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
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                path_green,
                event_depth,
                receiver_depth,
                dist_range,
                delta_dist,
                path_nd_without_Q,
                check_finished_tpts_table,
            )
    if jpype.isJVMStarted():
        jpype.shutdownJVM()

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
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in tqdm(group_list, desc="Computing Green's Func Lib"):
        for i in range(len(item)):
            # print("computing " + str(item[i]))
            item[i] = item[i] + [path_green, check_finished]
            call_qseis06(*item[i])
    if convert_pd2bin:
        convert_pd2bin_qseis06_all(path_green, remove_pd)


def create_grnlib_qseis06_parallel_single_node(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in tqdm(group_list, desc="Computing Green's Func Lib"):
        # print("computing " + str(item))
        for i in range(len(item)):
            item[i] = item[i] + [path_green, check_finished]
        pool = Pool()
        r = pool.starmap_async(call_qseis06, item)
        r.get()
        pool.close()
        pool.join()
    if convert_pd2bin:
        convert_pd2bin_qseis06_all(path_green, remove_pd)


def create_grnlib_qseis06_parallel_multi_nodes(path_green, check_finished=False):
    s = datetime.datetime.now()
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
    ind_group = rank // processes_num
    ind_para = rank - ind_group * processes_num
    print(rank, ind_group, ind_para)
    call_qseis06(
        event_depth=group_list[ind_group][ind_para][0],
        receiver_depth=group_list[ind_group][ind_para][1],
        n_group=group_list[ind_group][ind_para][2],
        order=group_list[ind_group][ind_para][3],
        path_green=path_green,
        check_finished=check_finished,
    )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def convert_pd2bin_qseis06_all(path_green, remove=False):
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
