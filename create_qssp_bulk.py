import os
import shutil
import datetime
import pickle
from multiprocessing import Pool

from mpi4py import MPI

from pygrnwang.create_qssp import mt_com_list, create_inp, create_dir, call_qssp2020
from pygrnwang.utils import group


def pre_process_spec(
    processes_num,
    event_depth_list,
    receiver_depth_list,
    path_green,
    path_bin,
    spec_time_window,
    time_window,
    time_reduction,
    dist_range,
    delta_dist_range,
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
    path_nd=None,
    earth_model_layer_num=None,
):
    path_bin_call = os.path.join(path_green, "qssp2020.bin")
    if not os.path.exists(path_bin_call):
        shutil.copy(path_bin, path_bin_call)
    item_list_spec = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            create_dir(event_depth, receiver_depth, path_green)
            create_inp(
                "spec",
                event_depth,
                receiver_depth,
                path_green,
                spec_time_window,
                time_window,
                time_reduction,
                dist_range,
                delta_dist_range,
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
                1,
                [0 for _ in range(11)],
                path_nd,
                earth_model_layer_num,
            )
            item_list_spec.append([event_depth, receiver_depth, "spec"])

    group_list_spec = group(item_list_spec, processes_num)
    with open(os.path.join(path_green, "group_list_spec.pkl"), "wb") as fw:
        pickle.dump(group_list_spec, fw)
    return group_list_spec


def pre_process_func(
    processes_num,
    event_depth_list,
    receiver_depth_list,
    path_green,
    path_bin,
    spec_time_window,
    time_window,
    time_reduction,
    dist_range,
    delta_dist_range,
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
    output_observables: list,
    path_nd=None,
    earth_model_layer_num=None,
):
    path_bin_call = os.path.join(path_green, "qssp2020.bin")
    if not os.path.exists(path_bin_call):
        shutil.copy(path_bin, path_bin_call)
    item_list_func = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            for mt_com in mt_com_list:
                create_inp(
                    mt_com,
                    event_depth,
                    receiver_depth,
                    path_green,
                    spec_time_window,
                    time_window,
                    time_reduction,
                    dist_range,
                    delta_dist_range,
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
                    0,
                    output_observables,
                    path_nd,
                    earth_model_layer_num,
                )
                item_list_func.append([event_depth, receiver_depth, mt_com])

    group_list_func = group(item_list_func, processes_num)
    with open(os.path.join(path_green, "group_list_func.pkl"), "wb") as fw:
        pickle.dump(group_list_func, fw)
    return group_list_func


def pre_process(
        processes_num,
        event_depth_list,
        receiver_depth_list,
        path_green,
        path_bin,
        spec_time_window,
        time_window,
        time_reduction,
        dist_range,
        delta_dist_range,
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
        output_observables,
        path_nd=None,
        earth_model_layer_num=None,
):
    pre_process_spec(
        processes_num,
        event_depth_list,
        receiver_depth_list,
        path_green,
        path_bin,
        spec_time_window,
        time_window,
        time_reduction,
        dist_range,
        delta_dist_range,
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
        path_nd,
        earth_model_layer_num,
    )
    pre_process_func(
        processes_num,
        event_depth_list,
        receiver_depth_list,
        path_green,
        path_bin,
        spec_time_window,
        time_window,
        time_reduction,
        dist_range,
        delta_dist_range,
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
        output_observables,
        path_nd,
        earth_model_layer_num,
    )


def create_grnlib_parallel_single_node(path_green):
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_spec.pkl"), "rb") as fr:
        group_list_spec = pickle.load(fr)
    for item in group_list_spec:
        # print(
        #     "computing event_dep %.1f km, receiver_dep %.1f km"
        #     % (item[0], item[1])
        # )
        print("computing spec lib " + str(item))
        for i in range(len(item)):
            item[i] = item[i] + [path_green]
        pool = Pool()
        r = pool.starmap_async(call_qssp2020, item)
        r.get()
        pool.close()
        pool.join()

    with open(os.path.join(path_green, "group_list_func.pkl"), "rb") as fr:
        group_list_func = pickle.load(fr)
    for item in group_list_func:
        # print(
        #     "computing event_dep %.1f km, receiver_dep %.1f km, %s"
        #     % (item[0], item[1], item[2])
        # )
        print("computing time lib " + str(item))
        for i in range(len(item)):
            item[i] = item[i] + [path_green]
        pool = Pool()
        r = pool.starmap_async(call_qssp2020, item)
        r.get()
        pool.close()
        pool.join()
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_spec_parallel_multi_nodes(path_green):
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
            "computing spec lib ind_group:%d rank:%d event_depth:%.1f receiver_depth:%.1f"
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
            )


def create_grnlib_func_parallel_multi_nodes(path_green):
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
            "computing time lib ind_group:%d rank:%d event_depth:%.1f receiver_depth:%.1f mt_com:%s"
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
            )


if __name__ == "__main__":
    pass
