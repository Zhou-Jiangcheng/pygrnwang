import os
import glob
import pickle
import shutil
import platform
import json
import datetime
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm
import jpype

try:
    from mpi4py import MPI
except:
    pass

from .create_qssp2020 import (
    mt_com_list,
    create_inp_qssp2020,
    create_dir_qssp2020,
    call_qssp2020,
    convert_pd2bin_qssp2020,
)
from .utils import group, convert_earth_model_nd2nd_without_Q
from .pytaup import create_tpts_table


# 新增：imap_unordered 的打包调用助手
def _call_qssp2020_star(args):
    return call_qssp2020(*args)


def pre_process_spec(
    processes_num,
    path_green,
    path_bin,
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
    path_bin,
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
    path_bin,
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
    print("Preprocessing")
    pre_process_spec(
        processes_num,
        path_green,
        path_bin,
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
        path_bin,
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
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                os.path.join(path_green, "GreenFunc"),
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
    s = datetime.datetime.now()
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
    s = datetime.datetime.now()
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
