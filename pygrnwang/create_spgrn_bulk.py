import os
import pickle
import shutil
import platform
import json
import datetime

try:
    from mpi4py import MPI
except:
    pass
from multiprocessing import Pool
from .create_spgrn import (
    create_dir_spgrn2020,
    create_inp_spgrn2020,
    call_spgrn2020,
)
from .read_green_info_spgrn2020 import read_green_info_spgrn2020
from .utils import group, convert_earth_model_nd2nd_without_Q


def pre_process_spgrn2020(
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
    if platform.system() == "Windows":
        path_bin_call = os.path.join(path_green, "spgrn2020.exe")
    else:
        path_bin_call = os.path.join(path_green, "spgrn2020.bin")
    if not os.path.exists(path_bin_call):
        shutil.copy(path_bin, path_bin_call)

    item_list = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            create_dir_spgrn2020(event_depth, receiver_depth, path_green)
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
        "path_green": path_green,
        "path_bin": path_bin,
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
        "sampling_num": round(time_window / sampling_interval + 1),
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
    green_info_dep = read_green_info_spgrn2020(
        str(
            os.path.join(
                path_green, "GreenFunc", "%.2f" % event_depth, "%.2f" % receiver_depth
            )
        ),
        event_depth,
    )
    green_info["dist_list"] = green_info_dep["dist_list"]
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)


def create_grnlib_spgrn2020_sequential(path_green, check_finished=False):
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


def create_grnlib_spgrn2020_parallel_single_node(path_green, check_finished=False):
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in group_list:
        print("computing " + str(item) + " km")
        for i in range(len(item)):
            item[i] = item[i] + [path_green, check_finished]
        pool = Pool()
        r = pool.starmap_async(call_spgrn2020, item)
        r.get()
        pool.close()
        pool.join()
    update_green_info_lib_json(path_green, group_list[0][0][0], group_list[0][0][1])
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_spgrn2020_parallel_multi_nodes(path_green, check_finished=False):
    s = datetime.datetime.now()
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
