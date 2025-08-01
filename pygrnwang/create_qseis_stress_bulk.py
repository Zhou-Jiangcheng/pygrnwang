import os
import shutil
import platform
import pickle
import json
import datetime
from multiprocessing import Pool

from tqdm import tqdm
from mpi4py import MPI
import jpype

from .create_qseis_stress import (
    create_dir_qseis_stress,
    create_inp_qseis_stress,
    call_qseis_stress,
    convert_pd2bin_qseis_stress,
)
from .pytaup import create_tpts_table
from .utils import group, convert_earth_model_nd2nd_without_Q


def pre_process_qseis_stress(
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
    print("Preprocessing")
    os.makedirs(path_green, exist_ok=True)
    if platform.system() == "Windows":
        path_bin_call = os.path.join(path_green, "qseis_stress.exe")
    else:
        path_bin_call = os.path.join(path_green, "qseis_stress.bin")
    shutil.copy(path_bin, path_bin_call)

    N_dist, N_dist_group = None, None
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            # print("creating green func dir, info, inp for event_depth=%.2f receiver_depth=%.2f" %
            #       (event_depth, receiver_depth))
            N_dist, N_dist_group = create_dir_qseis_stress(
                path_green=path_green,
                event_depth=event_depth,
                receiver_depth=receiver_depth,
                dist_range=dist_range,
                delta_dist=delta_dist,
                N_each_group=N_each_group,
            )
            create_inp_qseis_stress(
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
                inp_list.append([event_dep, receiver_dep, nn])
    inp_list_sorted = sorted(inp_list, key=lambda x: abs(x[0] - x[1]))
    group_list = group(inp_list_sorted, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore
    return group_list


def create_grnlib_qseis_stress_sequential(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in tqdm(group_list, desc="Computing dynamic stress"):
        for i in range(len(item)):
            # print("computing " + str(item[i]))
            item[i] = item[i] + [path_green, check_finished]
            call_qseis_stress(*item[i])
    if convert_pd2bin:
        convert_pd2bin_qseis_stress_all(path_green, remove_pd)


def create_grnlib_qseis_stress_parallel(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    for item in tqdm(group_list, desc="Computing dynamic stress"):
        # print("computing " + str(item))
        for i in range(len(item)):
            item[i] = item[i] + [path_green, check_finished]
        pool = Pool()
        r = pool.starmap_async(call_qseis_stress, item)
        r.get()
        pool.close()
        pool.join()
    if convert_pd2bin:
        convert_pd2bin_qseis_stress_all(path_green, remove_pd)


def convert_pd2bin_qseis_stress_all(path_green, remove=False):
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
                    convert_pd2bin_qseis_stress(
                        os.path.join(sub_dir, sub_sub_dir), remove=remove
                    )


def create_grnlib_qseis_stress_parallel_multi_nodes(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True
):
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
    call_qseis_stress(
        event_depth=group_list[ind_group][ind_para][0],
        receiver_depth=group_list[ind_group][ind_para][1],
        n_group=group_list[ind_group][ind_para][2],
        path_green=path_green,
        check_finished=check_finished,
    )
    if convert_pd2bin:
        convert_pd2bin_qseis_stress_all(path_green, remove_pd)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


if __name__ == "__main__":
    pass
