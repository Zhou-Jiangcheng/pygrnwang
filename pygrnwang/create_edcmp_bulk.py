import os
import shutil
import platform
import pickle
import json
import datetime
import math

import numpy as np
from tqdm import tqdm
from mpi4py import MPI
from multiprocessing import Pool

from .create_edcmp import create_inp_edcmp2, call_edcmp2
from .utils import group, cal_max_dist_from_2d_points
from .geo import convert_sub_faults_geo2ned, d2m


def pre_process_edcmp2(
        processes_num,
        path_green,
        path_bin,
        obs_depth_list,
        obs_x_range,
        obs_y_range,
        obs_delta_x,
        obs_delta_y,
        source_array,
        source_ref=None,
        obs_ref=None,
        layered=True,
        lam=30516224000,
        mu=33701888000,
):
    """
    Pre-processes input data for edcmp2, updates Green's function library settings,
    and prepares input files for each observation depth.

    :param processes_num: Number of parallel processes to be used.
    :param path_green: Directory path for the Green's function library.
    :param path_bin: Path to the edcmp2 binary file to be used if not present in path_green.
    :param obs_depth_list: List of observation depths (in km) for which inputs are to be created.
    :param obs_x_range: List or tuple [min, max] of observation x-coordinates (latitude in degrees).
    :param obs_y_range: List or tuple [min, max] of observation y-coordinates (longitude in degrees).
    :param obs_delta_x: Sampling interval in the x-direction (in degrees).
    :param obs_delta_y: Sampling interval in the y-direction (in degrees).
    :param source_array: NumPy array with each row containing source parameters:
                         [slip (m), lat (deg), lon (deg), depth (km), strike (deg), dip (deg), rake (deg)].
    :param source_ref: Optional reference coordinate for sources [lat, lon]. If provided,
                       source coordinates will be converted from geographic to NED coordinates.
    :param obs_ref: Optional reference coordinate for observations [lat, lon]. If provided,
                    the observation coordinate ranges and intervals are converted to km.
    :param layered: Boolean flag indicating whether the Green's function library is for a layered model.
    :param lam: Lamé parameter λ (default: 30516224000).
    :param mu: Shear modulus μ (default: 33701888000).
    :return: group_list_edcmp, a grouping of observation depths according to processes_num,
             which is also saved to a pickle file in path_green.
    """
    # Print status (commented out)
    # print("preprocessing edcmp2")

    # Create a copy of the source array.
    source_array_new = source_array.copy()

    # Determine the correct binary filename based on the operating system.
    if platform.system() == "Windows":
        path_bin_call = os.path.join(path_green, "edcmp2.exe")
    else:
        path_bin_call = os.path.join(path_green, "edcmp2.bin")

    # If the binary does not exist in the Green's function directory, copy it.
    if not os.path.exists(path_bin_call):
        shutil.copy(path_bin, path_bin_call)

    # Load the current Green's function library information.
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)

    # If a source reference is provided, convert source coordinates.
    if source_ref is not None:
        source_array_new[:, 3] = source_array_new[:, 3] * 1e3  # Convert depth from km to m.
        sources_ned = convert_sub_faults_geo2ned(
            source_array_new[:, 1:4], np.concatenate([source_ref, np.zeros(1)]), True
        )
        source_array_new[:, 1:4] = sources_ned / 1e3  # Convert NED coordinates back to km.

    # If an observation reference is provided, convert observation ranges and intervals.
    if obs_ref is not None:
        obs_x_range = list((np.array(obs_x_range) - obs_ref[0]) * d2m / 1e3)  # Convert x-range to km.
        obs_y_range = list((np.array(obs_y_range) - obs_ref[1]) * d2m / 1e3)  # Convert y-range to km.
        obs_delta_x = obs_delta_x * d2m / 1e3  # Convert x sampling interval to km.
        obs_delta_y = obs_delta_y * d2m / 1e3  # Convert y sampling interval to km.

    # Check if the observation grid exceeds the range supported by the Green's function library.
    nx = math.ceil((obs_x_range[1] - obs_x_range[0]) / obs_delta_x) + 1
    ny = math.ceil((obs_y_range[1] - obs_y_range[0]) / obs_delta_y) + 1
    obs_points_x, obs_points_y = np.meshgrid(
        np.linspace(obs_x_range[0], obs_x_range[1], nx),
        np.linspace(obs_y_range[0], obs_y_range[1], ny),
    )
    obs_points = np.concatenate(
        [np.array([obs_points_x.flatten()]).T, np.array([obs_points_y.flatten()]).T],
        axis=1,
    )
    max_distance = cal_max_dist_from_2d_points(source_array_new[:, 1:3], obs_points)
    if max_distance > green_info["grn_dist_range"][1]:
        raise ValueError(
            "\nThe max distance in obs is %f km,\n"
            "exceeds the max value %f km\n"
            "in the Green's function library!"
            % (max_distance, green_info["grn_dist_range"][1])
        )

    # For each observation depth, create the corresponding edcmp2 input file.
    for obs_depth in obs_depth_list:
        sub_sub_dir = str(os.path.join(path_green, "edcmp2", "%.2f" % obs_depth))
        os.makedirs(sub_sub_dir, exist_ok=True)
        create_inp_edcmp2(
            path_green,
            obs_depth,
            obs_x_range,
            obs_y_range,
            obs_delta_x,
            obs_delta_y,
            source_array_new,
            layered,
            lam,
            mu,
        )

    # Update the green_info dictionary with new observation and model parameters.
    green_info["obs_x_range"] = obs_x_range
    green_info["obs_y_range"] = obs_y_range
    green_info["obs_delta_x"] = obs_delta_x
    green_info["obs_delta_y"] = obs_delta_y
    green_info["source_ref"] = source_ref
    green_info["obs_ref"] = obs_ref
    green_info["layered"] = layered
    green_info["lam"] = lam
    green_info["mu"] = mu
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
            os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    # Group the observation depths into subgroups for parallel processing.
    group_list_edcmp = group(obs_depth_list, processes_num)
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "wb") as fw:
        pickle.dump(group_list_edcmp, fw)  # type: ignore

    return group_list_edcmp

def compute_static_stress_edcmp2_sequential(path_green, check_finished=False):
    # s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    for item in tqdm(group_list_edcmp,
                     desc="Computing static stress"):
        for i in range(len(item)):
            # print("computing " + str(item[i]) + " km")
            call_edcmp2(item[i], path_green, check_finished)
    # e = datetime.datetime.now()
    # print("run time:%s" % str(e - s))


def compute_static_stress_edcmp2_parallel_single_node(path_green, check_finished=False):
    # s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    for item in tqdm(group_list_edcmp,
                     desc="Computing static stress"):
        for i in range(len(item)):
            item[i] = [item[i]] + [path_green, check_finished]
        pool = Pool()
        r = pool.starmap_async(call_edcmp2, item)
        r.get()
        pool.close()
        pool.join()
    # e = datetime.datetime.now()
    # print("run time:" + str(e - s))


def create_grnlib_edcmp2_parallel_multi_nodes(path_green, check_finished=False):
    s = datetime.datetime.now()
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
        call_edcmp2(
            obs_depth=group_list_edcmp[ind_group][rank],
            path_green=path_green,
            check_finished=check_finished,
        )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))
