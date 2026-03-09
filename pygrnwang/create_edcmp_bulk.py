import os
import shutil
import platform
import pickle
import json
import datetime
import math
from multiprocessing import Pool

import numpy as np
from tqdm import tqdm

try:
    from mpi4py import MPI
except:
    pass

from .create_edcmp import create_inp_edcmp2, call_edcmp2, convert_edcmp2
from .utils import group
from .geo import d2km, convert_sub_faults_geo2ned, cal_max_dist_from_2d_points


def _call_edcmp2_star(args):
    return call_edcmp2(*args)


def pre_process_edcmp2(
    processes_num: int,
    path_green: str,
    grn_source_depth_range,
    grn_source_delta_depth,
    grn_dist_range,
    grn_delta_dist,
    obs_depth_list,
    output_observables=(1, 0, 0, 0),
    layered=True,
    lam=30516224000,
    mu=33701888000,
):
    """
    Pre-processes input data for edcmp2, updates Green's function library settings,
    and prepares input files for each observation depth.

    :param processes_num: Number of parallel processes to be used.
    :param path_green: Directory path for the Green's function library.
    :param layered: Boolean flag indicating whether the Green's function library is for a layered model.
    :param lam: Lamé parameter λ (default: 30516224000).
    :param mu: Shear modulus μ (default: 33701888000).
    :return: group_list_edcmp, a grouping of observation depths according to processes_num,
             which is also saved to a pickle file in path_green.
    """
    # Load the current Green's function library information.
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)

    item_list = []
    event_depth_list = np.arange(
        grn_source_depth_range[0],
        grn_source_depth_range[1] + grn_source_delta_depth,
        grn_source_delta_depth,
    )
    for event_depth in event_depth_list:
        for obs_depth in obs_depth_list:
            for mt_ind in range(5):
                create_inp_edcmp2(
                    path_green=path_green,
                    event_depth=event_depth,
                    obs_depth=obs_depth,
                    dist_range=grn_dist_range,
                    delta_dist=grn_delta_dist,
                    mt_ind=mt_ind,
                    output_observables=output_observables,
                    layered=layered,
                    lam=lam,
                    mu=mu,
                )
                item_list.append([event_depth, obs_depth, mt_ind])

    # Update the green_info dictionary with new observation and model parameters.
    green_info["layered"] = layered
    if layered:
        green_info["lam"] = None
        green_info["mu"] = None
    else:
        green_info["lam"] = lam
        green_info["mu"] = mu
    green_info["output_observables"] = output_observables
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    # Group the observation depths into subgroups for parallel processing.
    group_list_edcmp = group(item_list, processes_num)
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "wb") as fw:
        pickle.dump(group_list_edcmp, fw)  # type: ignore

    return group_list_edcmp


def create_grnlib_edcmp2_sequential(path_green, check_finished=False):
    # s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    for item in tqdm(group_list_edcmp, desc="Computing static stress"):
        for i in range(len(item)):
            # print("computing " + str(item[i]) + " km")
            item[i] = item[i] + [path_green, check_finished]
            call_edcmp2(*item[i])
    # e = datetime.datetime.now()
    # print("run time:%s" % str(e - s))


def create_grnlib_edcmp2_parallel(path_green, check_finished=False, 
                                  convert_bulk=True, remove=False):
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list_edcmp = pickle.load(fr)
    tasks = []
    for grp in group_list_edcmp:
        for d in grp:
            tasks.append(tuple(d + [path_green, check_finished]))

    processes = None
    try:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            processes = json.load(fr).get("processes_num", None)
    except Exception:
        pass

    with Pool(processes=processes) as pool:
        for _ in tqdm(
            pool.imap_unordered(_call_edcmp2_star, tasks, chunksize=1),
            total=len(tasks),
            desc="Computing static lib",
        ):
            pass
    if convert_bulk:
        convert_pd2np_edcmp2_all(path_green, remove=remove)
    e = datetime.datetime.now()
    return e - s


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
            event_depth=group_list_edcmp[ind_group][rank][0],
            obs_depth=group_list_edcmp[ind_group][rank][1],
            mt_ind=group_list_edcmp[ind_group][rank][2],
            path_green=path_green,
            check_finished=check_finished,
        )
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


_EDCMP_OUTPUT_TYPE_NAMES = ["disp", "strain", "stress", "tilt"]
_EDCMP_CHA_NUM = {"disp": 3, "strain": 6, "stress": 6, "tilt": 2}


def convert_pd2np_edcmp2_all(path_green, remove=False):
    print("converting ascii files to npy files")
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    grn_source_depth_range = green_info["grn_source_depth_range"]
    grn_source_delta_depth = green_info["grn_source_delta_depth"]
    event_depth_list = np.arange(
        grn_source_depth_range[0],
        grn_source_depth_range[1] + grn_source_delta_depth,
        grn_source_delta_depth,
    )
    obs_depth_list = green_info["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]

    grn_dist_range = green_info["grn_dist_range"]
    grn_dist_delta = green_info["grn_delta_dist"]
    n_dist = len(
        np.arange(grn_dist_range[0], grn_dist_range[1] + grn_dist_delta, grn_dist_delta)
    )

    output_observables = np.nonzero(np.array(green_info["output_observables"]))[0]
    n_dep = len(event_depth_list)
    n_obs = len(obs_depth_list)

    # Pre-allocate one bulk array per active output_type:
    # shape (n_dep, n_obs, 5, n_dist, cha_num)
    bulk_data = {}
    for o in output_observables:
        ot = _EDCMP_OUTPUT_TYPE_NAMES[int(o)]
        bulk_data[ot] = np.zeros(
            (n_dep, n_obs, 5, n_dist, _EDCMP_CHA_NUM[ot]), dtype=np.float32
        )

    for i, event_depth in enumerate(event_depth_list):
        for j, obs_depth in enumerate(obs_depth_list):
            for k in range(5):
                path_sub_dir = os.path.join(
                    path_green,
                    "edcmp2",
                    "%.2f" % event_depth,
                    "%.2f" % obs_depth,
                    "%d" % k,
                )
                for o in output_observables:
                    v_ijko = convert_edcmp2(
                        path_sub_dir=path_sub_dir,
                        output_type_ind=int(o),
                        remove=remove,
                    )
                    bulk_data[_EDCMP_OUTPUT_TYPE_NAMES[int(o)]][i, j, k] = v_ijko

    for ot, arr in bulk_data.items():
        out_path = os.path.join(path_green, "edcmp2_%s.npy" % ot)
        np.save(out_path, arr)
        print("Saved %s, shape: %s" % (os.path.basename(out_path), str(arr.shape)))
