import os
import pickle

from mpi4py import MPI
import numpy as np
import pandas as pd
from scipy.interpolate import LinearNDInterpolator

from pygrnwang.create_qssp import mt_com_list
from pygrnwang.focal_mechanism import convert_mt_axis

enz_list = ["ee", "en", "ez", "nn", "nz", "zz"]


def find_nearest(point_sta, points_geo):
    distances = np.sqrt(((points_geo - point_sta) ** 2).sum(axis=1))
    nearest_indice = np.argsort(distances)[0]
    nearest_point = points_geo[nearest_indice]
    return nearest_indice, nearest_point


def find_nearest4(point_sta, points_geo):
    distances = np.sqrt(((points_geo - point_sta) ** 2).sum(axis=1))
    nearest4_indices = np.argsort(distances)[:4]
    nearest4_points = points_geo[nearest4_indices]
    return nearest4_indices, nearest4_points


def interp4(points_4, values_vector_4, point_new):
    x = points_4[:, 0].flatten()
    y = points_4[:, 1].flatten()
    N_T = len(values_vector_4)
    v_new = np.zeros(N_T)
    for i in range(N_T):
        v = values_vector_4[i, :].flatten()
        f = LinearNDInterpolator(points_4, v)
        v = f(point_new[0], point_new[1])
        if np.isnan(v):
            nearest_indice = find_nearest(point_new, points_4)[0]
            v_new[i] = values_vector_4[i, nearest_indice]
        else:
            v_new[i] = v
    return v_new


def read_stress_tensor(
    path_green,
    event_depth,
    receiver_depth,
    points_green_geo,
    source,
    station,
    mt,
    interp=False,
) -> np.ndarray:
    """

    :param path_green:
    :param event_depth:
    :param receiver_depth:
    :param points_green_geo:
    :param source:
    :param station:
    :param mt:
    :param interp:
    :return: stress_enz, shape is (n*6)
    """
    point_sta = np.array(station[:2]) - np.array(source[:2])
    mt_rtp = convert_mt_axis(mt, "ned2rtp")
    stress_enz_rtp = []
    for i_rtp in range(6):
        path_func = str(
            os.path.join(
                path_green,
                "GreenFunc",
                "%.1f" % event_depth,
                "%.1f" % receiver_depth,
                mt_com_list[i_rtp],
                "",
            )
        )
        stress_all_1mt_com = []
        for enz_str in enz_list:
            if os.path.exists(os.path.join(path_func, "_stress_%s.npy" % enz_str)):
                stress_1enz_1mt_com = np.load(
                    os.path.join(path_func, "_stress_%s.npy" % enz_str)
                )
                stress_all_1mt_com.append(stress_1enz_1mt_com)
            else:
                stress_1enz_1mt_com = pd.read_csv(
                    os.path.join(path_func, "_stress_%s.dat" % enz_str), sep="\\s+"
                ).to_numpy()
                np.save(
                    os.path.join(path_func, "_stress_%s.npy" % enz_str),
                    stress_1enz_1mt_com,
                )
                stress_1enz_1mt_com = np.load(
                    os.path.join(path_func, "_stress_%s.npy" % enz_str)
                )
                stress_all_1mt_com.append(stress_1enz_1mt_com)
        N_T = stress_all_1mt_com[0].shape[0]
        if interp:
            nearest4_indices, nearest4_points = find_nearest4(
                point_sta, points_green_geo
            )
            stress_inerp_1mt_com = np.zeros((N_T, 6))
            for i_stress in range(6):
                stress_n4 = np.zeros((N_T, 4))
                for i_near in range(4):
                    stress_n4[:, i_near] = stress_all_1mt_com[i_stress][
                        :, nearest4_indices[i_near] + 1
                    ]
                stress_inerp_1mt_com[:, i_stress] = interp4(
                    nearest4_points, stress_n4, point_sta
                )
        else:
            nearest_indice, nearest_point = find_nearest(point_sta, points_green_geo)
            stress_inerp_1mt_com = np.zeros((N_T, 6))
            for i_stress in range(6):
                stress_inerp_1mt_com[:, i_stress] = stress_all_1mt_com[i_stress][
                    :, nearest_indice + 1
                ]
        stress_enz_rtp.append(stress_inerp_1mt_com)

    stress_enz = (
        stress_enz_rtp[0] * mt_rtp[0]
        + stress_enz_rtp[1] * mt_rtp[1]
        + stress_enz_rtp[2] * mt_rtp[2]
        + stress_enz_rtp[3] * mt_rtp[3]
        + stress_enz_rtp[4] * mt_rtp[4]
        + stress_enz_rtp[5] * mt_rtp[5]
    )
    return stress_enz


def convert_pd2np(path_green, event_dep, receiver_dep):
    for k in range(6):
        for l in range(6):
            path_dat = str(
                os.path.join(
                    path_green,
                    "GreenFunc",
                    "%.1f" % event_dep,
                    "%.1f" % receiver_dep,
                    mt_com_list[k],
                    "_stress_%s.dat" % enz_list[l],
                )
            )
            dat = pd.read_csv(path_dat, sep="\\s+").to_numpy()
            path_npy = path_dat[:-3] + "npy"
            np.save(path_npy, dat)


def convert_pd2np_single_thread(path_green, event_dep_list, receiver_dep_list):
    for i in range(len(event_dep_list)):
        for j in range(len(receiver_dep_list)):
            print(event_dep_list[i], receiver_dep_list[j])
            for k in range(6):
                for l in range(6):
                    path_dat = str(
                        os.path.join(
                            path_green,
                            "GreenFunc",
                            "%.1f" % event_dep_list[i],
                            "%.1f" % receiver_dep_list[j],
                            mt_com_list[k],
                            "_stress_%s.dat" % enz_list[l],
                        )
                    )
                    dat = pd.read_csv(path_dat, sep="\\s+").to_numpy()
                    path_npy = path_dat[:-3] + "npy"
                    np.save(path_npy, dat)


def convert_pd2np_mpi(path_green):
    with open(os.path.join(path_green, "group_list_func.pkl"), "rb") as fr:
        group_list_func = pickle.load(fr)
    N_all = 0
    for ind_group in range(len(group_list_func)):
        N_all = N_all + len(group_list_func[ind_group])
    for ind_group in range(len(group_list_func)):
        comm = MPI.COMM_WORLD
        processes_num = comm.Get_size()
        rank = comm.Get_rank()
        if processes_num != len(group_list_func[0]):
            raise ValueError(
                "processes_num is %d, item num in group is %d. \n"
                "Pleasse check the process num!"
                % (processes_num, len(group_list_func[0]))
            )
        print("ind_group:%d rank:%d" % (ind_group, rank))
        if ind_group * len(group_list_func[0]) + rank < N_all:
            convert_pd2np(
                path_green,
                group_list_func[ind_group][rank][0],
                group_list_func[ind_group][rank][1],
            )


if __name__ == "__main__":
    pass
