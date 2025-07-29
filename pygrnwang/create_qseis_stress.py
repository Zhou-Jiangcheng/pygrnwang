import os
import platform
import subprocess
import math

import numpy as np
import pandas as pd

from .utils import convert_earth_model_nd2inp
from .qseis_stress_inp import s as str_inp


def create_dir_qseis_stress(
    path_green,
    event_depth,
    receiver_depth,
    dist_range,
    delta_dist,
    N_each_group=500,
):
    sub_dir = str(
        os.path.join(path_green, "%.2f" % event_depth, "%.2f" % receiver_depth)
    )
    N_dist = math.ceil((dist_range[1] - dist_range[0]) / delta_dist) + 1
    N_dist_group = math.ceil(N_dist / N_each_group)
    for n in range(N_dist_group):
        os.makedirs(os.path.join(sub_dir, "%d_0" % n), exist_ok=True)
    return N_dist, N_dist_group


def create_inp_qseis_stress(
    path_green,
    event_depth,
    receiver_depth,
    dist_range,
    delta_dist,
    N_dist,
    N_dist_group,
    N_each_group,
    time_window,
    sampling_interval,
    slowness_int_algorithm=0,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=12,
    anti_alias=0.01,
    free_surface=True,
    wavelet_duration=4,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
):
    path_sub_dir = str(
        os.path.join(path_green, "%.2f" % event_depth, "%.2f" % receiver_depth)
    )
    # when receiver depth is not 0, change dists in inp file
    r_ratio = (6371 - receiver_depth) / 6371
    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]

    lines_earth = lines[207:-22]
    lines_end = lines[-22:]

    # SOURCE PARAMETERS
    lines[26] = "%.2f\n" % event_depth

    # RECEIVER PARAMETERS
    lines[42] = "%.2f\n" % receiver_depth
    lines[43] = "1 1\n"
    lines[46] = "%f %f %d\n" % (
        0.0,
        time_window,
        round(time_window / sampling_interval + 1),
    )
    lines[47] = "%d %f\n" % (1, time_reduction_velo)

    # WAVENUMBER INTEGRATION PARAMETERS
    lines[66] = "%d\n" % slowness_int_algorithm
    if slowness_window is not None:
        lines[67] = "%f %f %f %f\n" % (
            slowness_window[0],
            slowness_window[1],
            slowness_window[2],
            slowness_window[3],
        )
    else:
        lines[67] = "0.0 0.0 0.0 0.0\n"
    lines[68] = "%f\n" % wavenumber_sampling_rate
    lines[69] = "%f\n" % anti_alias

    # OPTIONS FOR PARTIAL SOLUTIONS
    if free_surface:
        lines[103] = "0\n"
    else:
        lines[103] = "1\n"

    # SOURCE TIME FUNCTION (WAVELET) PARAMETERS (Note 3)
    lines[124] = "%d %d\n" % (wavelet_duration, wavelet_type)

    # GLOBAL MODEL PARAMETERS (Note 5)
    if flat_earth_transform:
        lines[191] = "1\n"
    else:
        lines[191] = "0\n"

    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    else:
        lines_earth = lines_earth[:earth_model_layer_num]
    lines[200] = "%d\n" % earth_model_layer_num
    lines = lines[:207] + lines_earth + lines_end

    for n in range(N_dist_group - 1):
        lines[44] = "%d\n" % N_each_group
        lines[45] = "%f %f\n" % (
            (dist_range[0] + n * N_each_group * delta_dist) * r_ratio,
            (dist_range[0] + ((n + 1) * N_each_group - 1) * delta_dist) * r_ratio,
        )
        path_inp = os.path.join(path_sub_dir, "%d_0" % n, "grn.inp")
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    else:
        res = N_dist - (N_dist_group - 1) * N_each_group
        lines[44] = "%d\n" % res
        lines[45] = "%f %f\n" % (
            (dist_range[0] + (N_dist_group - 1) * N_each_group * delta_dist) * r_ratio,
            dist_range[1] * r_ratio,
        )
        path_inp = os.path.join(path_sub_dir, "%d_0" % (N_dist_group - 1), "grn.inp")
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    return path_inp


def call_qseis_stress(
    event_depth, receiver_depth, n_group, path_green, check_finished=False
):
    sub_sub_dir = str(
        os.path.join(
            path_green,
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            "%d_0" % n_group,
        )
    )
    os.chdir(sub_sub_dir)
    if check_finished and os.path.exists(os.path.join(sub_sub_dir, ".finished")):
        return None
    path_inp = str(os.path.join(sub_sub_dir, "grn.inp"))

    if platform.system() == "Windows":
        qssp_process = subprocess.Popen(
            [os.path.join(path_green, "qseis_stress.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        qssp_process.communicate(str.encode(path_inp))
    else:
        qssp_process = subprocess.Popen(
            [os.path.join(path_green, "qseis_stress.bin")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        qssp_process.communicate(str.encode(path_inp))
    with open(os.path.join(sub_sub_dir, ".finished"), "w") as fw:
        fw.writelines([])
        return None


def convert_pd2bin_qseis_stress(path_greenfunc, remove=False):
    for com in ["tr", "tz", "tv", "trr", "szz", "szr"]:
        ex_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "ex.%s" % com)), sep="\\s+"
        ).to_numpy()
        ss_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "ss.%s" % com)), sep="\\s+"
        ).to_numpy()
        ds_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "ds.%s" % com)), sep="\\s+"
        ).to_numpy()
        cl_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "cl.%s" % com)), sep="\\s+"
        ).to_numpy()
        time_series_com = np.concatenate(
            [
                ex_com[:, 1:],
                ss_com[:, 1:],
                ds_com[:, 1:],
                cl_com[:, 1:],
            ]
        )
        time_series_com = np.array(time_series_com, dtype=np.float32)
        np.save(os.path.join(path_greenfunc, "grn_%s.npy" % com), time_series_com)
        if remove:
            os.remove(os.path.join(path_greenfunc, "ex.%s" % com))
            os.remove(os.path.join(path_greenfunc, "ss.%s" % com))
            os.remove(os.path.join(path_greenfunc, "ds.%s" % com))
            os.remove(os.path.join(path_greenfunc, "cl.%s" % com))

    for com in ["tt", "ttr", "szt"]:
        ss_r = pd.read_csv(
            str(os.path.join(path_greenfunc, "ss.%s" % com)), sep="\\s+"
        ).to_numpy()
        ds_r = pd.read_csv(
            str(os.path.join(path_greenfunc, "ds.%s" % com)), sep="\\s+"
        ).to_numpy()
        time_series_com = np.concatenate(
            [
                ss_r[:, 1:],
                ds_r[:, 1:],
            ]
        )
        time_series_com = np.array(time_series_com, dtype=np.float32)
        np.save(os.path.join(path_greenfunc, "grn_%s.npy" % com), time_series_com)
        if remove:
            os.remove(os.path.join(path_greenfunc, "ss.%s" % com))
            os.remove(os.path.join(path_greenfunc, "ds.%s" % com))


if __name__ == "__main__":
    pass
