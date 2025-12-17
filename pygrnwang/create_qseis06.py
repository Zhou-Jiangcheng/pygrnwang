import os
import math

import numpy as np
import pandas as pd

from .qseis06inp import s as str_inp
from .utils import convert_earth_model_nd2inp, call_exe


def create_dir_qseis06(
    path_green,
    event_depth,
    receiver_depth,
    dist_range,
    delta_dist,
    N_each_group=100,
    order=0,
):
    sub_dir = str(
        os.path.join(path_green, "%.2f" % event_depth, "%.2f" % receiver_depth)
    )
    N_dist = math.ceil((dist_range[1] - dist_range[0]) / delta_dist) + 1
    N_dist_group = math.ceil(N_dist / N_each_group)
    for n in range(N_dist_group):
        if receiver_depth > 0:
            for o in range(2 * order + 1):
                os.makedirs(os.path.join(sub_dir, "%d_%d" % (n, o)), exist_ok=True)
        else:
            for o in range(order + 1):
                os.makedirs(os.path.join(sub_dir, "%d_%d" % (n, o)), exist_ok=True)
    return N_dist, N_dist_group


def create_inp_qseis06(
    path_sub_dir,
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
    order=0,
):
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
        path_inp = os.path.join(path_sub_dir, "%d_%d" % (n, order), "grn.inp")
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    else:
        res = N_dist - (N_dist_group - 1) * N_each_group
        lines[44] = "%d\n" % res
        lines[45] = "%f %f\n" % (
            (dist_range[0] + (N_dist_group - 1) * N_each_group * delta_dist) * r_ratio,
            dist_range[1] * r_ratio,
        )
        path_inp = os.path.join(
            path_sub_dir, "%d_%d" % (N_dist_group - 1, order), "grn.inp"
        )
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    return path_inp


def create_inp_qseis06_points(
    path_green,
    event_depth,
    receiver_depth,
    n_group,
    points,
    time_window,
    sampling_interval,
    slowness_int_algorithm=0,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=2,
    anti_alias=0.01,
    free_surface=True,
    wavelet_duration=4,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
    order=0,
):
    path_sub_dir = str(
        os.path.join(path_green, "%.2f" % event_depth, "%.2f" % receiver_depth)
    )
    # when receiver depth is not 0, change dists in inp file
    r_ratio = (6371 - receiver_depth) / 6371
    points = points * r_ratio

    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]

    lines_earth = lines[207:-22]
    lines_end = lines[-22:]

    # SOURCE PARAMETERS
    lines[26] = "%f\n" % event_depth

    # RECEIVER PARAMETERS
    lines[42] = "%f\n" % receiver_depth
    lines[43] = "0 1\n"
    lines[44] = "%d\n" % len(points)
    lines[45] = " ".join("%.4f" % _ for _ in points) + "\n"
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

    path_inp = os.path.join(path_sub_dir, "%d_%d" % (n_group, order), "grn.inp")
    with open(path_inp, "w") as fw:
        fw.writelines(lines)
    return path_inp


def call_qseis06(
    event_depth, receiver_depth, n_group, order, path_green, check_finished=False
):
    sub_sub_dir = str(
        os.path.join(
            path_green,
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            "%d_%d" % (n_group, order),
        )
    )
    os.chdir(sub_sub_dir)
    path_inp = os.path.join(sub_sub_dir, "grn.inp")
    path_finished = os.path.join(sub_sub_dir, ".finished")

    if (
        check_finished
        and os.path.exists(path_finished)
        and len(os.listdir(sub_sub_dir)) > 2
    ):
        with open(path_finished, "r") as fr:
            output = fr.readlines()
        return output

    output = call_exe(
        path_inp=path_inp,
        path_finished=path_finished,
        name="qseis06",
    )
    return output


def convert_pd2bin_qseis06(path_greenfunc, remove=False):
    for com in ["tr", "tz", "tv"]:
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

    for com in ["tt"]:
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
