import os
import platform
import subprocess
import math

import numpy as np
import pandas as pd

from .utils import convert_earth_model_nd2inp


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
        path_green,
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
        time_reduction_velo,
        wavenumber_sampling_rate=2,
        anti_alias=0.01,
        free_surface=True,
        wavelet_duration=4,
        flat_earth_transform=True,
        path_nd=None,
        earth_model_layer_num=None,
        order=0,
):
    path_inp = os.path.join(path_green, "qseis06.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from pygrnwang.qseis06inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]

    lines_earth = lines[207:-22]
    lines_end = lines[-22:]

    lines[26] = "%.2f\n" % event_depth
    lines[42] = "%.2f\n" % receiver_depth
    lines[43] = '1 1\n'

    lines[46] = "%.2f %.2f %d\n" % (
        0.0,
        time_window,
        2 ** (math.ceil(math.log(time_window / sampling_interval + 1, 2))),
    )
    lines[47] = "%d %.2f\n" % (1, time_reduction_velo)

    lines[68] = "%.2f\n" % wavenumber_sampling_rate
    lines[69] = "%.2f\n" % anti_alias

    if free_surface:
        lines[103] = "0\n"
    else:
        lines[103] = "1\n"
    lines[124] = "%d 1\n" % wavelet_duration

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
        lines[45] = "%.2f %.2f\n" % (
            dist_range[0] + n * N_each_group * delta_dist,
            dist_range[0] + ((n + 1) * N_each_group - 1) * delta_dist,
        )
        path_inp = os.path.join(path_sub_dir, "%d_%d" % (n, order), "grn.inp")
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    else:
        res = N_dist - (N_dist_group - 1) * N_each_group
        lines[44] = "%d\n" % res
        lines[45] = "%.2f %.2f\n" % (
            dist_range[0] + (N_dist_group - 1) * N_each_group * delta_dist,
            dist_range[1],
        )
        path_inp = os.path.join(
            path_sub_dir, "%d_%d" % (N_dist_group - 1, order), "grn.inp"
        )
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    return path_inp


def create_inp_qseis06_points(
        path_green,
        path_sub_dir,
        event_depth,
        receiver_depth,
        n_group,
        points,
        time_window,
        sampling_interval,
        time_reduction_velo,
        wavenumber_sampling_rate=2,
        anti_alias=0.01,
        free_surface=True,
        wavelet_duration=4,
        flat_earth_transform=True,
        path_nd=None,
        earth_model_layer_num=None,
        order=0,
):
    path_inp = os.path.join(path_green, "qseis06.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from pygrnwang.qseis06inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]

    lines_earth = lines[207:-22]
    lines_end = lines[-22:]

    lines[26] = "%f\n" % event_depth
    lines[42] = "%f\n" % receiver_depth

    lines[43] = '0 1\n'
    lines[44] = '%d\n' % len(points)
    lines[45] = " ".join("%.4f" % _ for _ in points)+"\n"

    lines[46] = "%.2f %.2f %d\n" % (
        0.0,
        time_window,
        2 ** (math.ceil(math.log(time_window / sampling_interval + 1, 2))),
    )
    lines[47] = "%d %.2f\n" % (1, time_reduction_velo)

    lines[68] = "%.2f\n" % wavenumber_sampling_rate
    lines[69] = "%.2f\n" % anti_alias

    if free_surface:
        lines[103] = "0\n"
    else:
        lines[103] = "1\n"
    lines[124] = "%d 1\n" % wavelet_duration

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

    path_inp = os.path.join(
        path_sub_dir, "%d_%d" % (n_group, order), "grn.inp"
    )
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
    if check_finished and os.path.exists(os.path.join(sub_sub_dir, ".finished")):
        return None
    path_inp = str(os.path.join(sub_sub_dir, "grn.inp"))

    if platform.system() == "Windows":
        qssp_process = subprocess.Popen(
            [os.path.join(path_green, "qseis06.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        qssp_process.communicate(str.encode(path_inp))
    else:
        qssp_process = subprocess.Popen(
            [os.path.join(path_green, "qseis06.bin")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        qssp_process.communicate(str.encode(path_inp))
    with open(os.path.join(sub_sub_dir, ".finished"), "w") as fw:
        fw.writelines([])


def convert_pd2bin_qseis06(sub_sub_dir):
    ex_z_df = pd.read_csv(os.path.join(sub_sub_dir, "ex.tz"), sep="\\s+")
    ex_r_df = pd.read_csv(os.path.join(sub_sub_dir, "ex.tr"), sep="\\s+")
    ss_z_df = pd.read_csv(os.path.join(sub_sub_dir, "ss.tz"), sep="\\s+")
    ss_r_df = pd.read_csv(os.path.join(sub_sub_dir, "ss.tr"), sep="\\s+")
    ss_t_df = pd.read_csv(os.path.join(sub_sub_dir, "ss.tt"), sep="\\s+")
    ds_z_df = pd.read_csv(os.path.join(sub_sub_dir, "ds.tz"), sep="\\s+")
    ds_r_df = pd.read_csv(os.path.join(sub_sub_dir, "ds.tr"), sep="\\s+")
    ds_t_df = pd.read_csv(os.path.join(sub_sub_dir, "ds.tt"), sep="\\s+")
    cl_z_df = pd.read_csv(os.path.join(sub_sub_dir, "cl.tz"), sep="\\s+")
    cl_r_df = pd.read_csv(os.path.join(sub_sub_dir, "cl.tr"), sep="\\s+")
    time_series = np.concatenate(
        [
            ex_z_df.values[:, 1:],
            ex_r_df.values[:, 1:],
            ss_z_df.values[:, 1:],
            ss_r_df.values[:, 1:],
            ss_t_df.values[:, 1:],
            ds_z_df.values[:, 1:],
            ds_r_df.values[:, 1:],
            ds_t_df.values[:, 1:],
            cl_z_df.values[:, 1:],
            cl_r_df.values[:, 1:],
        ]
    ).T
    time_series = np.array(time_series, dtype=np.float32)
    time_series.tofile(os.path.join(sub_sub_dir, "grn.bin"))


if __name__ == "__main__":
    pass
