import datetime
import os
import platform
import subprocess
import math

import numpy as np
import pandas as pd


def create_dir(event_depth, path_green, dist_range, delta_dist, N_each_group=100):
    sub_dir = str(os.path.join(path_green, "%.1f" % event_depth))
    N_dist = math.ceil((dist_range[1] - dist_range[0]) / delta_dist)
    N_dist_group = math.ceil(N_dist / N_each_group)
    if not os.path.exists(sub_dir):
        os.mkdir(sub_dir)
    sub_sub_dirs = []
    for n in range(N_dist_group):
        sub_sub_dir = os.path.join(sub_dir, "%d" % n)
        if not os.path.exists(sub_sub_dir):
            os.mkdir(sub_sub_dir)
        sub_sub_dirs.append(sub_sub_dir)
    return sub_dir, sub_sub_dirs, N_dist, N_dist_group


def create_greeninfo(
    event_depth,
    time_window,
    sampling_interval,
    dist_range,
    delta_dist,
    sub_dir,
    N_dist,
    N_dist_group,
    N_each_group,
):
    lines = ""
    lines += "event_depth: %f\n" % event_depth
    lines += "time_window: %f\n" % time_window
    sampling_num = 2 ** (math.ceil(math.log(time_window / sampling_interval, 2)))
    lines += "sampling_num: %d\n" % sampling_num
    lines += "sampling_interval: %f\n" % sampling_interval
    lines += "dist_range: %f %f\n" % (dist_range[0], dist_range[1])
    lines += "delta_dist: %f\n" % delta_dist
    lines += "N_dist: %d\n" % N_dist
    lines += "N_dist_group: %d\n" % N_dist_group
    lines += "N_each_group: %d\n" % N_each_group
    with open(os.path.join(sub_dir, "GreenInfo%.1f.dat" % event_depth), "w") as fw:
        fw.writelines(lines)


def create_inp(
    event_depth,
    path_green,
    time_window,
    sampling_interval,
    dist_range,
    delta_dist,
    sub_sub_dirs,
    N_dist,
    N_dist_group,
    isurf,
    rm_down=False,
    earth_model_layer_num=139,
    N_each_group=200,
    time_reduce_slowness=8,
):
    with open(os.path.join(path_green, "qseis06.inp"), "r") as fr:
        lines = fr.readlines()
    lines[26] = "%.1f\n" % event_depth
    lines[46] = "%.2f %.2f %d\n" % (
        0.0,
        time_window,
        2 ** (math.ceil(math.log(time_window / sampling_interval, 2))),
    )
    lines[47] = "%d %.1f\n" % (0, time_reduce_slowness)
    lines[103] = "%d\n" % isurf
    if rm_down and event_depth > 0.1:
        lines[105] = "1\n"
        lines[106] = "0.0 %.1f 2\n" % (event_depth - 0.1)
    else:
        lines[105] = "0\n"
        lines[106] = "#%.1f %.1f %d\n" % (0, 0, 2)
    lines[226] = "%d\n" % earth_model_layer_num

    for n in range(N_dist_group - 1):
        lines[44] = "%d\n" % N_each_group
        lines[45] = "%.2f %.2f\n" % (
            dist_range[0] + n * N_each_group * delta_dist,
            dist_range[0] + ((n + 1) * N_each_group - 1) * delta_dist,
        )
        path_inp = os.path.join(sub_sub_dirs[n], "%.1f_%d.inp" % (event_depth, n))
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    else:
        res = N_dist + 1 - (N_dist_group - 1) * N_each_group
        lines[44] = "%d\n" % res
        lines[45] = "%.2f %.2f\n" % (
            dist_range[0] + (N_dist_group - 1) * N_each_group * delta_dist,
            dist_range[0] + (N_dist - res) * delta_dist,
        )
        path_inp = os.path.join(
            sub_sub_dirs[-1], "%.1f_%d.inp" % (event_depth, N_dist_group - 1)
        )
        with open(path_inp, "w") as fw:
            fw.writelines(lines)


def call_qseis(event_depth, n_group, path_green):
    sub_sub_dir = str(os.path.join(path_green, "%.1f" % event_depth, "%d" % n_group))
    os.chdir(sub_sub_dir)
    path_inp = str(os.path.join(sub_sub_dir, "%.1f_%d.inp" % (event_depth, n_group)))

    if platform.system() == "Windows":
        qssp_process = subprocess.Popen(
            [os.path.join(path_green, "qseis06.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE)
        qssp_process.communicate(str.encode(path_inp))
    else:
        try:
            qssp_process = subprocess.Popen(
                [os.path.join(path_green, "qseis06.bin")],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
            qssp_process.communicate(str.encode(path_inp))
        except Exception as e:
            print(e)
            raise ("this system is not supported yet, \
                please compile the source code of qseis06, \
                copy and replace the qseis06.bin file ")
    convert2bin(sub_sub_dir)


def convert2bin(sub_sub_dir):
    print(sub_sub_dir)
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
            ex_z_df.values,
            ex_r_df.values,
            ss_z_df.values,
            ss_r_df.values,
            ss_t_df.values,
            ds_z_df.values,
            ds_r_df.values,
            ds_t_df.values,
            cl_z_df.values,
            cl_r_df.values,
        ]
    ).T
    print(time_series.shape)
    time_series = np.array(time_series, dtype=np.float32)
    time_series.tofile(os.path.join(sub_sub_dir, "grn.dat"))


def create_grnlib(
    event_depth,
    path_green,
    time_window,
    sampling_interval,
    dist_range,
    delta_dist,
    isurf,
    rm_down,
    earth_model_layer_num=139,
    N_each_group=100,
    time_reduce_slowness=8,
):
    print("creating green func lib event_depth=%d" % event_depth)
    s = datetime.datetime.now()
    sub_dir, sub_sub_dirs, N_dist, N_dist_group = create_dir(
        event_depth, path_green, dist_range, delta_dist, N_each_group
    )
    create_greeninfo(
        event_depth,
        time_window,
        sampling_interval,
        dist_range,
        delta_dist,
        sub_dir,
        N_dist,
        N_dist_group,
        N_each_group,
    )
    create_inp(
        event_depth,
        path_green,
        time_window,
        sampling_interval,
        dist_range,
        delta_dist,
        sub_sub_dirs,
        N_dist,
        N_dist_group,
        isurf,
        rm_down,
        earth_model_layer_num,
        N_each_group,
        time_reduce_slowness,
    )
    for n in range(N_dist_group):
        print(n)
        call_qseis(event_depth, n, path_green)
        # convert2bin(sub_sub_dirs[n])
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))
    print("done")


if __name__ == "__main__":
    pass
