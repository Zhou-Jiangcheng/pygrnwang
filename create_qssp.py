import datetime
import os
import platform
import subprocess
import shutil

import numpy as np

from pygrnwang.utils import convert_earth_model_nd2inp


m2d = 8.993216059187304e-06


mt_com_list = ["mrr", "mtt", "mpp", "mrt", "mrp", "mtp"]


def create_dir(event_depth, receiver_depth, path_green):
    # 新建文件夹并返回路径,若文件夹已存在则直接返回路径
    # if not os.path.exists(os.path.join(path_green, "GreenFunc")):
    #     os.mkdir(os.path.join(path_green, "GreenFunc"))
    # if not os.path.exists(os.path.join(path_green, "GreenSpec")):
    #     os.mkdir(os.path.join(path_green, "GreenSpec"))
    path_func = os.path.join(
        path_green, "GreenFunc", "%.1f" % event_depth, "%.1f" % receiver_depth, ""
    )
    os.makedirs(path_func, exist_ok=True)
    for mt_com in mt_com_list:
        os.makedirs(os.path.join(str(path_func), mt_com), exist_ok=True)

    path_spec = os.path.join(
        path_green, "GreenSpec", "%.1f" % event_depth, "%.1f" % receiver_depth, ""
    )
    os.makedirs(path_spec, exist_ok=True)
    return path_func, path_spec


def create_points(dist_range, delta_dist):
    """
    :param dist_range: km
    :param delta_dist: km
    """
    rmax = dist_range[1]
    rmin = dist_range[0]
    d = delta_dist
    n = round(np.ceil(2 * rmax / d))
    points = []
    for i in range(n):
        for j in range(n):
            points.append([i * d, j * d])
    points = np.array(points)
    points = points - rmax
    r_points = np.sqrt(np.sum(points ** 2, axis=1))
    points = points[(r_points <= rmax) & (r_points >= rmin)]
    points = points * 1e3 * m2d
    return points


def create_locs(points, time_reduction):
    lines_locs = []
    for i in range(len(points)):
        loc = "%f %f 'Loc%d' %f\n" % (
            points[i, 0], points[i, 1], i, time_reduction)
        lines_locs.append(loc)
    return lines_locs


def create_inp(
        mt_com,
        event_depth,
        receiver_depth,
        path_green,
        spec_time_window,
        time_window,
        time_reduction,
        dist_range,
        delta_dist_range,
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
        cal_gf,
        output_observables: list,
        path_nd=None,
        earth_model_layer_num=None,
):
    path_func = str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.1f" % event_depth,
            "%.1f" % receiver_depth,
            mt_com,
            "",
        )
    )
    path_spec = str(
        os.path.join(
            path_green,
            "GreenSpec",
            "%.1f" % event_depth,
            "%.1f" % receiver_depth,
            "",
        )
    )
    path_inp = os.path.join(path_green, "qssp2020.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from qssp2020inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]
    last_line = [lines[-1]]
    lines_earth_head = lines[141:155]
    lines_earth = lines[155:-1]
    lines = lines[:140]  # cutoff locs

    lines[24] = "%.1f\n" % receiver_depth

    lines[52] = "%.2f  %.2f\n" % (spec_time_window, sampling_interval)
    lines[53] = "%.2f\n" % max_frequency
    lines[54] = "%.2f\n" % max_slowness
    lines[55] = "%.2f\n" % anti_alias
    lines[56] = "%d %.2f %.2f\n" % (
        turning_point_filter,
        turning_point_d1,
        turning_point_d2,
    )
    lines[57] = "6371.0 %d\n" % free_surface_filter

    lines[65] = "%f %d\n" % (gravity_fc, gravity_harmonic)

    lines[75] = "%d %d %d %d\n" % (
        cal_sph, cal_tor, min_harmonic, max_harmonic)

    lines[86] = "1 0.0 '%s'\n" % path_spec
    lines[87] = "%.1f 'Green_%.1fkm' %d\n" % (event_depth, event_depth, cal_gf)

    lines[111] = "1 1\n"
    mt = [0 for _ in range(6)]
    if mt_com in mt_com_list:
        ind_mt = mt_com_list.index(mt_com)
        mt[ind_mt] = 1
    lines[112] = "1.0 %.1f %.1f %.1f %.1f %.1f %.1f 0.0 0.0 %.1f 0.0 0.0\n" % (
        mt[0],
        mt[1],
        mt[2],
        mt[3],
        mt[4],
        mt[5],
        event_depth,
    )

    lines[134] = "%d %d %d %d %d %d %d %d %d %d %d\n" % (
        output_observables[0],
        output_observables[1],
        output_observables[2],
        output_observables[3],
        output_observables[4],
        output_observables[5],
        output_observables[6],
        output_observables[7],
        output_observables[8],
        output_observables[9],
        output_observables[10],
    )
    lines[135] = "'%s'\n" % path_func
    lines[136] = "%f\n" % time_window
    lines[137] = "0 0 0\n"
    lines[138] = "0 %f\n" % max_slowness

    points = create_points(dist_range=dist_range, delta_dist=delta_dist_range)
    n_locs = len(points)
    lines[139] = "%s\n" % n_locs
    lines_locs = create_locs(points=points, time_reduction=time_reduction)
    lines = lines + lines_locs

    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    lines_earth_head[7] = "%d  0\n" % earth_model_layer_num
    if mt_com == "spec":
        path_inp = os.path.join(
            path_spec,
            "%.1f_%.1f_spec.inp" % (event_depth, receiver_depth),
        )
    elif mt_com in mt_com_list:
        ind_mt = mt_com_list.index(mt_com)
        path_inp = os.path.join(
            path_func,
            "%.1f_%.1f_%s.inp" % (
                event_depth, receiver_depth, mt_com_list[ind_mt]),
        )
    else:
        return ValueError("mt_com wrong!!!")
    with open(path_inp, "w") as fw:
        fw.writelines(lines + lines_earth_head + lines_earth + last_line)
    return path_inp


def call_qssp2020(event_depth, receiver_depth, mt_com, path_green):
    os.chdir(path_green)
    if mt_com == "spec":
        path_inp = str(
            os.path.join(
                path_green,
                "GreenSpec",
                "%.1f" % event_depth,
                "%.1f" % receiver_depth,
                "%.1f_%.1f_spec.inp" % (event_depth, receiver_depth),
            )
        ).replace("'", "")
    elif mt_com in mt_com_list:
        path_inp = str(
            os.path.join(
                path_green,
                "GreenFunc",
                "%.1f" % event_depth,
                "%.1f" % receiver_depth,
                mt_com,
                "%.1f_%.1f_%s.inp" % (event_depth, receiver_depth, mt_com),
            )
        ).replace("'", "")
    else:
        return ValueError("mt_com wrong!!!")
    path_bin_call = os.path.join(path_green, "qssp2020.bin")
    qssp_process = subprocess.Popen(
        [path_bin_call], stdin=subprocess.PIPE, stdout=subprocess.PIPE
    )
    qssp_process.communicate(str.encode(path_inp))
    # if platform.system() == "Windows":
    #     spgrn_process = subprocess.Popen(
    #         [os.path.join(path_green, "qssp2020.exe")],
    #         stdin=subprocess.PIPE,
    #         stdout=subprocess.PIPE)
    #     spgrn_process.communicate(str.encode(path_inp))
    # else:
    #     try:
    #         spgrn_process = subprocess.Popen(
    #             [os.path.join(path_green, "qssp2020.bin")],
    #             stdin=subprocess.PIPE,
    #             stdout=subprocess.PIPE)
    #         spgrn_process.communicate(str.encode(path_inp))
    #     except Exception as e:
    #         print(e)
    #         raise ("this system is not supported yet, \
    #             please compile the source code of qssp2020, \
    #             copy and replace the qssp2020.bin file ")


def create_grnlib(
        event_depth,
        receiver_depth,
        path_green,
        path_bin,
        spec_time_window,
        time_window,
        time_reduction,
        dist_range,
        delta_dist_range,
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
        output_observables: list,
        path_nd=None,
        earth_model_layer_num=None,
):
    print("creating green func lib recv_depth=%d" % event_depth)
    s = datetime.datetime.now()
    create_dir(event_depth, receiver_depth, path_green)
    create_inp(
        "spec",
        event_depth,
        receiver_depth,
        path_green,
        spec_time_window,
        time_window,
        time_reduction,
        dist_range,
        delta_dist_range,
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
        1,
        [0 for _ in range(11)],
        path_nd,
        earth_model_layer_num,
    )
    for mt_com in mt_com_list:
        create_inp(
            mt_com,
            event_depth,
            receiver_depth,
            path_green,
            spec_time_window,
            time_window,
            time_reduction,
            dist_range,
            delta_dist_range,
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
            0,
            output_observables,
            path_nd,
            earth_model_layer_num,
        )
    path_bin_call = os.path.join(path_green, "qssp2020.bin")
    if not os.path.exists(path_bin_call):
        shutil.copy(path_bin, path_bin_call)
    call_qssp2020(event_depth, receiver_depth, "spec", path_green)
    for mt_com in mt_com_list:
        call_qssp2020(event_depth, receiver_depth, mt_com, path_green)
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))
    print("done")


if __name__ == "__main__":
    output_observables_ = [0 for _ in range(11)]
    output_observables_[5] = 1
    create_grnlib(
        event_depth=10,
        receiver_depth=20,
        path_green="/e/qssp2020_green_lib_d10",
        path_bin="/home/zjc/python_works/pygrnwang/qssp2020.bin",
        spec_time_window=255,
        time_window=255,
        time_reduction=-10,
        dist_range=[0, 320],
        delta_dist_range=10,
        sampling_interval=1,
        max_frequency=0.2,
        max_slowness=0.4,
        anti_alias=0.01,
        turning_point_filter=0,
        turning_point_d1=0,
        turning_point_d2=0,
        free_surface_filter=1,
        gravity_fc=0,
        gravity_harmonic=0,
        cal_sph=1,
        cal_tor=1,
        min_harmonic=6000,
        max_harmonic=6000,
        output_observables=output_observables_,
        path_nd="/home/zjc/python_works/pygrnwang/turkey.nd",
        earth_model_layer_num=7,
    )
