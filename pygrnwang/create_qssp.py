import os

import numpy as np
import pandas as pd

from .geo import d2m
from .utils import convert_earth_model_nd2inp, call_exe

mt_com_list = ["mrr", "mtt", "mpp", "mrt", "mrp", "mtp"]
output_type_list = [
    "disp",
    "velo",
    "acce",
    "rota",
    "rota_rate",
    "stress",
    "stress_rate",
    "strain",
    "strain_rate",
    "gravitation",
    "gravimeter",
]


def create_dir_qssp2020(event_depth, receiver_depth, path_green):
    path_func = os.path.join(
        path_green, "GreenFunc", "%.2f" % event_depth, "%.2f" % receiver_depth, ""
    )
    os.makedirs(path_func, exist_ok=True)
    for mt_com in mt_com_list:
        os.makedirs(os.path.join(str(path_func), mt_com), exist_ok=True)

    path_spec = os.path.join(
        path_green, "GreenSpec", "%.2f" % event_depth, "%.2f" % receiver_depth, ""
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
    points_lat = np.linspace(rmin, rmax, round(np.ceil((rmax - rmin) / delta_dist)) + 1)
    points_lat = points_lat * 1e3 / d2m
    points = np.concatenate([points_lat, np.zeros_like(points_lat)])
    points = points.reshape((2, -1)).T
    return points


def create_locs(points, time_reduction):
    lines_locs = []
    for i in range(len(points)):
        loc = "%f %f 'Loc%d' %f\n" % (points[i, 0], points[i, 1], i, time_reduction)
        lines_locs.append(loc)
    return lines_locs


def create_inp_qssp2020(
    mt_com,
    path_green,
    event_depth,
    receiver_depth,
    spec_time_window,
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
    source_radius,
    cal_gf,
    source_duration,
    output_observables: list,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    path_func = str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            mt_com,
            "",
        )
    )
    path_spec = str(
        os.path.join(
            path_green,
            "GreenSpec",
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            "",
        )
    )
    path_inp = os.path.join(path_green, "qssp2020.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from .qssp2020inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]
    last_line = [lines[-1]]
    lines_earth_head = lines[141:155]
    lines_earth = lines[155:-1]
    lines = lines[:141]  # cutoff locs

    lines[24] = "%.2f\n" % receiver_depth

    lines[52] = "%f  %f\n" % (spec_time_window, sampling_interval)
    lines[53] = "%f\n" % max_frequency
    lines[54] = "%f\n" % max_slowness
    lines[55] = "%f\n" % anti_alias
    lines[56] = "%d %f %f\n" % (
        turning_point_filter,
        turning_point_d1,
        turning_point_d2,
    )
    lines[57] = "6371.0 %d\n" % free_surface_filter

    lines[65] = "%f %d\n" % (gravity_fc, gravity_harmonic)

    lines[75] = "%d %d %d %d\n" % (cal_sph, cal_tor, min_harmonic, max_harmonic)

    lines[86] = "1 %f '%s'\n" % (source_radius, path_spec)
    lines[87] = "%.2f 'Green_%.2fkm' %d\n" % (event_depth, event_depth, cal_gf)

    lines[111] = "1 1\n"
    mt = [0 for _ in range(6)]
    if mt_com in mt_com_list:
        ind_mt = mt_com_list.index(mt_com)
        mt[ind_mt] = 1
    lines[112] = (
        "1.0 "
        + " ".join("%f" % mt[_] for _ in range(6))
        + " 0.0 0.0 %.2f 0.0 %f\n" % (event_depth, source_duration)
    )
    lines[135] = " ".join(["%d" % output_observables[_] for _ in range(11)]) + "\n"
    lines[136] = "'%s'\n" % path_func
    lines[137] = "%f\n" % time_window
    lines[138] = "0 0 0\n"
    lines[139] = "0 %f\n" % max_slowness

    points = create_points(dist_range=dist_range, delta_dist=delta_dist)
    n_locs = len(points)
    lines[140] = "%s\n" % n_locs
    lines_locs = create_locs(points=points, time_reduction=time_reduction)
    lines = lines + lines_locs

    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    lines_earth = lines_earth[:earth_model_layer_num]
    lines_earth_head[7] = "%d  %d\n" % (earth_model_layer_num, physical_dispersion)
    if mt_com == "spec":
        path_inp = os.path.join(
            path_spec,
            "spec.inp",
        )
    elif mt_com in mt_com_list:
        ind_mt = mt_com_list.index(mt_com)
        path_inp = os.path.join(
            path_func,
            "%s.inp" % mt_com_list[ind_mt],
        )
    else:
        return ValueError("mt_com wrong!!!")
    with open(path_inp, "w") as fw:
        fw.writelines(lines + lines_earth_head + lines_earth + last_line)
    return path_inp


def call_qssp2020(
    event_depth, receiver_depth, mt_com, path_green, check_finished=False
):
    os.chdir(path_green)
    if mt_com == "spec":
        path_inp = str(
            os.path.join(
                path_green,
                "GreenSpec",
                "%.2f" % event_depth,
                "%.2f" % receiver_depth,
                "spec.inp",
            )
        ).replace("'", "")
    elif mt_com in mt_com_list:
        path_inp = str(
            os.path.join(
                path_green,
                "GreenFunc",
                "%.2f" % event_depth,
                "%.2f" % receiver_depth,
                mt_com,
                "%s.inp" % mt_com,
            )
        ).replace("'", "")
    else:
        return ValueError("mt_com wrong!!!")

    sub_dir = os.path.dirname(path_inp)
    path_finished = os.path.join(sub_dir, ".finished")
    if (
        check_finished
        and os.path.exists(path_finished)
        and len(os.listdir(sub_dir)) > 2
    ):
        with open(path_finished, "r") as fr:
            output = fr.readlines()
        return output

    output = call_exe(
        path_green=path_green,
        path_inp=path_inp,
        path_finished=path_finished,
        name="qssp2020",
    )
    return output


def convert_pd2bin_qssp2020(path_green, event_depth, receiver_depth, output_type_ind):
    if output_type_list[output_type_ind] == "gravimeter":
        enz_list = [""]
    elif output_type_list[output_type_ind] in [
        "disp",
        "velo",
        "acce",
        "rota",
        "rota_rate",
        "gravitation",
    ]:
        enz_list = ["_e", "_n", "_z"]
    elif output_type_list[output_type_ind] in [
        "stress",
        "stress_rate",
        "strain",
        "strain_rate",
    ]:
        enz_list = ["_ee", "_en", "_ez", "_nn", "_nz", "_zz"]
    else:
        raise ValueError("output_type wrong")
    for k in range(6):
        for l in range(len(enz_list)):
            path_dat = str(
                os.path.join(
                    path_green,
                    "GreenFunc",
                    "%.2f" % event_depth,
                    "%.2f" % receiver_depth,
                    mt_com_list[k],
                    "_%s%s.dat" % (output_type_list[output_type_ind], enz_list[l]),
                )
            )
            dat = pd.read_csv(path_dat, sep="\\s+").to_numpy()
            dat = np.array(dat, dtype=np.float32)[:, 1:].T
            dat.tofile(os.path.join(path_dat[:-4], "%s.bin" % path_dat[:-4]))


if __name__ == "__main__":
    pass
