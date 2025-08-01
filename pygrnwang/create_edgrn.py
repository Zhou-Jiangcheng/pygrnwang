import os
import math

from .edgrn2inp import s as str_inp
from .utils import convert_earth_model_nd2inp, call_exe


def create_inp_edgrn2(
    path_green,
    obs_depth,
    grn_dist_range,
    grn_delta_dist,
    grn_source_depth_range,
    grn_delta_source_depth,
    wavenumber_sampling_rate=12,
    path_nd=None,
    earth_model_layer_num=None,
):
    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]
    lines_earth = lines[84:-1]
    lines_end = lines[-1]
    lines = lines[:84]

    lines[39] = "%.2f\n" % (obs_depth * 1e3)
    n_dist = math.ceil((grn_dist_range[1] - grn_dist_range[0]) / grn_delta_dist) + 1
    lines[40] = "%d %f %f\n" % (
        n_dist,
        grn_dist_range[0] * 1e3,
        grn_dist_range[1] * 1e3,
    )
    n_source_depth = (
        math.ceil(
            (grn_source_depth_range[1] - grn_source_depth_range[0])
            / grn_delta_source_depth
        )
        + 1
    )
    lines[41] = "%d %f %f\n" % (
        n_source_depth,
        grn_source_depth_range[0] * 1e3,
        grn_source_depth_range[1] * 1e3,
    )
    lines[51] = "%d\n" % wavenumber_sampling_rate
    path_edgrn_obs_dep = str(os.path.join(path_green, "edgrn2", "%.2f" % obs_depth, ""))
    lines[61] = "'%s' 'edgrn.ss' 'edgrn.ds' 'edgrn.cl'\n" % (path_edgrn_obs_dep)

    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    for i in range(len(lines_earth)):
        temp = lines_earth[i].split()
        lines_earth[i] = (
            temp[0] + " " + " ".join("%.4f" % (float(_) * 1e3) for _ in temp[1:-2])
        )
        lines_earth[i] = lines_earth[i] + "\n"
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    else:
        lines_earth = lines_earth[:earth_model_layer_num]
    lines[80] = "%d\n" % earth_model_layer_num

    lines = lines + lines_earth + [lines_end]
    with open(os.path.join(path_edgrn_obs_dep, "grn.inp"), "w") as fw:
        fw.writelines(lines)


def call_edgrn2(obs_depth, path_green, check_finished=False):
    os.chdir(path_green)
    sub_sub_dir = str(os.path.join("edgrn2", "%.2f" % obs_depth))
    if (
        check_finished
        and os.path.exists(os.path.join(sub_sub_dir, ".finished"))
        and len(os.listdir(sub_sub_dir)) > 2
    ):
        return None
    path_inp = str(os.path.join(sub_sub_dir, "grn.inp"))
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
        path_green=path_green,
        path_inp=path_inp,
        path_finished=path_finished,
        name="edgrn2",
    )
    return output
