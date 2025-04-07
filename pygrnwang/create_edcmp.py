import os
import platform
import subprocess
import math


def create_inp_edcmp2(
        path_green,
        obs_depth,
        obs_x_range,
        obs_y_range,
        obs_delta_x,
        obs_delta_y,
        source_array_edcmp,
        layered=True,
        lam=30516224000,
        mu=33701888000,
):
    """

    :param path_green:
    :param obs_depth: km
    :param obs_x_range: km
    :param obs_y_range: km
    :param obs_delta_x: km
    :param obs_delta_y: km
    :param source_array_edcmp:
     [slip(m), x(km), y(km), z(km),
     strike(deg), dip(deg), rake(deg),
     sub_len_strike(km), sub_len_dip(km)]
    :param layered:
    :param lam:
    :param mu:
    :return:
    """
    path_inp = os.path.join(path_green, "edcmp2.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from .edcmp2inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]

    lines_after_sources = lines[97:]
    lines = lines[:96]

    nx = math.ceil((obs_x_range[1] - obs_x_range[0]) / obs_delta_x) + 1
    lines[45] = "%d %f %f\n" % (
        nx,
        obs_x_range[0] * 1e3,
        obs_x_range[1] * 1e3,
    )
    ny = math.ceil((obs_y_range[1] - obs_y_range[0]) / obs_delta_y) + 1
    lines[46] = "%d %f %f\n" % (
        ny,
        obs_y_range[0] * 1e3,
        obs_y_range[1] * 1e3,
    )
    path_edcmp_obs_dep = str(os.path.join(path_green, "edcmp2", "%.2f" % obs_depth, ""))
    lines[59] = "'%s'\n" % path_edcmp_obs_dep
    n_sources = len(source_array_edcmp)
    lines[91] = "%d\n" % n_sources

    lines_sources = []
    for i in range(len(source_array_edcmp)):
        si = "%d %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f %.2f\n" % (
            i + 1,
            source_array_edcmp[i, 0],
            source_array_edcmp[i, 1] * 1e3,
            source_array_edcmp[i, 2] * 1e3,
            source_array_edcmp[i, 3] * 1e3,
            source_array_edcmp[i, 7] * 1e3,
            source_array_edcmp[i, 8] * 1e3,
            source_array_edcmp[i, 4],
            source_array_edcmp[i, 5],
            source_array_edcmp[i, 6],
        )
        lines_sources.append(si)

    if layered:
        lines_after_sources[32] = "1\n"
        path_edgrn = os.path.join(path_green, "edgrn2", "%.2f" % obs_depth, "")
        lines_after_sources[33] = "'%s' 'edgrn.ss' 'edgrn.ds' 'edgrn.cl'\n" % (
            path_edgrn
        )
    else:
        lines_after_sources[32] = "0\n"
        lines_after_sources[33] = "%f %f %f\n" % (obs_depth * 1e3, lam, mu)

    lines = lines + lines_sources + lines_after_sources
    with open(os.path.join(path_edcmp_obs_dep, "grn.inp"), "w") as fw:
        fw.writelines(lines)


def call_edcmp2(obs_depth, path_green, check_finished=False):
    sub_sub_dir = str(os.path.join(path_green, "edcmp2", "%.2f" % obs_depth))
    if check_finished and os.path.exists(os.path.join(sub_sub_dir, ".finished")) and\
        len(os.listdir(sub_sub_dir))>2:
        return None
    path_inp = str(os.path.join(sub_sub_dir, "grn.inp"))

    if platform.system() == "Windows":
        edcmp_process = subprocess.Popen(
            [os.path.join(path_green, "edcmp2.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        edcmp_process.communicate(str.encode(path_inp))
    else:
        edcmp_process = subprocess.Popen(
            [os.path.join(path_green, "edcmp2.bin")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        edcmp_process.communicate(str.encode(path_inp))
    with open(os.path.join(sub_sub_dir, ".finished"), "w") as fw:
        fw.writelines([])
