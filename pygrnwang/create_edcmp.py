import os

import numpy as np
import pandas as pd

from .edcmp2inp import s as str_inp
from .utils import call_exe

fm_base_list = (
    (315.0, 90.0, 0.0),  # [1,0,0,-1,0,0] m1
    (0.0, 90.0, 0.0),  # [0,1,0,0,0,0] mne
    (0.0, 0.0, 180.0),  # [0,0,1,0,0,0] mnd
    (180.0, 45.0, -90.0),  # [0,0,0,1,0,-1] m2
    (0.0, 0.0, 90.0),  # [0,0,0,0,1,0] med
)
output_name_list = ("disp", "strain", "stress", "tilt")


def create_inp_edcmp2(
    path_green: str,
    event_depth: float,
    obs_depth: float,
    dist_range: tuple[float, float],
    delta_dist: float,
    mt_ind: int,
    output_observables: tuple[int, int, int, int],
    layered: bool = True,
    lam: float = 30516224000.0,
    mu: float = 33701888000.0,
):
    """
    :param path_green:
    :param event_depth: km
    :param obs_depth: km
    :param dist_range: km
    :param delta_dist: km
    :param mt_ind: index of pure double-couple mt
    :param output_observables: disp, strain, stress, tilt
    :param layered:
    :param lam:
    :param mu:
    :return:
    """
    fm = fm_base_list[mt_ind]
    path_sub_dir = str(
        os.path.join(
            path_green,
            "edcmp2",
            "%.2f" % event_depth,
            "%.2f" % obs_depth,
            "%d" % mt_ind,
            "",
        )
    )
    os.makedirs(path_sub_dir, exist_ok=True)
    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]
    lines_after_sources = lines[97:]
    lines = lines[:96]

    lines[44] = "1\n"
    nx = len(np.arange(dist_range[0], dist_range[1] + delta_dist, delta_dist))
    lines[45] = "%d\n" % nx
    lines[46] = "(%f,0) (%f,0)\n" % (
        dist_range[0] * 1e3,
        dist_range[1] * 1e3,
    )

    lines[59] = "'%s'\n" % path_sub_dir
    lines[60] = "%d %d %d %d\n" % (
        output_observables[0],
        output_observables[1],
        output_observables[2],
        output_observables[3],
    )
    lines[91] = "1\n"

    lines_sources = [
        "1 1 0 0 %f 1 1 %.1f %.1f %.1f\n" % (event_depth * 1e3, fm[0], fm[1], fm[2])
    ]
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
    with open(os.path.join(path_sub_dir, "grn.inp"), "w") as fw:
        fw.writelines(lines)


def call_edcmp2(event_depth, obs_depth, mt_ind, path_green, check_finished=False):
    os.chdir(path_green)
    sub_sub_dir = str(
        os.path.join(
            path_green,
            "edcmp2",
            "%.2f" % event_depth,
            "%.2f" % obs_depth,
            "%d" % mt_ind,
            "",
        )
    )
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
        name="edcmp2",
    )
    return output


def convert_edcmp2(path_sub_dir, output_type_ind, remove=False):
    fname_df = str(
        os.path.join(path_sub_dir, "hs.%s" % output_name_list[output_type_ind])
    )
    df = pd.read_csv(
        fname_df,
        skiprows=3,
        sep="\\s+",
        header=None,
    )
    values_raw = df.to_numpy()[:, 2:]
    np.save(
        str(os.path.join(path_sub_dir, "%s.npy" % output_name_list[output_type_ind])),
        values_raw,
    )
    if remove:
        os.remove(fname_df)
    return values_raw
