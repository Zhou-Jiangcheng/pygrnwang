import os

import numpy as np

from .spgrn2020inp import s as str_inp
from .read_green_info_spgrn import read_green_info_spgrn
from .utils import call_exe, convert_earth_model_nd2inp, check_file_size


def create_dir_spgrn(event_depth, receiver_depth, path_green):
    os.makedirs(os.path.join(path_green, "GreenFunc"), exist_ok=True)
    os.makedirs(os.path.join(path_green, "GreenSpec"), exist_ok=True)
    path_func = os.path.join(
        path_green, "GreenFunc", "%.2f" % event_depth, "%.2f" % receiver_depth, ""
    )
    path_spec = os.path.join(
        path_green, "GreenSpec", "%.2f" % event_depth, "%.2f" % receiver_depth, ""
    )
    os.makedirs(path_func, exist_ok=True)
    os.makedirs(path_spec, exist_ok=True)
    return path_func, path_spec


def create_inp_spgrn2020(
    path_green,
    event_depth,
    receiver_depth,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    source_radius,
    cal_gf,
    time_window,
    green_before_p,
    source_duration,
    dist_range,
    delta_dist_range,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    path_func = str(
        os.path.join(
            path_green, "GreenFunc", "%.2f" % event_depth, "%.2f" % receiver_depth, ""
        )
    )
    path_spec = str(
        os.path.join(
            path_green, "GreenSpec", "%.2f" % event_depth, "%.2f" % receiver_depth, ""
        )
    )

    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]
    last_line = [lines[-1]]
    lines_earth = lines[113:-1]
    lines = lines[:113]  # cutoff earth model
    lines[25] = "%.2f\n" % receiver_depth
    lines[41] = "%f  %f\n" % (spec_time_window, sampling_interval)
    lines[42] = "%f\n" % max_frequency
    lines[43] = "%f\n" % max_slowness
    lines[44] = "%f\n" % anti_alias
    lines[52] = "%f %d\n" % (gravity_fc, gravity_harmonic)
    lines[60] = "%d %d\n" % (cal_sph, cal_tor)
    lines[71] = '"%s"\n' % path_spec
    lines[73] = '%.2f  %.2f  "grn_d%.2f"  %d\n' % (
        event_depth,
        source_radius,
        event_depth,
        cal_gf,
    )
    lines[90] = '"%s"\n' % path_func
    lines[91] = '"GreenInfo%.2f.dat"\n' % event_depth
    lines[94] = "%f  %f\n" % (time_window, sampling_interval)
    lines[95] = "%f\n" % -green_before_p
    lines[96] = "%f\n" % source_duration
    lines[98] = "%f  %f  %f  %f\n" % (
        dist_range[0],
        dist_range[1],
        delta_dist_range[0],
        delta_dist_range[1],
    )
    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    lines[106] = "%d  %d\n" % (earth_model_layer_num, physical_dispersion)
    path_inp = os.path.join(path_func, "grn.inp")
    with open(path_inp, "w") as fw:
        fw.writelines(lines + lines_earth + last_line)
    return path_inp


def call_spgrn2020(event_depth, receiver_depth, path_green, check_finished=False):
    # print(event_depth, receiver_depth, path_green, check_finished)
    sub_sub_dir = str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
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
        with open(path_finished, "r", encoding="utf-8") as fr:
            output = fr.readlines()
        return output

    output = call_exe(
        path_inp=path_inp,
        path_finished=path_finished,
        name="spgrn2020",
    )
    return output


def check_output_spgrn(path_func, event_depth, tables=(), expected_info=None,
                       check_values=False):
    """List the problems with the output of one SPGRN job.

    Parameters
    ----------
    path_func : str
        Job directory, path_green/GreenFunc/event_depth/receiver_depth.
    event_depth : float
        Source depth in km, used in the output file names.
    tables : sequence of str, optional
        Names of Fortran-record travel-time tables the job writes (SPGRN2020:
        tptable.dat, tstable.dat); each holds 12 bytes per distance plus two
        4-byte record markers. Default: ().
    expected_info : dict or None, optional
        Library metadata whose dist_list and samples_num the job must match.
        Default: None.
    check_values : bool, optional
        Also read the Green's functions and report NaN or infinite values. Default: False.

    Returns
    -------
    problems : list of str
        Empty when GreenInfo is readable and grn_d holds every distance.
    """
    path_info = os.path.join(path_func, "GreenInfo%.2f.dat" % event_depth)
    if not os.path.exists(path_info):
        return ["%s is missing" % path_info]
    try:
        info = read_green_info_spgrn(path_func, event_depth)
    except (ValueError, IndexError) as exc:
        return ["%s cannot be read: %s" % (path_info, exc)]
    n_dist = len(info["dist_list"])
    # per distance: 3 header values, then 10 records of samples between markers
    n_each = 3 + (2 + info["samples_num"]) * 10
    path_grn = os.path.join(path_func, "grn_d%.2f" % event_depth)
    problems = []
    problem = check_file_size(path_grn, n_dist * n_each * 4)
    if problem:
        problems.append(problem)
    elif check_values:
        data = np.fromfile(path_grn, dtype=np.float32).reshape(n_dist, n_each)
        samples = data[:, 3:].reshape(n_dist, 10, -1)[:, :, 1:-1]
        if not np.all(np.isfinite(samples)):
            problems.append("%s contains NaN or infinite values" % path_grn)
    for name in tables:
        problem = check_file_size(os.path.join(path_func, name), 8 + 12 * n_dist)
        problems += [problem] if problem else []
    if expected_info is not None and (
        info["dist_list"] != expected_info["dist_list"]
        or info["samples_num"] != expected_info["samples_num"]
    ):
        problems.append(
            "%s has other distances or samples than green_lib_info.json" % path_info
        )
    return problems


if __name__ == "__main__":
    pass
