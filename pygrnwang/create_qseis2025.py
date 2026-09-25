import os
import math

import numpy as np
import pandas as pd

from .utils import (
    convert_earth_model_nd2inp,
    call_exe,
    check_ascii_table,
    check_file_size,
    write_bin_atomic,
)
from .qseis2025inp import s as str_inp

# nrmax in fortran_src_codes/qseis2025_src/qsglobal.h: distances per input file
QSEIS2025_NRMAX = 101
# The Fortran arrays are static, so every qseis2025 process commits about
# 1.3 GiB (measured) whatever the grid size; the margin covers the Python side.
QSEIS2025_MEMORY_PER_JOB_GB = 1.4

# output files for each of the five output_observables flags
# (disp/velo, volume, strain, stress, rotation):
# P-SV components for ex, ss, ds, cl and SH components for ss, ds
PSV_COMS_QSEIS2025 = [["tz", "tr"], ["tv"], ["ezz", "ezr", "err", "ett"],
                      ["szz", "szr", "srr", "stt"], ["ot"]]
SH_COMS_QSEIS2025 = [["tt"], [], ["ezt", "ert"], ["szt", "srt"], ["oz", "or"]]
PSV_STYPES = ["ex", "ss", "ds", "cl"]
SH_STYPES = ["ss", "ds"]


def create_dir_qseis2025(
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


def create_inp_qseis2025(
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
    output_observables,
    slowness_int_algorithm=0,
    eps_estimate_wavenumber=1e-6,
    source_radius_ratio=0.05,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=12,
    anti_alias=0.01,
    free_surface=0,
    wavelet_duration=4,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
):
    dist_range = dist_range.copy()
    output_observables = output_observables.copy()
    path_sub_dir = str(
        os.path.join(path_green, "%.2f" % event_depth, "%.2f" % receiver_depth)
    )
    # when receiver depth is not 0, change dists in inp file
    if flat_earth_transform:
        r_ratio = (6371 - receiver_depth) / 6371
    else:
        r_ratio = 1
    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]

    lines_earth = lines[233:-22]
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
    lines[73] = "%d\n" % slowness_int_algorithm
    lines[74] = "%f %f\n" % (eps_estimate_wavenumber, source_radius_ratio)
    if slowness_window is not None:
        lines[75] = "%f %f %f %f\n" % (
            slowness_window[0],
            slowness_window[1],
            slowness_window[2],
            slowness_window[3],
        )
    else:
        lines[75] = "0.0 0.0 0.0 0.0\n"
    lines[76] = "%f\n" % wavenumber_sampling_rate
    lines[77] = "%f\n" % anti_alias

    # OPTIONS FOR PARTIAL SOLUTIONS
    lines[111] = "%d\n" % free_surface

    # SOURCE TIME FUNCTION (WAVELET) PARAMETERS (Note 3)
    lines[132] = "%d %d\n" % (wavelet_duration, wavelet_type)

    # OUTPUT FILES FOR GREEN'S FUNCTIONS (Note 4)
    lines[181] = " ".join(["%d" % output_observables[_] for _ in range(5)]) + "\n"

    # GLOBAL MODEL PARAMETERS (Note 5)
    if flat_earth_transform:
        lines[217] = "1\n"
    else:
        lines[217] = "0\n"

    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    else:
        lines_earth = lines_earth[:earth_model_layer_num]
    lines[226] = "%d\n" % earth_model_layer_num
    lines = lines[:227] + lines_earth + lines_end

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
            # keep the spacing at delta_dist: the last point may pass dist_range[1]
            (dist_range[0] + (N_dist - 1) * delta_dist) * r_ratio,
        )
        path_inp = os.path.join(path_sub_dir, "%d_0" % (N_dist_group - 1), "grn.inp")
        with open(path_inp, "w") as fw:
            fw.writelines(lines)
    return path_inp


def call_qseis2025(
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
        name="qseis2025",
    )
    return output


def read_output_observables_qseis2025(path_inp):
    """Read the five output_observables flags from a generated grn.inp."""
    with open(path_inp, "r") as fr:
        lines = fr.readlines()
    # create_inp_qseis2025 writes the flags on this line
    flags = lines[181].split() if len(lines) > 181 else []
    if len(flags) != 5 or any(flag not in ("0", "1") for flag in flags):
        raise ValueError("Cannot read output_observables from %s" % path_inp)
    return [int(flag) for flag in flags]


def qseis2025_output_coms(output_observables):
    """Return (component, source types) for the files of the selected observables."""
    coms = []
    for ind, selected in enumerate(output_observables):
        if selected:
            coms += [(com, PSV_STYPES) for com in PSV_COMS_QSEIS2025[ind]]
            coms += [(com, SH_STYPES) for com in SH_COMS_QSEIS2025[ind]]
    return coms


def check_output_qseis(path_greenfunc, coms, sampling_num, n_dist, check_values=False):
    """List the problems with the output of one QSEIS job.

    Parameters
    ----------
    path_greenfunc : str
        Job directory.
    coms : list of tuple
        (component, source types) of every file the job must write.
    sampling_num : int
        Samples per trace.
    n_dist : int
        Distances computed by the job.
    check_values : bool, optional
        Also read binary files and report NaN or infinite values. Default: False.

    Returns
    -------
    problems : list of str
        Empty when every component exists, as a binary file of the size the
        readers expect or as complete ASCII files.
    """
    problems = []
    for com, stypes in coms:
        path_bin = os.path.join(path_greenfunc, "grn_%s.bin" % com)
        if os.path.exists(path_bin):
            problem = check_file_size(
                path_bin, len(stypes) * sampling_num * n_dist * 4, check_values
            )
            problems += [problem] if problem else []
            continue
        for stype in stypes:
            # a header line, then one row per sample: time and one value per distance
            problem = check_ascii_table(
                os.path.join(path_greenfunc, "%s.%s" % (stype, com)),
                sampling_num,
                n_dist + 1,
            )
            problems += [problem] if problem else []
    return problems


def check_output_qseis2025(
    path_greenfunc, output_observables, sampling_num, n_dist, check_values=False
):
    """List the problems with the output of one qseis2025 job (see check_output_qseis)."""
    return check_output_qseis(
        path_greenfunc,
        qseis2025_output_coms(output_observables),
        sampling_num,
        n_dist,
        check_values,
    )


def convert_qseis_ascii(path_greenfunc, coms, remove=False):
    """Convert the ASCII output of one QSEIS job to grn_<component>.bin files.

    Components without any ASCII file are skipped; ASCII files are removed
    only after their binary file is written.
    """
    for com, stypes in coms:
        paths_ascii = [
            os.path.join(path_greenfunc, "%s.%s" % (stype, com)) for stype in stypes
        ]
        exist = [os.path.exists(path) for path in paths_ascii]
        if not any(exist):
            continue
        if not all(exist):
            raise ValueError(
                "Incomplete output in %s: %s missing"
                % (
                    path_greenfunc,
                    ", ".join(
                        os.path.basename(path)
                        for path, ok in zip(paths_ascii, exist)
                        if not ok
                    ),
                )
            )
        time_series_com = [
            pd.read_csv(path, sep="\\s+").to_numpy()[:, 1:] for path in paths_ascii
        ]
        if len({data.shape for data in time_series_com}) != 1:
            raise ValueError(
                "Output files for %s in %s have different sizes; "
                "the job did not finish" % (com, path_greenfunc)
            )
        write_bin_atomic(
            np.concatenate(time_series_com, dtype=np.float32).T,
            os.path.join(path_greenfunc, "grn_%s.bin" % com),
        )
        if remove:
            for path in paths_ascii:
                os.remove(path)


def convert_pd2bin_qseis2025(path_greenfunc, remove=False):
    convert_qseis_ascii(path_greenfunc, qseis2025_output_coms([1] * 5), remove)


if __name__ == "__main__":
    pass
