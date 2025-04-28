import os
import subprocess
import platform

from .utils import (
    convert_earth_model_nd2inp
)


def create_dir_spgrn2020(event_depth, receiver_depth, path_green):
    # 新建文件夹并返回路径,若文件夹已存在则直接返回路径
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
    path_inp = os.path.join(path_green, "spgrn2020.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from .spgrn2020inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]
    last_line = [lines[-1]]
    lines_earth = lines[113:-1]
    lines = lines[:113]  # cutoff earth model
    lines[25] = "%.2f\n" % receiver_depth
    lines[41] = "%.2f  %.2f\n" % (spec_time_window, sampling_interval)
    lines[42] = "%.2f\n" % max_frequency
    lines[43] = "%.2f\n" % max_slowness
    lines[44] = "%.2f\n" % anti_alias
    lines[52] = "%.2f %d\n" % (gravity_fc, gravity_harmonic)
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
    lines[94] = "%.2f  %.2f\n" % (time_window, sampling_interval)
    lines[95] = "%.2f\n" % -green_before_p
    lines[96] = "%.2f\n" % source_duration
    lines[98] = "%.2f  %.2f  %.2f  %.2f\n" % (
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
    # os.chdir(path_green)
    path_inp = str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            "grn.inp",
        )
    )

    sub_dir = os.path.dirname(path_inp)
    path_finished = os.path.join(sub_dir, ".finished")
    if check_finished and os.path.exists(path_finished):
        return None

    if platform.system() == "Windows":
        spgrn_process = subprocess.Popen(
            [os.path.join(path_green, "spgrn2020.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        spgrn_process.communicate(str.encode(path_inp))
    else:
        spgrn_process = subprocess.Popen(
            [os.path.join(path_green, "spgrn2020.bin")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        spgrn_process.communicate(str.encode(path_inp))

    with open(path_finished, "w") as fw:
        fw.writelines([])


if __name__ == "__main__":
    pass
