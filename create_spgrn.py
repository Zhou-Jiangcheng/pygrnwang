import os
import platform
import shutil
import subprocess
import datetime

from pygrnwang.utils import convert_earth_model_nd2inp

d2km = 111.19492664455873


def create_dir(event_depth, path_green):
    # 新建文件夹并返回路径,若文件夹已存在则直接返回路径
    if not os.path.exists(os.path.join(path_green, "GreenFunc")):
        os.mkdir(os.path.join(path_green, "GreenFunc"))
    if not os.path.exists(os.path.join(path_green, "GreenSpec")):
        os.mkdir(os.path.join(path_green, "GreenSpec"))
    path_func = os.path.join(path_green, "GreenFunc", "%.1f" % event_depth, "")
    path_spec = os.path.join(path_green, "GreenSpec", "%.1f" % event_depth, "")
    if not os.path.exists(path_func):
        os.mkdir(path_func)
    if not os.path.exists(path_spec):
        os.mkdir(path_spec)
    return path_func, path_spec


def create_inp(
    event_depth,
    path_green,
    spec_time_window,
    time_window,
    sampling_interval,
    before_tp,
    dist_range,
    delta_dist_range,
    max_frequency,
    path_nd=None,
    earth_model_layer_num=None,
):
    path_func = str(os.path.join(path_green, "GreenFunc", "%.1f" % event_depth, ""))
    path_spec = str(os.path.join(path_green, "GreenSpec", "%.1f" % event_depth, ""))
    path_inp = os.path.join(path_green, "spgrn2020.inp")
    if os.path.exists(path_inp):
        with open(path_inp, "r") as fr:
            lines = fr.readlines()
    else:
        from pygrnwang.spgrn2020inp import s

        lines = s.split("\n")
        lines = [line + "\n" for line in lines]
    last_line = [lines[-1]]
    lines_earth = lines[113:-1]
    lines = lines[:113]  # cutoff earth model
    lines[41] = "%.2f  %.2f\n" % (spec_time_window, sampling_interval)
    lines[42] = "%.2f\n" % max_frequency
    lines[71] = '"%s"\n' % path_spec
    lines[73] = '%.1f  0.00  "grn_d%.1f"  1\n' % (event_depth, event_depth)
    lines[90] = '"%s"\n' % path_func
    lines[91] = '"GreenInfo%.1f.dat"\n' % event_depth
    lines[94] = "%.2f  %.2f\n" % (time_window, sampling_interval)
    lines[95] = "%.1f\n" % -before_tp
    lines[98] = "%.1f  %.1f  %.1f  %.1f\n" % (
        dist_range[0] * d2km,
        dist_range[1] * d2km,
        delta_dist_range[0],
        delta_dist_range[1],
    )
    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    lines[106] = "%d  0\n" % earth_model_layer_num
    path_inp = os.path.join(path_func, "%.1f.inp" % event_depth)
    with open(path_inp, "w") as fw:
        fw.writelines(lines + lines_earth + last_line)
    return path_inp


def call_spgrn(event_depth, path_green):
    os.chdir(path_green)
    path_inp = str(
        os.path.join(
            path_green, "GreenFunc", "%.1f" % event_depth, "%.1f.inp" % event_depth
        )
    )

    if platform.system() == "Windows":
        spgrn_process = subprocess.Popen(
            [os.path.join(path_green, "spgrn2020.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE)
        spgrn_process.communicate(str.encode(path_inp))
    else:
        try:
            spgrn_process = subprocess.Popen(
                [os.path.join(path_green, "spgrn2020.bin")],
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
            spgrn_process.communicate(str.encode(path_inp))
        except Exception as e:
            print(e)
            raise ("this system is not supported yet, \
                please compile the source code of spgrn2020, \
                copy and replace the spgrn2020.bin file ")


def create_grnlib(
    event_depth,
    path_green,
    path_bin,
    spec_time_window,
    time_window,
    sampling_interval,
    before_tp,
    dist_range,
    delta_dist_range,
    max_frequency,
    path_nd=None,
    earth_model_layer_num=None,
):
    print("creating green func lib recv_depth=%d" % event_depth)
    s = datetime.datetime.now()
    create_dir(event_depth, path_green)
    create_inp(
        event_depth,
        path_green,
        spec_time_window,
        time_window,
        sampling_interval,
        before_tp,
        dist_range,
        delta_dist_range,
        max_frequency,
        path_nd,
        earth_model_layer_num,
    )
    path_bin_call = os.path.join(path_green, "spgrn2020.bin")
    if not os.path.exists(path_bin_call):
        shutil.copy(path_bin, path_bin_call)
    call_spgrn(event_depth, path_green)
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))
    print("done")


if __name__ == "__main__":
    pass
