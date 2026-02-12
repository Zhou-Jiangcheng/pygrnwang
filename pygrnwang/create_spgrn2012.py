import os

from .spgrn2012inp import s as str_inp
from .utils import call_exe, convert_earth_model_nd2inp


def create_inp_spgrn2012(
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
        t0,
        v0,
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
    lines_earth = lines[110:-1]
    lines = lines[:110]  # cutoff earth model
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
    lines[89] = '"%s"\n' % path_func
    lines[90] = '"GreenInfo%.2f.dat"\n' % event_depth
    lines[91] = "%f  %f\n" % (time_window, sampling_interval)
    lines[92] = "%f %f\n" % (t0, v0)
    lines[93] = "%f\n" % source_duration
    lines[95] = "%f  %f  %f  %f\n" % (
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
    lines[103] = "%d  %d\n" % (earth_model_layer_num, physical_dispersion)
    path_inp = os.path.join(path_func, "grn.inp")
    with open(path_inp, "w") as fw:
        fw.writelines(lines + lines_earth + last_line)
    return path_inp


def call_spgrn2012(event_depth, receiver_depth, path_green, check_finished=False):
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
        with open(path_finished, "r") as fr:
            output = fr.readlines()
        return output

    output = call_exe(
        path_inp=path_inp,
        path_finished=path_finished,
        name="spgrn2012",
    )


if __name__ == "__main__":
    pass
