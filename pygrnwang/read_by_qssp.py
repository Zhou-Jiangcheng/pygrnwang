import os
import hashlib
import subprocess
import platform
import json

import numpy as np
import pandas as pd

from .qssp2020inp import s as str_inp
from .signal_process import resample
from .utils import convert_earth_model_nd2inp
from .focal_mechanism import check_convert_fm, convert_mt_axis


def hash_read_pars(
    event_lat,
    event_lon,
    event_depth,
    receiver_lat,
    receiver_lon,
    receiver_depth,
):
    data = (
        str(float(event_lat)).encode()
        + str(float(event_lon)).encode()
        + str(float(event_depth)).encode()
        + str(float(receiver_lat)).encode()
        + str(float(receiver_lon)).encode()
        + str(float(receiver_depth)).encode()
    )
    return hashlib.md5(data).hexdigest()


def create_inp_qssp2020_read(
    read_name,
    path_green,
    event_lat,
    event_lon,
    event_depth,
    receiver_lat,
    receiver_lon,
    receiver_depth,
    focal_mechanism,
    green_info,
):
    spec_time_window = green_info["spec_time_window"]
    sampling_interval = green_info["sampling_interval"]
    max_frequency = green_info["max_frequency"]
    max_slowness = green_info["max_slowness"]
    anti_alias = green_info["anti_alias"]
    turning_point_filter = green_info["turning_point_filter"]
    turning_point_d1 = green_info["turning_point_d1"]
    turning_point_d2 = green_info["turning_point_d2"]
    free_surface_filter = green_info["free_surface_filter"]
    gravity_fc = green_info["gravity_fc"]
    gravity_harmonic = green_info["gravity_harmonic"]
    cal_sph = green_info["cal_sph"]
    cal_tor = green_info["cal_tor"]
    min_harmonic = green_info["min_harmonic"]
    max_harmonic = green_info["max_harmonic"]
    source_radius = green_info["source_radius"]
    source_duration = green_info["source_duration"]
    output_observables = green_info["output_observables"]
    time_window = green_info["time_window"]
    time_reduction = green_info["time_reduction"]
    path_nd = green_info["path_nd"]
    earth_model_layer_num = green_info["earth_model_layer_num"]
    physical_dispersion = green_info["physical_dispersion"]

    path_spec = str(
        os.path.join(
            path_green,
            "GreenSpec",
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            "",
        )
    )
    if read_name == "hash":
        read_name = hash_read_pars(
            event_lat,
            event_lon,
            event_depth,
            receiver_lat,
            receiver_lon,
            receiver_depth,
        )
    path_func = str(os.path.join(path_green, "read", read_name, ""))
    os.makedirs(path_func, exist_ok=True)

    lines = str_inp.split("\n")
    lines = [line + "\n" for line in lines]
    last_line = [lines[-1]]
    lines_earth_head = lines[141:155]
    lines_earth = lines[155:-1]
    lines = lines[:140]  # cutoff locs

    lines[24] = "%.2f\n" % receiver_depth

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

    lines[75] = "%d %d %d %d\n" % (cal_sph, cal_tor, min_harmonic, max_harmonic)

    lines[86] = "1 %f '%s'\n" % (source_radius, path_spec)
    lines[87] = "%.2f 'Green_%.2fkm' %d\n" % (event_depth, event_depth, 0)

    lines[111] = "1 1\n"
    mt_ned = check_convert_fm(focal_mechanism)
    mt = convert_mt_axis(mt=mt_ned, convert_flag="ned2rtp")
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

    lines_receiver = [
        "1\n",
        "%.4f %.4f 'read' %.2f\n" % (receiver_lat, receiver_lon, time_reduction),
    ]

    lines = lines + lines_receiver
    path_inp = os.path.join(path_func, "read.inp")

    if path_nd is not None:
        lines_earth = convert_earth_model_nd2inp(
            path_nd=path_nd, path_output="earth_model.dat"
        )
    if earth_model_layer_num is None:
        earth_model_layer_num = len(lines_earth)
    lines_earth = lines_earth[:earth_model_layer_num]
    lines_earth_head[7] = "%d  %d\n" % (earth_model_layer_num, physical_dispersion)

    with open(path_inp, "w") as fw:
        fw.writelines(lines + lines_earth_head + lines_earth + last_line)
    return path_inp


def call_qssp2020_read(path_green, path_inp, check_finished=False):
    os.chdir(path_green)
    sub_dir = os.path.dirname(path_inp)
    hash_v = os.path.basename(sub_dir)
    path_inp = os.path.join(".", "read", hash_v, "read.inp").replace("'", "")
    path_finished = os.path.join(sub_dir, ".finished")
    if check_finished and os.path.exists(path_finished):
        return None

    if platform.system() == "Windows":
        spgrn_process = subprocess.Popen(
            [os.path.join(path_green, "qssp2020.exe")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        spgrn_process.communicate(str.encode(path_inp))
    else:
        spgrn_process = subprocess.Popen(
            [os.path.join(path_green, "qssp2020.bin")],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
        )
        spgrn_process.communicate(str.encode(path_inp))

    with open(path_finished, "w") as fw:
        fw.writelines([])


def read_by_qssp(
    path_green,
    event_lat,
    event_lon,
    event_depth_km,
    receiver_lat,
    receiver_lon,
    receiver_depth_km,
    focal_mechanism,
    srate,
    read_name="hash",
    output_type="disp",
    green_info=None,
    check_finished=False,
):
    """ """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    srate_grn = 1 / green_info["sampling_interval"]
    sampling_num = (
        round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    )
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]
    if not isinstance(grn_dep_list, list):
        grn_dep_source = grn_dep_list
    else:
        grn_dep_source = grn_dep_list[
            np.argmin(np.abs(event_depth_km - np.array(grn_dep_list)))
        ]
    if not isinstance(grn_receiver_list, list):
        grn_dep_receiver = grn_receiver_list
    else:
        grn_dep_receiver = grn_receiver_list[
            np.argmin(np.abs(receiver_depth_km - np.array(grn_receiver_list)))
        ]

    path_inp = create_inp_qssp2020_read(
        read_name=read_name,
        path_green=path_green,
        event_lat=event_lat,
        event_lon=event_lon,
        event_depth=grn_dep_source,
        receiver_lat=receiver_lat,
        receiver_lon=receiver_lon,
        receiver_depth=grn_dep_receiver,
        focal_mechanism=focal_mechanism,
        green_info=green_info,
    )
    call_qssp2020_read(
        path_green=path_green, path_inp=path_inp, check_finished=check_finished
    )

    if output_type == "gravimeter":
        enz_list = [""]
    elif output_type in ["disp", "velo", "acce", "rota", "rota_rate", "gravitation"]:
        enz_list = ["_e", "_n", "_z"]
    elif output_type in ["stress", "stress_rate", "strain", "strain_rate"]:
        enz_list = ["_ee", "_en", "_ez", "_nn", "_nz", "_zz"]
    else:
        raise ValueError(
            "output_type must in  disp | velo | acce | strain | strain_rate | "
            "stress | stress_rate | rotation | rotation_rate | gravitation | gravimeter"
        )

    path_func = os.path.dirname(path_inp)
    data_enz = np.zeros((sampling_num, len(enz_list)))
    for i_enz in range(len(enz_list)):
        df = pd.read_csv(
            os.path.join(path_func, "_%s%s.dat" % (output_type, enz_list[i_enz])),
            sep="\\s+",
        )
        data_enz[:, i_enz] = df["read"].values

    seismograms = data_enz.T

    conv_shift = round(green_info["source_duration"] * srate_grn / 2)
    if conv_shift != 0:
        seismograms = np.roll(seismograms, -conv_shift)
        seismograms[:, -conv_shift:] = 0

    seismograms_resample = np.zeros(
        (seismograms.shape[0], round(sampling_num * srate / srate_grn))
    )
    for i in range(seismograms.shape[0]):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
        )

    return seismograms_resample
