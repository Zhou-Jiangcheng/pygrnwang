import os
import json

import numpy as np
import scipy.signal as signal
import pandas as pd

from .utils import shift_green2real_tpts
from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz
from .pytaup import cal_first_p_s
from .signal_process import resample


def read_time_series_qseis06_bin(path_greenfunc, start_count):
    time_series_list = []
    for com in ["tr", "tz"]:
        time_series_com = np.load(os.path.join(path_greenfunc, "grn_%s.npy" % com))
        time_series_list.append(time_series_com[:, start_count].reshape(4, -1))
    for com in ["tt"]:
        time_series_com = np.load(os.path.join(path_greenfunc, "grn_%s.npy" % com))
        time_series_list.append(time_series_com[:, start_count].reshape(2, -1))
    return time_series_list


def read_time_series_qseis06_ascii(path_greenfunc, start_count):
    time_seris_list = []
    start_count = start_count + 1
    for com in ["tr", "tz"]:
        ex_com = pd.read_csv(
            os.path.join(path_greenfunc, "ex.%s" % com), sep="\\s+"
        ).to_numpy()  # type:ignore
        ss_com = pd.read_csv(
            os.path.join(path_greenfunc, "ss.%s" % com), sep="\\s+"
        ).to_numpy()  # type:ignore
        ds_com = pd.read_csv(
            os.path.join(path_greenfunc, "ds.%s" % com), sep="\\s+"
        ).to_numpy()  # type:ignore
        cl_com = pd.read_csv(
            os.path.join(path_greenfunc, "cl.%s" % com), sep="\\s+"
        ).to_numpy()  # type:ignore
        time_series_com = np.concatenate(
            [
                ex_com[:, start_count],
                ss_com[:, start_count],
                ds_com[:, start_count],
                cl_com[:, start_count],
            ]
        ).reshape(4, -1)
        time_series_com = np.array(time_series_com, dtype=np.float32)
        time_seris_list.append(time_series_com)

    for com in ["tt"]:
        ss_r = pd.read_csv(
            os.path.join(path_greenfunc, "ss.%s" % com), sep="\\s+"
        ).to_numpy()  # type:ignore
        ds_r = pd.read_csv(
            os.path.join(path_greenfunc, "ds.%s" % com), sep="\\s+"
        ).to_numpy()  # type:ignore
        time_series_com = np.concatenate(
            [
                ss_r[:, start_count],
                ds_r[:, start_count],
            ]
        ).reshape(2, -1)
        time_series_com = np.array(time_series_com, dtype=np.float32)
        time_seris_list.append(time_series_com)
    return time_seris_list


def synthesize_qseis06(time_series_list, m1, m2):
    r = (
        time_series_list[0][0] * m1[0]
        + time_series_list[0][1] * m1[1]
        + time_series_list[0][2] * m1[2]
        + time_series_list[0][3] * m1[3]
    )
    z = (
        time_series_list[1][0] * m1[0]
        + time_series_list[1][1] * m1[1]
        + time_series_list[1][2] * m1[2]
        + time_series_list[1][3] * m1[3]
    )
    t = time_series_list[2][0] * m2[0] + time_series_list[2][1] * m2[1]
    return r, t, z


def seek_qseis06(
    path_green,
    event_depth_km,
    receiver_depth_km,
    az_deg,
    dist_km,
    focal_mechanism,
    srate,
    output_type="disp",
    rotate=True,
    before_p=None,
    pad_zeros=False,
    shift=False,
    only_seismograms=True,
    model_name="ak135fc",
    green_info=None,
):
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    srate_grn = 1 / green_info["sampling_interval"]
    sampling_num = green_info["sampling_num"]
    time_reduction_velo = green_info["time_reduction_velo"]
    dist_range = green_info["dist_range"]
    delta_dist = green_info["delta_dist"]
    num_each_group = green_info["N_each_group"]
    wavelet_type = green_info["wavelet_type"]
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

    path_greenfunc = str(
        os.path.join(path_green, "%.2f" % grn_dep_source, "%.2f" % grn_dep_receiver)
    )
    ind = round((dist_km - dist_range[0]) / delta_dist)
    ind_group = ind // num_each_group
    green_dist = dist_range[0] + ind * delta_dist
    start_count = ind - ind_group * num_each_group
    if time_reduction_velo != 0:
        time_reduction = green_dist / time_reduction_velo
    else:
        time_reduction = 0
    path_greenfunc_sub = os.path.join(path_greenfunc, "%d_0" % ind_group)

    if os.path.exists(os.path.join(path_greenfunc_sub, "grn_tz.npy")):
        time_series_list = read_time_series_qseis06_bin(
            path_greenfunc=path_greenfunc_sub, start_count=start_count
        )
    else:
        time_series_list = read_time_series_qseis06_ascii(
            path_greenfunc=path_greenfunc_sub, start_count=start_count
        )

    # r,t,z
    [M11, M12, M13, M22, M23, M33] = check_convert_fm(focal_mechanism=focal_mechanism)

    exp = (M11 + M22 + M33) / 3
    ss1 = M12
    ss2 = (M11 - M22) / 2
    ds1 = M13
    ds2 = M23
    clvd = (-M11 / 2 - M22 / 2 + M33) / (3 / 2)

    az_rad = np.deg2rad(az_deg)
    sin_az, cos_az = np.sin(az_rad), np.cos(az_rad)
    sin_2az, cos_2az = np.sin(2 * az_rad), np.cos(2 * az_rad)

    m1 = [exp, ss1 * sin_2az + ss2 * cos_2az, ds1 * cos_az + ds2 * sin_az, clvd]
    m2 = [ss1 * cos_2az - ss2 * sin_2az, ds1 * sin_az - ds2 * cos_az]

    r, t, z = synthesize_qseis06(time_series_list=time_series_list, m1=m1, m2=m2)

    if rotate:
        seismograms = rotate_rtz_to_enz(az_in_deg=az_deg, r=r, t=-t, z=-z)
    else:
        seismograms = np.array([r, -t, -z])

    tpts_table = None
    if (before_p is not None) or shift:
        first_p_grn, first_s_grn = cal_first_p_s(
            event_depth_km=event_depth_km,
            receiver_depth_km=receiver_depth_km,
            dist_km=green_dist,
            model_name=model_name,
        )
        tpts_table = {"p_onset": first_p_grn, "s_onset": first_s_grn}

    ts_count = 0
    if before_p is not None:
        ts_count = round(
            (tpts_table["p_onset"] - time_reduction - before_p) * srate_grn
        )
    if pad_zeros:
        if before_p is not None:
            raise ValueError("can not set before_p and pad_zeros together")
        ts_count = round(-time_reduction * srate_grn)
    seismograms = np.roll(seismograms, -ts_count)
    if ts_count > 0:
        seismograms[:, -ts_count:] = 0
    elif ts_count < 0:
        seismograms[:, :-ts_count] = 0

    first_p = None
    first_s = None
    if shift:
        seismograms, first_p, first_s = shift_green2real_tpts(
            seismograms=seismograms,
            tpts_table=tpts_table,
            srate=srate_grn,
            green_before_p=tpts_table["p_onset"] - time_reduction,
            event_depth_km=event_depth_km,
            dist_in_km=dist_km,
            receiver_depth_km=receiver_depth_km,
            model_name=model_name,
        )

    # conv_shift = round(green_info["wavelet_duration"] / 2)
    # if conv_shift != 0:
    #     seismograms = np.roll(seismograms, -conv_shift)
    #     seismograms[:, -conv_shift:] = 0

    seismograms_resample = np.zeros((3, round(sampling_num * srate / srate_grn)))
    for i in range(3):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
        )

    if wavelet_type == 1 and output_type == "disp":
        seismograms_resample = np.cumsum(seismograms_resample, axis=1) / srate
    elif wavelet_type == 2 and output_type == "velo":
        seismograms_resample = (
            signal.convolve(
                seismograms_resample.T,
                np.array([1, -1])[:, None],
                mode="same",
                method="auto",
            ).T
            / srate
        )
    if only_seismograms:
        return seismograms_resample
    else:
        return (
            seismograms_resample,
            tpts_table,
            first_p,
            first_s,
            grn_dep_source,
            grn_dep_receiver,
            green_dist,
        )


if __name__ == "__main__":
    pass
