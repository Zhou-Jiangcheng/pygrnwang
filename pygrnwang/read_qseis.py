import os
import json

import numpy as np
import pandas as pd

from .utils import shift_green2real_tpts
from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz, cal_first_p_s
from .signal_process import resample


def read_time_series_qseis06_bin(path_bin, start_count, sampling_num):
    fr = open(path_bin, "rb")
    time_series = np.fromfile(
        file=fr,
        dtype=np.float32,
        count=sampling_num * 10,
        offset=start_count * sampling_num * 10 * 4,
    )
    time_series = time_series.reshape(10, sampling_num)
    fr.close()
    return time_series


def read_time_series_qseis06_ascii(path_greenfunc, start_count):
    start_count = start_count + 1
    ex_z_df = pd.read_csv(os.path.join(path_greenfunc, "ex.tz"), sep="\\s+")
    ex_r_df = pd.read_csv(os.path.join(path_greenfunc, "ex.tr"), sep="\\s+")
    ss_z_df = pd.read_csv(os.path.join(path_greenfunc, "ss.tz"), sep="\\s+")
    ss_r_df = pd.read_csv(os.path.join(path_greenfunc, "ss.tr"), sep="\\s+")
    ss_t_df = pd.read_csv(os.path.join(path_greenfunc, "ss.tt"), sep="\\s+")
    ds_z_df = pd.read_csv(os.path.join(path_greenfunc, "ds.tz"), sep="\\s+")
    ds_r_df = pd.read_csv(os.path.join(path_greenfunc, "ds.tr"), sep="\\s+")
    ds_t_df = pd.read_csv(os.path.join(path_greenfunc, "ds.tt"), sep="\\s+")
    cl_z_df = pd.read_csv(os.path.join(path_greenfunc, "cl.tz"), sep="\\s+")
    cl_r_df = pd.read_csv(os.path.join(path_greenfunc, "cl.tr"), sep="\\s+")
    time_series = np.concatenate(
        [
            ex_z_df.iloc[:, start_count].values,
            ex_r_df.iloc[:, start_count].values,
            ss_z_df.iloc[:, start_count].values,
            ss_r_df.iloc[:, start_count].values,
            ss_t_df.iloc[:, start_count].values,
            ds_z_df.iloc[:, start_count].values,
            ds_r_df.iloc[:, start_count].values,
            ds_t_df.iloc[:, start_count].values,
            cl_z_df.iloc[:, start_count].values,
            cl_r_df.iloc[:, start_count].values,
        ]
    ).T
    time_series = time_series.reshape(10, len(ex_z_df))
    return time_series


def synthesize_qseis06(az_in_deg, time_series, focal_mechanism):
    """
    Note 4:
    Double-Couple   m11/ m22/ m33/ m12/ m23/ m31  Azimuth_Factor_(tz,tr,tv)/(tt)
    ============================================================================
    explosion       1.0/ 1.0/ 1.0/ -- / -- / --       1.0         /   0.0
    strike-slip     -- / -- / -- / 1.0/ -- / --       sin(2*azi)  /   cos(2*azi)
                    1.0/-1.0/ -- / -- / -- / --       cos(2*azi)  /  -sin(2*azi)
    dip-slip        -- / -- / -- / -- / -- / 1.0      cos(azi)    /   sin(azi)
                    -- / -- / -- / -- / 1.0/ --       sin(azi)    /  -cos(azi)
    clvd           -0.5/-0.5/ 1.0/ -- / -- / --       1.0         /   0.0
    ============================================================================
    Single-Force    fx / fy / fz                  Azimuth_Factor_(tz,tr,tv)/(tt)
    ============================================================================
    fz              -- / -- / 1.0                        1.0      /   0.0
    fx              1.0/ -- / --                         cos(azi) /   sin(azi)
    fy              -- / 1.0/ --                         sin(azi) /  -cos(azi)
    ============================================================================
    """
    [M11, M12, M13, M22, M23, M33] = check_convert_fm(focal_mechanism=focal_mechanism)
    exp = (M11 + M22 + M33) / 3
    ss1 = M12
    ss2 = (M11 - M22) / 2
    ds1 = M13
    ds2 = M23
    clvd = (-M11 / 2 - M22 / 2 + M33) / (3 / 2)

    az = np.deg2rad(az_in_deg)
    sin_az, cos_az = np.sin(az), np.cos(az)
    sin_2az, cos_2az = np.sin(2 * az), np.cos(2 * az)
    m1 = [exp, ss1 * sin_2az + ss2 * cos_2az, ds1 * cos_az + ds2 * sin_az, clvd]
    m2 = [ss1 * cos_2az - ss2 * sin_2az, ds1 * sin_az - ds2 * cos_az]
    z = (
        time_series[0] * m1[0]
        + time_series[2] * m1[1]
        + time_series[5] * m1[2]
        + time_series[8] * m1[3]
    )
    r = (
        time_series[1] * m1[0]
        + time_series[3] * m1[1]
        + time_series[6] * m1[2]
        + time_series[9] * m1[3]
    )
    t = time_series[4] * m2[0] + time_series[7] * m2[1]

    seismograms = np.array([r, -t, -z])
    return seismograms


def seek_qseis06(
    path_green,
    event_depth_km,
    receiver_depth_km,
    az_deg,
    dist_km,
    focal_mechanism,
    srate,
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
    sampling_num = green_info['sampling_num']
    time_reduction_velo = green_info["time_reduction_velo"]
    dist_range = green_info["dist_range"]
    delta_dist = green_info["delta_dist"]
    num_each_group = green_info["N_each_group"]
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]
    if not isinstance(grn_dep_list, list):
        grn_dep = grn_dep_list
    else:
        grn_dep = grn_dep_list[
            np.argmin(np.abs(event_depth_km - np.array(grn_dep_list)))]
    if not isinstance(grn_receiver_list, list):
        grn_receiver = grn_receiver_list
    else:
        grn_receiver = grn_receiver_list[
            np.argmin(np.abs(receiver_depth_km - np.array(grn_receiver_list)))]

    path_greenfunc = str(
        os.path.join(path_green, "%.2f" % grn_dep, "%.2f" % grn_receiver)
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
    path_bin = os.path.join(path_greenfunc, "grn.bin")

    if os.path.exists(path_bin):
        time_series = read_time_series_qseis06_bin(
            path_bin=path_bin,
            start_count=start_count,
            sampling_num=sampling_num,
        )
    else:
        time_series = read_time_series_qseis06_ascii(
            path_greenfunc=path_greenfunc_sub, start_count=start_count
        )

    # r,t,z
    seismograms = synthesize_qseis06(
        az_in_deg=az_deg, time_series=time_series, focal_mechanism=focal_mechanism
    )
    if rotate:
        seismograms = rotate_rtz_to_enz(
            az_in_deg=az_deg, r=seismograms[0], t=seismograms[1], z=seismograms[2]
        )

    tpts_table = None
    if (before_p is not None) or shift or pad_zeros:
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

    conv_shift = round(green_info['wavelet_duration'] / 2)
    if conv_shift != 0:
        seismograms = np.roll(seismograms, -conv_shift)
        seismograms[:, -conv_shift:] = 0

    seismograms_resample = np.zeros((3, round(sampling_num * srate / srate_grn)))
    for i in range(3):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
        )

    if only_seismograms:
        return seismograms_resample
    else:
        return seismograms_resample, tpts_table, first_p, first_s, grn_dep, green_dist


if __name__ == "__main__":
    pass
