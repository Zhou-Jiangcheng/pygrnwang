import os
import json

import numpy as np
import pandas as pd
import scipy.signal as signal

from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz, rotate_symmetric_tensor_series
from .pytaup import read_tpts_table
from .utils import read_layerd_material, shift_green2real_tpts
from .signal_process import resample

# Define output type lists
one_com_list = ["volume"]
three_com_list = ["disp", "velo"]
rota_com_list = ["rota", "rota_rate"]
six_com_list = ["strain", "strain_rate", "stress", "stress_rate"]


def get_outfile_name_list(output_type):
    if output_type == "volume":
        name_list_psv = ["tv"]
        name_list_sh = []
    elif output_type == "disp" or output_type == "velo":
        name_list_psv = ["tz", "tr"]
        name_list_sh = ["tt"]
    elif output_type == "strain" or output_type == "strain_rate":
        name_list_psv = ["ezz", "ezr", "err", "ett"]
        name_list_sh = ["ezt", "ert"]
    elif output_type == "stress" or output_type == "stress_rate":
        name_list_psv = ["szz", "szr", "srr", "stt"]
        name_list_sh = ["szt", "srt"]
    elif output_type == "rota" or output_type == "rota_rate":
        name_list_psv = ["ot"]
        name_list_sh = ["oz", "or"]
    else:
        raise ValueError(
            "output_type must in  disp | velo | strain | strain_rate | "
            "stress | stress_rate | rota | rota_rate"
        )
    return name_list_psv, name_list_sh


def read_time_series_qseis2025_ascii(
        path_greenfunc, start_count, output_type="disp"
):
    time_seris_list = []
    start_count = start_count + 1
    name_list_psv, name_list_sh = get_outfile_name_list(output_type)
    for com in name_list_psv:
        ex_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "ex.%s" % com)), sep="\\s+"
        ).to_numpy()  # type:ignore
        ss_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "ss.%s" % com)), sep="\\s+"
        ).to_numpy()  # type:ignore
        ds_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "ds.%s" % com)), sep="\\s+"
        ).to_numpy()  # type:ignore
        cl_com = pd.read_csv(
            str(os.path.join(path_greenfunc, "cl.%s" % com)), sep="\\s+"
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
    for com in name_list_sh:
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
    return name_list_psv, name_list_sh, time_seris_list


def read_time_series_qseis2025_bin(
        path_greenfunc,
        start_count,
        output_type,
):
    time_series_list = []
    name_list_psv, name_list_sh = get_outfile_name_list(output_type)
    for com in name_list_psv:
        time_series_com = np.load(os.path.join(path_greenfunc, "grn_%s.npy" % com))
        time_series_list.append(time_series_com[:, start_count].reshape(4, -1))
    for com in name_list_sh:
        time_series_com = np.load(os.path.join(path_greenfunc, "grn_%s.npy" % com))
        time_series_list.append(time_series_com[:, start_count].reshape(2, -1))
    return name_list_psv, name_list_sh, time_series_list


def synthesize_rzv(time_series, m1):
    # ex,ss,ds,cl
    rzv = (
            time_series[0] * m1[0]
            + time_series[1] * m1[1]
            + time_series[2] * m1[2]
            + time_series[3] * m1[3]
    )
    return rzv


def synthesize_t(time_series, m2):
    t = time_series[0] * m2[0] + time_series[1] * m2[1]
    return t


def seek_qseis2025(
        path_green,
        event_depth_km,
        receiver_depth_km,
        az_deg,
        dist_km,
        focal_mechanism,
        srate,
        output_type,
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
    grn_dist = dist_range[0] + ind * delta_dist
    start_count = ind - ind_group * num_each_group

    path_greenfunc_sub = os.path.join(path_greenfunc, "%d_0" % ind_group)
    if os.path.exists(os.path.join(path_greenfunc_sub, "grn_szt.npy")):
        name_list_psv, name_list_sh, time_series_list = (
            read_time_series_qseis2025_bin(
                path_greenfunc=path_greenfunc_sub,
                start_count=start_count,
                output_type=output_type,
            )
        )
    else:
        name_list_psv, name_list_sh, time_series_list = (
            read_time_series_qseis2025_ascii(
                path_greenfunc=path_greenfunc_sub,
                start_count=start_count,
                output_type=output_type,
            )
        )
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

    az_rad = np.deg2rad(az_deg)
    sin_az, cos_az = np.sin(az_rad), np.cos(az_rad)
    sin_2az, cos_2az = np.sin(2 * az_rad), np.cos(2 * az_rad)

    m1 = [exp, ss1 * sin_2az + ss2 * cos_2az, ds1 * cos_az + ds2 * sin_az, clvd]
    m2 = [ss1 * cos_2az - ss2 * sin_2az, ds1 * sin_az - ds2 * cos_az]

    # p-sv
    # ["tv,
    # "tr", "tz",
    # "ezz", "ezr", "err", "ett",
    # "szz", "szr", "srr", "stt",
    # "ot"]
    # sh
    # ["tt",
    # "ezt", "ert",
    # "szt", "srt",
    # "oz", "or"]
    if output_type in one_com_list:
        uv = synthesize_rzv(time_series=time_series_list[0], m1=m1)
        seismograms = np.array(uv)
    elif output_type in three_com_list:
        uz = -synthesize_rzv(time_series=time_series_list[0], m1=m1)
        ur = synthesize_rzv(time_series=time_series_list[1], m1=m1)
        ut = -synthesize_t(time_series=time_series_list[2], m2=m2)
        if rotate:
            seismograms = rotate_rtz_to_enz(az_deg, r=ur, t=ut, z=uz)
        else:
            seismograms = np.array([ur, ut, uz])
    elif output_type in rota_com_list:
        omega_t = -synthesize_rzv(time_series=time_series_list[0], m1=m1)
        omega_z = -synthesize_t(time_series=time_series_list[1], m2=m2)
        omega_r = synthesize_t(time_series=time_series_list[2], m2=m2)
        if rotate:
            seismograms = rotate_rtz_to_enz(az_deg, r=omega_r, t=omega_t, z=omega_z)
        else:
            seismograms = np.array([omega_r, omega_t, omega_z])
    elif output_type in six_com_list:
        s_zz = synthesize_rzv(time_series=time_series_list[0], m1=m1)
        s_zr = synthesize_rzv(time_series=time_series_list[1], m1=m1)
        s_rr = synthesize_rzv(time_series=time_series_list[2], m1=m1)
        s_tt = synthesize_rzv(time_series=time_series_list[3], m1=m1)
        s_zt = synthesize_t(time_series=time_series_list[4], m2=m2)
        s_rt = synthesize_t(time_series=time_series_list[5], m2=m2)

        sigma_tensor = np.array([s_tt, s_rt, -s_zt, s_rr, -s_zr, s_zz])
        if rotate:
            seismograms = rotate_symmetric_tensor_series(sigma_tensor.T, az_rad).T
        else:
            seismograms = sigma_tensor
    else:
        raise ValueError(
            "output_type must in  disp | velo | strain | strain_rate | "
            "stress | stress_rate | rota | rota_rate"
        )

    tpts_table = None
    if (before_p is not None) or shift or pad_zeros:
        first_p_grn, first_s_grn = read_tpts_table(
            path_green=path_green,
            event_depth_km=grn_dep_source,
            receiver_depth_km=grn_dep_receiver,
            ind=ind,
        )
        tpts_table = {"p_onset": first_p_grn, "s_onset": first_s_grn}

    if time_reduction_velo != 0:
        time_reduction = grn_dist / time_reduction_velo
    else:
        time_reduction = 0

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

    seismograms_resample = np.zeros(
        (len(seismograms), round(sampling_num * srate / srate_grn))
    )
    for i in range(len(seismograms)):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
        )
    if (wavelet_type == 1) and ("rate" not in output_type):
        seismograms_resample = np.cumsum(seismograms_resample, axis=1) / srate
    elif (wavelet_type == 2) and (("rate" in output_type) or (output_type == "velo")):
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
            grn_dist,
        )


if __name__ == "__main__":
    pass
