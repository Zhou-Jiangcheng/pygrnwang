import os
import json

import numpy as np


from .read_qseis import (
    read_time_series_qseis06_bin,
    read_time_series_qseis06_ascii,
    synthesize_qseis06,
)
from .utils import (
    shift_green2real_tpts,
    read_nd,
)
from .create_qseis_bulk import create_order_ind
from .focal_mechanism import check_convert_fm
from .geo import rotate_symmetric_tensor_series
from .pytaup import cal_first_p_s
from .signal_process import resample


def rotate_strain_polar2cartesian(theta, e_rr, e_rt, e_tt):
    theta = np.deg2rad(theta)
    cos_theta = np.cos(theta)
    sin_theta = np.sin(theta)
    cos_theta_sq = cos_theta**2
    sin_theta_sq = sin_theta**2
    sin_time_cos_theta = sin_theta * cos_theta
    e_xx = e_rr * cos_theta_sq + e_tt * sin_theta_sq - 2 * e_rt * sin_time_cos_theta
    e_yy = e_rr * sin_theta_sq + e_tt * cos_theta_sq + 2 * e_rt * sin_time_cos_theta
    e_xy = e_rt * (cos_theta_sq - sin_theta_sq) + (e_tt - e_rr) * sin_time_cos_theta
    return e_xx, e_xy, e_yy


def diff_central_1order(v_array, diff_accu_order):
    if diff_accu_order == 2:
        return (v_array[2] - v_array[0]) / 2
    elif diff_accu_order == 4:
        return (
            1 / 12 * v_array[0]
            - 2 / 3 * v_array[1]
            + 2 / 3 * v_array[3]
            - 1 / 12 * v_array[4]
        )
    elif diff_accu_order == 6:
        return (
            -1 / 60 * v_array[0]
            + 3 / 20 * v_array[1]
            - 3 / 4 * v_array[2]
            + 3 / 4 * v_array[4]
            - 3 / 20 * v_array[5]
            + 1 / 60 * v_array[6]
        )
    elif diff_accu_order == 8:
        return (
            1 / 280 * v_array[0]
            - 4 / 105 * v_array[1]
            + 1 / 5 * v_array[2]
            - 4 / 5 * v_array[3]
            + 4 / 5 * v_array[5]
            - 1 / 5 * v_array[6]
            + 4 / 105 * v_array[7]
            - 1 / 280 * v_array[8]
        )
    else:
        raise ValueError("diff_accu_order must be in [2,4,6,8]")


def seek_qseis06_strain_rate_diff(
    path_green,
    event_depth_km,
    receiver_depth_km,
    az_deg,
    dist_km,
    focal_mechanism,
    srate,
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
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]
    diff_accu_order = green_info["diff_accu_order"]
    k_dr = green_info["k_dr"]
    dz = green_info["dz"] * 1e3
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
    dr = grn_dist * k_dr * 1e3
    start_count = ind - ind_group * num_each_group

    if time_reduction_velo != 0:
        time_reduction = grn_dist / time_reduction_velo
    else:
        time_reduction = 0

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
    pm1_paz = [
        exp,
        ss1 * 2 * cos_2az + ss2 * (-2) * sin_2az,
        ds1 * (-1) * sin_az + ds2 * cos_az,
        clvd,
    ]
    pm2_paz = [
        ss1 * (-2) * sin_2az - ss2 * 2 * cos_2az,
        ds1 * cos_az - ds2 * (-1) * sin_az,
    ]

    # partial u_i / partial r
    velo_raw_dr = np.zeros((3 * (diff_accu_order + 1), sampling_num))
    for order in range(diff_accu_order + 1):
        order_ind = create_order_ind(order, diff_accu_order)
        path_greenfunc_order = str(
            os.path.join(path_greenfunc, "%d_%d" % (ind_group, order_ind))
        )
        path_bin = os.path.join(path_greenfunc_order, "grn_tz.npy")
        if os.path.exists(path_bin):
            time_series_list = read_time_series_qseis06_bin(
                path_greenfunc=path_greenfunc_order, start_count=start_count
            )
        else:
            time_series_list = read_time_series_qseis06_ascii(
                path_greenfunc=path_greenfunc_order, start_count=start_count
            )
        velo_dr = synthesize_qseis06(time_series_list=time_series_list, m1=m1, m2=m2)
        velo_raw_dr[order * 3 : (order + 1) * 3, :] = velo_dr
    pur_pr = diff_central_1order(velo_raw_dr[0:-1:3, :], diff_accu_order) / dr
    put_pr = diff_central_1order(velo_raw_dr[1:-1:3, :], diff_accu_order) / dr
    puz_pr = (
        diff_central_1order(velo_raw_dr[2 : len(velo_raw_dr) : 3, :], diff_accu_order)
        / dr
    )
    np.save("/home/zjc/Desktop/put_pr_diff.npy", put_pr)
    np.save("/home/zjc/Desktop/pur_pr_diff.npy", pur_pr)

    # partial u_i / partial theta
    path_greenfunc_order = str(os.path.join(path_greenfunc, "%d_%d" % (ind_group, 0)))
    path_bin = os.path.join(path_greenfunc_order, "grn_tz.npy")
    if os.path.exists(path_bin):
        time_series_list = read_time_series_qseis06_bin(
            path_greenfunc=path_greenfunc_order, start_count=start_count
        )
    else:
        time_series_list = read_time_series_qseis06_ascii(
            path_greenfunc=path_greenfunc_order, start_count=start_count
        )
    velo_dt = synthesize_qseis06(
        time_series_list=time_series_list, m1=pm1_paz, m2=pm2_paz
    )
    pur_pt = velo_dt[0]
    put_pt = velo_dt[1]
    puz_pt = velo_dt[2]

    # calculate err, ert, ett
    ur0 = velo_raw_dr[diff_accu_order // 2, :]
    ut0 = velo_raw_dr[diff_accu_order // 2 + 1, :]
    r = grn_dist * 1e3  # m
    e_rr = pur_pr
    e_tt = (put_pt + ur0) / r
    e_rt = 1 / 2 * (pur_pt / r + put_pr - ut0 / r)

    # partial u_i / partial z
    # calculate erz, etz, ezz
    if receiver_depth_km > 0:
        velo_raw_dz = np.zeros((3 * (diff_accu_order + 1), sampling_num))
        for order in range(diff_accu_order + 1):
            if order == diff_accu_order // 2:
                continue
            order_ind = create_order_ind(order, diff_accu_order) + diff_accu_order
            path_greenfunc_order = str(
                os.path.join(path_greenfunc, "%d_%d" % (ind_group, order_ind))
            )
            path_bin = os.path.join(path_greenfunc_order, "grn_tz.npy")
            if os.path.exists(path_bin):
                time_series_list = read_time_series_qseis06_bin(
                    path_greenfunc=path_greenfunc_order, start_count=start_count
                )
            else:
                time_series_list = read_time_series_qseis06_ascii(
                    path_greenfunc=path_greenfunc_order, start_count=start_count
                )
            velo_dz = synthesize_qseis06(
                time_series_list=time_series_list, m1=m1, m2=m2
            )
            velo_raw_dz[order * 3 : (order + 1) * 3, :] = velo_dz
        pur_pz = diff_central_1order(velo_raw_dz[0:-1:3, :], diff_accu_order) / dz
        put_pz = diff_central_1order(velo_raw_dz[1:-1:3, :], diff_accu_order) / dz
        puz_pz = (
            diff_central_1order(
                velo_raw_dz[2 : len(velo_raw_dr) : 3, :], diff_accu_order
            )
            / dz
        )
        e_rz = 1 / 2 * (pur_pz + puz_pr)
        e_tz = 1 / 2 * (put_pz + puz_pt / r)
        e_zz = puz_pz
    else:
        path_nd_noQ = green_info["path_nd_without_Q"]
        nd_model = read_nd(path_nd_noQ)
        depth_array = nd_model[:, 0]
        ind = np.argmin(np.abs(depth_array - receiver_depth_km))
        vp = nd_model[ind, 1] * 1e3
        vs = nd_model[ind, 2] * 1e3
        rho = nd_model[ind, 3] * 1e3
        mu = rho * vs**2
        lam = rho * vp**2 - 2 * mu
        e_rz = np.zeros_like(e_rr)
        e_tz = np.zeros_like(e_rr)
        e_zz = -lam / (lam + 2 * mu) * (e_rr + e_tt)

    strain_tensor = np.array([e_tt, e_rt, -e_tz, e_rr, -e_rz, e_zz])
    seismograms = rotate_symmetric_tensor_series(strain_tensor.T, np.deg2rad(-az_deg)).T

    tpts_table = None
    if (before_p is not None) or shift:
        first_p_grn, first_s_grn = cal_first_p_s(
            event_depth_km=event_depth_km,
            receiver_depth_km=receiver_depth_km,
            dist_km=grn_dist,
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

    seismograms_resample = np.zeros((6, round(sampling_num * srate / srate_grn)))
    for i in range(6):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
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


def convert_strain2stress(strain, lam, mu):
    # ee en ez nn nz zz
    stress = np.zeros_like(strain)
    s = strain[0, :] + strain[3, :] + strain[5, :]
    stress[0, :] = lam * s + 2 * mu * strain[0, :]
    stress[1, :] = 2 * mu * strain[1, :]
    stress[2, :] = 2 * mu * strain[2, :]
    stress[3, :] = lam * s + 2 * mu * strain[3, :]
    stress[4, :] = 2 * mu * strain[4, :]
    stress[5, :] = lam * s + 2 * mu * strain[5, :]
    return stress


def seek_qseis06_stress_rate_diff(
    path_green,
    event_depth_km,
    receiver_depth_km,
    az_deg,
    dist_km,
    focal_mechanism,
    srate,
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
    sampling_num = (
        round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    )
    (
        strain_rate,
        tpts_table,
        first_p,
        first_s,
        grn_dep_source,
        grn_dep_receiver,
        grn_dist,
    ) = seek_qseis06_strain_rate_diff(
        path_green,
        event_depth_km,
        receiver_depth_km,
        az_deg,
        dist_km,
        focal_mechanism,
        srate_grn,
        before_p,
        pad_zeros,
        shift,
        False,
        model_name,
        green_info,
    )

    path_nd_noQ = green_info["path_nd_without_Q"]
    nd_model = read_nd(path_nd_noQ)
    depth_array = nd_model[:, 0]
    ind = np.argmin(np.abs(depth_array - receiver_depth_km))
    vp = nd_model[ind, 1] * 1e3
    vs = nd_model[ind, 2] * 1e3
    rho = nd_model[ind, 3] * 1e3
    mu = rho * vs**2
    lam = rho * vp**2 - 2 * mu
    seismograms = convert_strain2stress(strain_rate, lam, mu)

    seismograms_resample = np.zeros((6, round(sampling_num * srate / srate_grn)))
    for i in range(6):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
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
