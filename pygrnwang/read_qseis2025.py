import os
import json
from typing import Union

import numpy as np
import pandas as pd
import scipy.signal as signal

from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz, rotate_symmetric_tensor_series
from .utils import shift_green2real_tpts, read_tpts_table
from .signal_process import resample

# Define output type lists
one_com_list = ["volume"]
three_com_list = ["disp", "velo", "acce"]
rota_com_list = ["rota", "rota_rate"]
six_com_list = ["strain", "strain_rate", "stress", "stress_rate"]


def get_outfile_name_list(output_type):
    if output_type == "volume":
        name_list_psv = ["tv"]
        name_list_sh = []
    elif output_type == "disp" or output_type == "velo" or output_type == "acce":
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


def read_time_series_qseis2025_ascii(path_greenfunc, start_count, output_type="disp"):
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


def get_sorted_grid_params(target, grid_list):
    """
    Helper to find neighbors and weight for 1D interpolation on a sorted list.
    Returns: val_low, val_high, weight_high
    """
    arr = np.array(grid_list)
    if len(arr) == 0:
        raise ValueError("Grid list is empty")
    if len(arr) == 1:
        return arr[0], arr[0], 0.0

    # Handle out of bounds by clamping
    if target <= arr[0]:
        return arr[0], arr[0], 0.0
    if target >= arr[-1]:
        return arr[-1], arr[-1], 0.0

    idx = np.searchsorted(arr, target)
    # arr[idx-1] <= target <= arr[idx]
    val_low = arr[idx - 1]
    val_high = arr[idx]
    weight = (target - val_low) / (val_high - val_low)
    return val_low, val_high, weight


def seek_qseis2025(
        path_green: str,
        event_depth_km: float,
        receiver_depth_km: float,
        az_deg: float,
        dist_km: float,
        focal_mechanism: Union[np.ndarray, list],
        srate: float,
        output_type: str = "disp",
        rotate: bool = True,
        before_p: Union[float, None] = None,
        pad_zeros: bool = False,
        shift: bool = False,
        only_seismograms: bool = True,
        model_name: str = "ak135fc",
        green_info: Union[dict, None] = None,
        interpolate_type: int = 0,
):
    """
    Read synthetic seismograms.

    :param path_green: Root directory of the data.
    :param event_depth_km: Event depth in km.
    :param receiver_depth_km: Receiver depth in km.
    :param az_deg: Azimuth in degrees.
    :param dist_km: Epicentral distance in km.
    :param focal_mechanism: [strike, dip, rake] or [M11, M12, M13, M22, M23, M33].
    :param srate: Sampling rate in Hz.
    :param output_type: disp | velo | acce | strain | strain_rate |
            stress | stress_rate | rota | rota_rate.
    :param before_p: Time before P-wave.
    :param pad_zeros: Pad with zeros.
    :param shift: Shift seismograms based on tpts.
    :param rotate: Rotate rtz2ned.
    :param only_seismograms: Return only seismograms.
    :param model_name: Model name.
    :param green_info: Green's function library info.
    :param interpolate_type:
            0 for nearest neighbor,
            1 for trilinear interpolation (Source Depth, Receiver Depth, Distance).
    :return: (
            seismograms_resample,
            tpts_table,
            first_p,
            first_s,
            grn_dep_source,
            grn_dep_receiver,
            grn_dist,
        )
    """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    srate_grn = 1 / green_info["sampling_interval"]
    sampling_num = green_info["sampling_num"]
    time_reduction_velo = green_info["time_reduction_velo"]
    dist_range = green_info["grn_dist_range"]
    delta_dist = green_info["grn_delta_dist"]
    num_each_group = green_info["N_each_group"]
    wavelet_type = green_info["wavelet_type"]
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]

    # --- 1. Identify Nearest Neighbors (Used for Metadata) ---
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

    # Calculate global distance index (float for interp, int for nearest)
    dist_min = dist_range[0]
    float_ind = (dist_km - dist_min) / delta_dist
    nearest_indice = max(0, round(float_ind))
    grn_dist = dist_range[0] + nearest_indice * delta_dist

    # --- 2. Helper to fetch raw data for a specific (source, receiver, distance) ---
    def fetch_raw_green_data(src_depth, rec_depth, dist_idx):
        path_greenfunc = str(
            os.path.join(path_green, "%.2f" % src_depth, "%.2f" % rec_depth)
        )
        # Ensure index is within bounds
        dist_idx = int(max(0, dist_idx))

        ind_group = dist_idx // num_each_group
        start_count = dist_idx - ind_group * num_each_group

        path_greenfunc_sub = os.path.join(path_greenfunc, "%d_0" % ind_group)
        if os.path.exists(os.path.join(path_greenfunc_sub, "grn_szt.npy")):
            _, _, ts_list = read_time_series_qseis2025_bin(
                path_greenfunc=path_greenfunc_sub,
                start_count=start_count,
                output_type=output_type,
            )
        else:
            _, _, ts_list = read_time_series_qseis2025_ascii(
                path_greenfunc=path_greenfunc_sub,
                start_count=start_count,
                output_type=output_type,
            )
        return ts_list

    # --- 3. Retrieve Data based on Interpolation Type ---
    if interpolate_type == 0:
        # === Type 0: Nearest Neighbor ===
        time_series_list = fetch_raw_green_data(
            grn_dep_source, grn_dep_receiver, nearest_indice
        )

    else:
        # === Type 1: Trilinear Interpolation (Source, Receiver, Distance) ===

        # A. Source Depth Interpolation Parameters
        if not isinstance(grn_dep_list, list):
            d_src_low, d_src_high, w_src = grn_dep_list, grn_dep_list, 0.0
        else:
            d_src_low, d_src_high, w_src = get_sorted_grid_params(
                event_depth_km, grn_dep_list
            )

        # B. Receiver Depth Interpolation Parameters (NEW)
        if not isinstance(grn_receiver_list, list):
            d_rec_low, d_rec_high, w_rec = grn_receiver_list, grn_receiver_list, 0.0
        else:
            d_rec_low, d_rec_high, w_rec = get_sorted_grid_params(
                receiver_depth_km, grn_receiver_list
            )

        # C. Distance Interpolation Parameters
        ind_low = int(np.floor(float_ind))
        if ind_low < 0:
            ind_low = 0
        w_dist = float_ind - ind_low
        if w_dist < 1e-4:
            w_dist = 0.0

        # D. Nested Interpolation Helpers

        # D1. Interpolate Distance (Innermost)
        def get_dist_interp_data(src_depth, rec_depth):
            raw_d0 = fetch_raw_green_data(src_depth, rec_depth, ind_low)
            if w_dist > 0:
                raw_d1 = fetch_raw_green_data(src_depth, rec_depth, ind_low + 1)
                return [
                    (1 - w_dist) * arr0 + w_dist * arr1
                    for arr0, arr1 in zip(raw_d0, raw_d1)
                ]
            else:
                return raw_d0

        # D2. Interpolate Receiver Depth (Middle)
        def get_rec_interp_data(src_depth):
            data_r0 = get_dist_interp_data(src_depth, d_rec_low)
            if w_rec > 1e-4 and d_rec_high != d_rec_low:
                data_r1 = get_dist_interp_data(src_depth, d_rec_high)
                return [
                    (1 - w_rec) * arr0 + w_rec * arr1
                    for arr0, arr1 in zip(data_r0, data_r1)
                ]
            else:
                return data_r0

        # E. Interpolate Source Depth (Outermost)
        data_low_src = get_rec_interp_data(d_src_low)
        if w_src > 1e-4 and d_src_high != d_src_low:
            data_high_src = get_rec_interp_data(d_src_high)
            time_series_list = [
                (1 - w_src) * arr0 + w_src * arr1
                for arr0, arr1 in zip(data_low_src, data_high_src)
            ]
        else:
            time_series_list = data_low_src

    # --- 4. Synthesize Seismograms (Common Logic) ---
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

    # --- 5. Post-Processing (TPTS, Time Shift, Resample) ---
    tpts_table = None
    if (before_p is not None) or shift or pad_zeros:
        first_p_grn, first_s_grn = read_tpts_table(
            path_green=path_green,
            event_depth_km=grn_dep_source,
            receiver_depth_km=grn_dep_receiver,
            ind=nearest_indice,
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
    elif (wavelet_type == 1) and (output_type == 'acce'):
        seismograms_resample = (
                signal.convolve(
                    seismograms_resample.T,
                    np.array([1, -1])[:, None],
                    mode="same",
                    method="auto",
                ).T
                / srate
        )
    elif (wavelet_type == 2) and (output_type == 'acce'):
        seismograms_resample = (
                signal.convolve(
                    seismograms_resample.T,
                    np.array([1, -2, 1])[:, None],
                    mode="same",
                    method="auto",
                ).T
                / (srate * srate)
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
