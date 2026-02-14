import os
import json
from typing import Union

import numpy as np
import pandas as pd

from .create_qssp2020 import mt_com_list
from .utils import (
    shift_green2real_tpts,
    create_rotate_z_mat,
    rotate_symmetric_tensor_series, read_tpts_table,
)
from .focal_mechanism import (
    convert_mt_axis,
    tensor2full_tensor_matrix,
    check_convert_fm,
)
from .geo import rotate_rtz_to_enz
from .signal_process import resample, filter_butter

# Define output type lists
one_com_list = ["gravimeter"]
three_com_list = ["disp", "velo", "acce", "rota", "rota_rate", "gravitation"]
six_com_list = ["strain", "strain_rate", "stress", "stress_rate"]


def read_time_series_qssp2020(path_bin, ind, sampling_num):
    fr = open(path_bin, "rb")
    time_series = np.fromfile(
        file=fr, dtype=np.float32, count=sampling_num, offset=ind * sampling_num * 4
    )
    fr.close()
    return time_series


def seek_raw_qssp2020(
    path_green,
    event_depth,
    receiver_depth,
    dist,
    output_type="disp",
    green_info=None,
):
    """
    read raw data at stations located on the straight line from the origin to the north
    moment tensor in rtp axis
    data in enz axis
    """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    nt_cut = round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    if output_type == "gravimeter":
        enz_list = [""]
        data_all = np.zeros((nt_cut, 6))
    elif output_type in ["disp", "velo", "acce", "rota", "rota_rate", "gravitation"]:
        enz_list = ["_e", "_n", "_z"]
        data_all = np.zeros((nt_cut, 18))
    elif output_type in ["stress", "stress_rate", "strain", "strain_rate"]:
        enz_list = ["_ee", "_en", "_ez", "_nn", "_nz", "_zz"]
        data_all = np.zeros((nt_cut, 36))
    else:
        raise ValueError(
            "output_type must in  disp | velo | acce | strain | strain_rate | "
            "stress | stress_rate | rotation | rotation_rate | gravitation | gravimeter"
        )

    for i_rtp in range(6):
        path_func = str(
            os.path.join(
                path_green,
                "GreenFunc",
                "%.2f" % event_depth,
                "%.2f" % receiver_depth,
                mt_com_list[i_rtp],
                "",
            )
        )
        for i_enz in range(len(enz_list)):
            dist_range = green_info["grn_dist_range"]
            delta_dist = green_info["grn_delta_dist"]
            ind = round((dist - dist_range[0]) / delta_dist)
            path_bin = os.path.join(
                path_func, "_%s%s.bin" % (output_type, enz_list[i_enz])
            )
            # path_bin = ''
            if os.path.exists(path_bin):
                data_enz = read_time_series_qssp2020(
                    path_bin=path_bin, ind=ind, sampling_num=nt_cut
                )
                data_all[:, i_enz + len(enz_list) * i_rtp] = data_enz
            else:
                data_enz = pd.read_csv(
                    os.path.join(
                        path_func, "_%s%s.dat" % (output_type, enz_list[i_enz])
                    ),
                    sep="\\s+",
                ).to_numpy()[:, 1:]
                data_all[:, i_enz + len(enz_list) * i_rtp] = data_enz[:, ind]
    return data_all


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


def seek_qssp2020(
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
    freq_band=None,
    butter_order: int = 4,
    zero_phase: bool = False,
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
    :param output_type:
        one_com - "gravimeter"
        three_com - "disp", "velo", "acce", "rota", "rota_rate", "gravitation"
        six_com - "strain", "strain_rate", "stress", "stress_rate"
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
    :param freq_band: Frequency band for bandpass filter [low_freq, high_freq] in Hz.
            Use None or [None, None] for no filtering (default).
            Use [low_freq, None] for highpass, [None, high_freq] for lowpass.
    :param butter_order: Order of Butterworth filter (default: 4).
    :param zero_phase: Whether to use zero-phase filtering (default: False).
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
    sampling_num = (
        round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    )
    dist_range = green_info["grn_dist_range"]
    delta_dist = green_info["grn_delta_dist"]
    time_reduction = green_info["time_reduction"]
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

    # Nearest distance logic for metadata
    dist_min = dist_range[0]
    float_ind = (dist_km - dist_min) / delta_dist
    nearest_indice = round(float_ind)
    grn_dist = dist_range[0] + nearest_indice * delta_dist

    # --- 2. Retrieve Raw Green's Function Data based on Interpolation Type ---

    if interpolate_type == 0:
        # === Type 0: Nearest Neighbor (Fastest) ===
        raw_final = seek_raw_qssp2020(
            path_green,
            grn_dep_source,
            grn_dep_receiver,
            grn_dist,
            output_type,
            green_info,
        )

    else:
        # === Type 1: Trilinear Interpolation (Src Depth, Rec Depth, Distance) ===

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

        dist_low_km = dist_min + ind_low * delta_dist
        dist_high_km = dist_min + (ind_low + 1) * delta_dist

        # D. Helper to fetch raw
        def fetch_raw(d_src, d_rec, d_dist):
            return seek_raw_qssp2020(
                path_green, d_src, d_rec, d_dist, output_type, green_info
            )

        # E. Nested Interpolation Logic
        # E1. Interpolate Distance (Innermost)
        def get_dist_interp_data(src_depth, rec_depth):
            raw_d0 = fetch_raw(src_depth, rec_depth, dist_low_km)
            if w_dist > 0:
                raw_d1 = fetch_raw(src_depth, rec_depth, dist_high_km)
                return (1 - w_dist) * raw_d0 + w_dist * raw_d1
            else:
                return raw_d0

        # E2. Interpolate Receiver Depth (Middle)
        def get_rec_interp_data(src_depth):
            data_r0 = get_dist_interp_data(src_depth, d_rec_low)
            if w_rec > 1e-4 and d_rec_high != d_rec_low:
                data_r1 = get_dist_interp_data(src_depth, d_rec_high)
                return (1 - w_rec) * data_r0 + w_rec * data_r1
            else:
                return data_r0

        # E3. Interpolate Source Depth (Outermost)
        data_low_src = get_rec_interp_data(d_src_low)
        if w_src > 1e-4 and d_src_high != d_src_low:
            data_high_src = get_rec_interp_data(d_src_high)
            raw_final = (1 - w_src) * data_low_src + w_src * data_high_src
        else:
            raw_final = data_low_src

    # --- 3. Prepare MT and Rotation (Common) ---
    gamma = np.deg2rad(az_deg)
    A_rotate = create_rotate_z_mat(gamma=gamma)
    focal_mechanism = check_convert_fm(focal_mechanism)
    mt_ned_full = tensor2full_tensor_matrix(mt=focal_mechanism, flag="ned")
    mt_rotate = A_rotate.T @ mt_ned_full @ A_rotate
    mt_list = np.array(
        [
            mt_rotate[0, 0],
            mt_rotate[0, 1],
            mt_rotate[0, 2],
            mt_rotate[1, 1],
            mt_rotate[1, 2],
            mt_rotate[2, 2],
        ]
    )
    mt_rtp = convert_mt_axis(mt_list, "ned2rtp")

    # --- 4. Process Components (Using the interpolated raw_final data) ---
    if output_type in one_com_list:
        u = raw_final
        u_enz_green_north = np.zeros((sampling_num, 1))
        for i_rtp in range(6):
            u_enz_green_north[:, 0] += mt_rtp[i_rtp] * u[:, i_rtp]
        seismograms = u_enz_green_north.T

    elif output_type in three_com_list:
        u_all = raw_final
        u_enz_green_north = np.zeros((sampling_num, 3))
        for i_rtp in range(6):
            for i_enz in range(3):
                u_enz_green_north[:, i_enz] += (
                    mt_rtp[i_rtp] * u_all[:, i_enz + 3 * i_rtp]
                )
        u_rtz = np.zeros((sampling_num, 3))
        u_rtz[:, 0] = u_enz_green_north[:, 1]
        u_rtz[:, 1] = -u_enz_green_north[:, 0]
        u_rtz[:, 2] = u_enz_green_north[:, 2]
        seismograms = u_rtz.T

        if rotate:
            seismograms = rotate_rtz_to_enz(
                az_deg, r=u_rtz[:, 0], t=u_rtz[:, 1], z=u_rtz[:, 2]
            )

    elif output_type in six_com_list:
        epsilon_all = raw_final
        epsilon_enz_green_north = np.tensordot(
            epsilon_all.reshape(sampling_num, 6, 6), mt_rtp, axes=(1, 0)
        )
        if rotate:
            seismograms = rotate_symmetric_tensor_series(
                epsilon_enz_green_north, gamma
            ).T
        else:
            seismograms = epsilon_enz_green_north.T
    else:
        raise ValueError(
            "output_type must be one of: disp, velo, acce, rota, rota_rate, gravimeter, gravitation, strain, strain_rate, stress, stress_rate"
        )

    # --- 5. Post-processing (Time shifting, resampling) ---
    # TPTS tables are read from the nearest neighbor grid point
    tpts_table = None
    if (before_p is not None) or shift or pad_zeros:
        grn_first_p, grn_first_s = read_tpts_table(
            path_green=os.path.join(path_green, "GreenFunc"),
            event_depth_km=grn_dep_source,
            receiver_depth_km=grn_dep_receiver,
            ind=nearest_indice,
        )
        tpts_table = {"p_onset": grn_first_p, "s_onset": grn_first_s}

    # Apply bandpass filter
    if freq_band is not None and (freq_band[0] is not None or freq_band[1] is not None):
        for i in range(len(seismograms)):
            seismograms[i] = filter_butter(
                seismograms[i], srate_grn, freq_band, butter_order, zero_phase
            )
    
    ts_count = 0
    if before_p is not None:
        ts_count = round(
            (tpts_table["p_onset"] - time_reduction - before_p) * srate_grn
        )
    if pad_zeros and time_reduction != 0:
        if before_p is not None:
            raise ValueError("can not set before_p and pad_zeros together")
        ts_count = round(-time_reduction * srate_grn)
        before_p = tpts_table["p_onset"]
    seismograms = np.roll(seismograms, -ts_count)
    if ts_count > 0:
        seismograms[:, -ts_count:] = 0
    elif ts_count < 0:
        seismograms[:, :-ts_count] = 0

    if shift:
        seismograms, first_p, first_s = shift_green2real_tpts(
            seismograms=seismograms,
            tpts_table=tpts_table,
            srate=srate_grn,
            green_before_p=before_p,
            event_depth_km=event_depth_km,
            dist_in_km=dist_km,
            receiver_depth_km=receiver_depth_km,
            model_name=model_name,
        )
    else:
        first_p = None
        first_s = None

    seismograms_resample = np.zeros(
        (seismograms.shape[0], round(sampling_num * srate / srate_grn))
    )
    for i in range(seismograms.shape[0]):
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
