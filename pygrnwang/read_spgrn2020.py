import os
import json
from typing import Union

import numpy as np
from scipy import signal

from .read_tpts_table import read_tpts_table
from .utils import shift_green2real_tpts
from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz
from .signal_process import resample, filter_butter


def synthesize_spgrn(az_in_deg, time_series, focal_mechanism):
    [M11, M12, M13, M22, M23, M33] = check_convert_fm(focal_mechanism=focal_mechanism)
    # expl.,strike-slip,dip-slip,clvd
    exp = (M11 + M22 + M33) / 3
    ss1 = -M12
    ss2 = (M11 - M22) / 2
    ds1 = M13
    ds2 = -M23
    clvd = M33 - exp

    az = np.deg2rad(180 - az_in_deg)
    # az = np.deg2rad(az_in_deg)
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

    seismograms = np.array([r, t, z])
    return seismograms


def read_spgrn_data_by_index(path_grn_data, dist_index, green_info):
    """Helper to read raw data at a specific distance index."""
    N_T = round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    # Structure based on original code: 3 header floats + (2 + N_T) floats * 10 components
    length_each = 3 + (2 + N_T) * 10
    start_count = dist_index * length_each

    with open(path_grn_data, "rb") as fr:
        time_series = np.fromfile(
            file=fr, dtype=np.float32, count=length_each, offset=start_count * 4
        )

    # Reshape logic from original code
    time_series = time_series[3:].reshape(10, 2 + N_T)
    time_series = time_series[:, 1:-1]
    return time_series


def get_sorted_grid_params(target, grid_list):
    """
    Helper to find neighbors and weight for 1D interpolation on a sorted list.
    Returns: val_low, val_high, weight_high, idx_low
    """
    arr = np.array(grid_list)
    if len(arr) == 0:
        raise ValueError("Grid list is empty")
    if len(arr) == 1:
        return arr[0], arr[0], 0.0, 0

    # Handle out of bounds by clamping
    if target <= arr[0]:
        return arr[0], arr[0], 0.0, 0
    if target >= arr[-1]:
        return arr[-1], arr[-1], 0.0, len(arr) - 1

    idx = np.searchsorted(arr, target)
    # arr[idx-1] <= target <= arr[idx]
    val_low = arr[idx - 1]
    val_high = arr[idx]

    if val_high == val_low:
        weight = 0.0
    else:
        weight = (target - val_low) / (val_high - val_low)

    return val_low, val_high, weight, idx - 1


def seek_spgrn2020(
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
    :param output_type: 'disp', 'velo', 'acce'.
    :param srate: Sampling rate in Hz.
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
    sampling_num = green_info["samples_num"]
    green_before_p = green_info["green_before_p"]
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]
    dist_list = green_info["dist_list"]

    # --- 1. Identify Nearest Neighbors (Used for Metadata) ---
    # These are still needed for reading tpts tables which typically don't support interpolation well
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
    nearest_dist_idx = np.argmin(np.abs(np.array(dist_list) - dist_km))
    grn_dist = dist_list[nearest_dist_idx]

    # --- 2. Retrieve Data based on Interpolation Type ---

    if interpolate_type == 0:
        # === Type 0: Nearest Neighbor ===
        path_greenfunc = str(
            os.path.join(
                path_green,
                "GreenFunc",
                "%.2f" % grn_dep_source,
                "%.2f" % grn_dep_receiver,
            )
        )
        path_grn_data = os.path.join(path_greenfunc, "grn_d%.2f" % grn_dep_source)
        time_series = read_spgrn_data_by_index(
            path_grn_data, nearest_dist_idx, green_info
        )

    else:
        # === Type 1: Trilinear Interpolation (Source Depth & Receiver Depth & Distance) ===

        # A. Source Depth Interpolation Parameters
        if not isinstance(grn_dep_list, list):
            d_src_low, d_src_high, w_src, _ = grn_dep_list, grn_dep_list, 0.0, 0
        else:
            d_src_low, d_src_high, w_src, _ = get_sorted_grid_params(
                event_depth_km, grn_dep_list
            )

        # B. Receiver Depth Interpolation Parameters (NEW)
        if not isinstance(grn_receiver_list, list):
            d_rec_low, d_rec_high, w_rec, _ = (
                grn_receiver_list,
                grn_receiver_list,
                0.0,
                0,
            )
        else:
            d_rec_low, d_rec_high, w_rec, _ = get_sorted_grid_params(
                receiver_depth_km, grn_receiver_list
            )

        # C. Distance Interpolation Parameters
        _, _, w_dist, dist_idx_low = get_sorted_grid_params(dist_km, dist_list)
        dist_idx_high = min(dist_idx_low + 1, len(dist_list) - 1)

        if dist_idx_low == dist_idx_high:
            w_dist = 0.0

        # D. Helper to fetch and interpolate distance for a specific (Source, Receiver) pair
        def fetch_distance_layer(src_depth, rec_depth):
            path_gf_layer = str(
                os.path.join(
                    path_green,
                    "GreenFunc",
                    "%.2f" % src_depth,
                    "%.2f" % rec_depth,
                )
            )
            path_data_layer = os.path.join(path_gf_layer, "grn_d%.2f" % src_depth)

            # Read low distance
            data_d0 = read_spgrn_data_by_index(
                path_data_layer, dist_idx_low, green_info
            )

            if w_dist > 1e-4:
                # Read high distance
                data_d1 = read_spgrn_data_by_index(
                    path_data_layer, dist_idx_high, green_info
                )
                return (1 - w_dist) * data_d0 + w_dist * data_d1
            else:
                return data_d0

        # E. Helper to fetch and interpolate Receiver Depth for a specific Source Depth
        def fetch_receiver_layer(src_depth):
            # 1. Fetch Low Receiver Depth
            data_r0 = fetch_distance_layer(src_depth, d_rec_low)

            # 2. Fetch High Receiver Depth (if needed)
            if w_rec > 1e-4 and d_rec_high != d_rec_low:
                data_r1 = fetch_distance_layer(src_depth, d_rec_high)
                return (1 - w_rec) * data_r0 + w_rec * data_r1
            else:
                return data_r0

        # F. Perform Source Depth Interpolation (Top Level)
        # 1. Low Source Depth
        ts_low_src = fetch_receiver_layer(d_src_low)

        # 2. High Source Depth (if needed)
        if w_src > 1e-4 and d_src_high != d_src_low:
            ts_high_src = fetch_receiver_layer(d_src_high)
            time_series = (1 - w_src) * ts_low_src + w_src * ts_high_src
        else:
            time_series = ts_low_src

    # --- 3. Synthesize Seismograms ---
    # r,t,z
    seismograms = synthesize_spgrn(
        az_in_deg=az_deg, time_series=time_series, focal_mechanism=focal_mechanism
    )
    if rotate:
        seismograms = rotate_rtz_to_enz(
            az_in_deg=az_deg, r=seismograms[0], t=seismograms[1], z=seismograms[2]
        )[:]

    # --- 4. Post-processing (Time shifting, resampling) ---
    # Get metadata from the nearest neighbor path (Metadata uses Nearest Neighbor)
    path_greenfunc_meta = str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % grn_dep_source,
            "%.2f" % grn_dep_receiver,
        )
    )
    tpts_table = read_tpts_table(
        path_greenfunc=path_greenfunc_meta, dist_in_km=dist_km, green_info=green_info
    )

    # Apply bandpass filter
    if freq_band is not None and (freq_band[0] is not None or freq_band[1] is not None):
        for i in range(len(seismograms)):
            seismograms[i] = filter_butter(
                seismograms[i], srate_grn, freq_band, butter_order, zero_phase
            )
    
    ts_count = 0
    if before_p is not None:
        ts_count = round((green_before_p - before_p) * srate_grn)
    if pad_zeros:
        if before_p is not None:
            raise ValueError("can not set before_p and pad_zeros together")
        ts_count = round((green_before_p - tpts_table["p_onset"]) * srate_grn)
        before_p = tpts_table["p_onset"]
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
            green_before_p=before_p,
            event_depth_km=event_depth_km,
            dist_in_km=dist_km,
            model_name=model_name,
        )

    len_after_resample = round(sampling_num * srate / srate_grn)
    seismograms_resample = np.zeros((3, len_after_resample))
    for i in range(3):
        seismograms_resample[i] = resample(
            seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
        )[:len_after_resample]

    if output_type == "disp":
        seismograms_resample = np.cumsum(seismograms_resample, axis=1) / srate
    elif output_type == "acce":
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
