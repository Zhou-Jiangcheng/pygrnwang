import os
import json

import numpy as np
import pandas as pd

from .create_qssp import mt_com_list
from .utils import (
    shift_green2real_tpts,
    create_rotate_z_mat,
    rotate_symmetric_tensor_series,
)
from .focal_mechanism import convert_mt_axis, mt2full_mt_matrix, check_convert_fm
from .geo import rotate_rtz_to_enz
from .signal_process import resample

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


def read_tpts_table_qssp2020(path_green, event_depth_km, receiver_depth_km, ind):
    fr_tp = open(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "tp_table.bin",
        ),
        "rb",
    )
    tp = np.fromfile(file=fr_tp, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_tp.close()

    fr_ts = open(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "ts_table.bin",
        ),
        "rb",
    )
    ts = np.fromfile(file=fr_ts, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_ts.close()
    return float(tp), float(ts)


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
    :param path_green:
    :param event_depth: km
    :param receiver_depth: km
    :param dist: km
    :param output_type: str, 'output_type must in  disp | velo | acce | strain | strain_rate | '
            'stress | stress_rate | rotation | rotation_rate | gravitation | gravimeter'
    :param green_info:
    """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    N_T = round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    if output_type == "gravimeter":
        enz_list = [""]
        data_all = np.zeros((N_T, 6))
    elif output_type in ["disp", "velo", "acce", "rota", "rota_rate", "gravitation"]:
        enz_list = ["_e", "_n", "_z"]
        data_all = np.zeros((N_T, 18))
    elif output_type in ["stress", "stress_rate", "strain", "strain_rate"]:
        enz_list = ["_ee", "_en", "_ez", "_nn", "_nz", "_zz"]
        data_all = np.zeros((N_T, 36))
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
            dist_range = green_info["dist_range"]
            delta_dist = green_info["delta_dist"]
            ind = round((dist - dist_range[0]) / delta_dist)
            path_bin = os.path.join(
                path_func, "_%s%s.bin" % (output_type, enz_list[i_enz])
            )
            # path_bin = ''
            if os.path.exists(path_bin):
                data_enz = read_time_series_qssp2020(
                    path_bin=path_bin, ind=ind, sampling_num=N_T
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


def seek_qssp2020(
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
    rotate=True,
    only_seismograms=True,
    output_type="disp",
    model_name="ak135",
    green_info=None,
):
    """
    Read seismic data with either one, three, or six components based on output_type.

    Depending on output_type, the function reads:
      - one component (e.g., gravimeter),
      - three components (e.g., displacement, velocity, acceleration, rotation, rotation rate, gravitation),
      - or six components (e.g., strain, strain rate, stress, stress rate).

    Parameters:
        path_green (str): Root directory of the data.
        event_depth_km (float): Event depth in km.
        receiver_depth_km (float): Receiver depth in km.
        az_deg (float): Azimuth in degrees.
        dist_km (float): Epicentral distance in km.
        focal_mechanism (ndarray/list):
                strike, dip, rake(deg), length=3 or
                [M11, M12, M13, M22, M23, M33], in NED coordinates, length=6.
        srate (float): Sampling rate in Hz.
        before_p (float, optional): Time before the P-wave arrival in seconds.
        pad_zeros (bool, optional): Whether to pad with zeros.
        shift (bool, optional): Whether to shift the seismograms.
        rotate (bool, optional): Whether to perform rotation from rtz2ned.
        only_seismograms (bool, optional): If True, only return seismograms.
        output_type (str, optional): Output type string. For example, 'disp','velo','acce',
                                     'rota','rota_rate','gravimeter','gravitation',
                                     'strain','strain_rate','stress','stress_rate'.
        model_name (str, optional): Model name string. If not provided, a default is set based on output_type:
                                    - For one or three component types: "ak135"
                                    - For six component types: "ak135" (or adjust as needed)
        green_info (dict, optional): Information of Green's function libarary
    Returns:
        If only_seismograms is True:
            seismograms (ndarray): For one-component: (1 x N_T),
                                   for three-component: (3 x N_T),
                                   for six-component: (6 x N_T).
        Otherwise, returns a tuple:
            (seismograms, tpts_table, first_p, first_s, green_dist)
            where tpts_table is a dict with keys "p_onset" and "s_onset".
    """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    srate_grn = 1 / green_info["sampling_interval"]
    sampling_num = (
        round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    )
    dist_range = green_info["dist_range"]
    delta_dist = green_info["delta_dist"]
    time_reduction = green_info["time_reduction"]
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]
    if not isinstance(grn_dep_list, list):
        grn_dep = grn_dep_list
    else:
        grn_dep = grn_dep_list[
            np.argmin(np.abs(event_depth_km - np.array(grn_dep_list)))
        ]
    if not isinstance(grn_receiver_list, list):
        grn_receiver = grn_receiver_list
    else:
        grn_receiver = grn_receiver_list[
            np.argmin(np.abs(receiver_depth_km - np.array(grn_receiver_list)))
        ]

    # Common: calculate rotation matrix and convert moment tensor
    gamma = np.deg2rad(az_deg)
    A_rotate = create_rotate_z_mat(gamma=gamma)
    focal_mechanism = check_convert_fm(focal_mechanism)
    mt_ned_full = mt2full_mt_matrix(mt=focal_mechanism, flag="ned")
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

    # Select processing branch based on output_type
    if output_type in one_com_list:
        # One-component branch (e.g., gravimeter)
        u = seek_raw_qssp2020(
            path_green, grn_dep, grn_receiver, dist_km, output_type, green_info
        )
        u_enz_green_north = np.zeros((sampling_num, 1))
        for i_rtp in range(6):
            u_enz_green_north[:, 0] += mt_rtp[i_rtp] * u[:, i_rtp]
        seismograms = u_enz_green_north.T

    elif output_type in three_com_list:
        # Three-component branch: process data (e.g., displacement, velocity, acceleration)
        u_all = seek_raw_qssp2020(
            path_green, grn_dep, grn_receiver, dist_km, output_type, green_info
        )
        u_enz_green_north = np.zeros((sampling_num, 3))
        for i_rtp in range(6):
            for i_enz in range(3):
                u_enz_green_north[:, i_enz] += (
                    mt_rtp[i_rtp] * u_all[:, i_enz + 3 * i_rtp]
                )
        # Convert from green_north to RTZ coordinates
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
        # Six-component branch: process data (e.g., strain, strain rate, stress, stress rate)
        epsilon_all = seek_raw_qssp2020(
            path_green, grn_dep, grn_receiver, dist_km, output_type, green_info
        )
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

    # Compute the Green's function distance closest to the input distance.
    nearest_indice = round((dist_km - dist_range[0]) / delta_dist)
    green_dist = dist_range[0] + nearest_indice * delta_dist

    # Process P-wave arrival times if before_p, shift, or pad_zeros is set.
    tpts_table = None
    if (before_p is not None) or shift or pad_zeros:
        grn_first_p, grn_first_s = read_tpts_table_qssp2020(
            path_green=path_green,
            event_depth_km=grn_dep,
            receiver_depth_km=grn_receiver,
            ind=nearest_indice,
        )
        tpts_table = {"p_onset": grn_first_p, "s_onset": grn_first_s}

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

    # Shift the seismograms if required.
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
            receiver_depth_km=receiver_depth_km,
            model_name=model_name,
        )

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

    # Return results.
    if only_seismograms:
        return seismograms_resample
    else:
        return (
            seismograms_resample,
            tpts_table,
            first_p,
            first_s,
            grn_dep,
            grn_receiver,
            green_dist,
        )


if __name__ == "__main__":
    pass
