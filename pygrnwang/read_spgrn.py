import os
import json

import numpy as np
from scipy import signal

from .read_tpts_table import read_tpts_table
from .utils import shift_green2real_tpts
from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz
from .signal_process import resample


def synthesize_spgrn2020(az_in_deg, time_series, focal_mechanism):
    """
    Kanamori, H., & Rivera, L. (2008).
    Source inversion of W phase: speeding up seismic tsunami warning.
    Geophysical Journal International, 175(1), 222–238.
    https://doi.org/10.1111/j.1365-246X.2008.03887.x
    """
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


def read_time_series_spgrn2020(path_grn_data, dist_in_km, green_info):
    # dist_green, dist_green_num = find_nearest_dichotomy(
    #     dist_in_km, green_info["dist_list"]
    # )
    dist_green_num = np.argmin(np.abs(np.array(green_info["dist_list"]) - dist_in_km))
    dist_green = green_info["dist_list"][dist_green_num]
    N_T = round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    length_each = 3 + (2 + N_T) * 10
    start_count = dist_green_num * length_each

    fr = open(path_grn_data, "rb")
    time_series = np.fromfile(
        file=fr, dtype=np.float32, count=length_each, offset=start_count * 4
    )
    time_series = time_series[3:].reshape(10, 2 + N_T)
    time_series = time_series[:, 1:-1]
    fr.close()
    # for i in range(10):
    #     plt.figure()
    #     plt.plot(time_series[i][:500])
    #     plt.show()
    return time_series, dist_green


def seek_spgrn2020(
    path_green,
    output_type,
    event_depth_km,
    receiver_depth_km,
    az_deg,
    dist_km,
    focal_mechanism,
    srate,
    rotate=True,
    before_p=20,
    pad_zeros=False,
    shift=False,
    only_seismograms=True,
    model_name="ak135",
    green_info=None,
):
    """

    :param path_green:
    :param output_type: 'disp'/'velo'/'acce'
    :param event_depth_km:
    :param receiver_depth_km:
    :param az_deg:
    :param dist_km:
    :param focal_mechanism:
    :param srate:
    :param rotate:
    :param before_p:
    :param pad_zeros:
    :param shift:
    :param only_seismograms: only return seismograms
    :param model_name: earth model name in obspy TaupModel
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
        os.path.join(
            path_green, "GreenFunc", "%.2f" % grn_dep_source, "%.2f" % grn_dep_receiver
        )
    )

    # green_info = read_green_info_spgrn2020(
    #     path_greenfunc=path_greenfunc, green_depth=grn_dep_source
    # )
    tpts_table = read_tpts_table(
        path_greenfunc=path_greenfunc, dist_in_km=dist_km, green_info=green_info
    )

    path_grn_data = os.path.join(path_greenfunc, "grn_d%.2f" % grn_dep_source)
    time_series, grn_dist = read_time_series_spgrn2020(
        path_grn_data=path_grn_data, dist_in_km=dist_km, green_info=green_info
    )
    # r,t,z
    seismograms = synthesize_spgrn2020(
        az_in_deg=az_deg, time_series=time_series, focal_mechanism=focal_mechanism
    )
    if rotate:
        seismograms = rotate_rtz_to_enz(
            az_in_deg=az_deg, r=seismograms[0], t=seismograms[1], z=seismograms[2]
        )[:]

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

    # conv_shift = round(green_info["source_duration"] * srate_grn / 2)
    # if conv_shift != 0:
    #     seismograms = np.roll(seismograms, -conv_shift)
    #     seismograms[:, -conv_shift:] = 0

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
