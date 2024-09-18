import os

import matplotlib.pyplot as plt
import numpy as np

from pygrnwang.read_green_info import read_green_info_spgrn2020
from pygrnwang.read_tpts_table import read_tpts_table
from pygrnwang.utils import *
from pygrnwang.signal_process import linear_interp, resample

d2km = 111.19492664455873


def find_green_depth_list(path_greenfunc):
    green_depth_list: list = os.listdir(path_greenfunc)
    for i in range(len(green_depth_list)):
        green_depth_list[i] = float(green_depth_list[i])
    green_depth_list.sort()
    return green_depth_list


def read_time_series(path_grn_data, dist_in_km, green_info):
    dist_green, dist_green_num = find_nearest_dichotomy(
        dist_in_km, green_info["dist_list"]
    )
    length_each = 3 + (1 + green_info["samples_num"] + 1) * 10
    start_count = dist_green_num * length_each

    fr = open(path_grn_data, "rb")
    time_series = np.fromfile(
        file=fr, dtype=np.float32, count=length_each, sep="", offset=start_count * 4
    )
    time_series = time_series[3:].reshape(10, green_info["samples_num"] + 2)
    time_series = time_series[:, 1:-1]
    fr.close()
    # for i in range(10):
    #     plt.figure()
    #     plt.plot(time_series[i][:500])
    #     plt.show()
    return time_series, dist_green


def synthesize(az_in_deg, time_series, focal_mechanism):
    [M11, M12, M13, M22, M23, M33] = check_convert_fm(focal_mechanism=focal_mechanism)
    # expl.,strike-slip,dip-slip,clvd
    exp = (M11 + M22 + M33) / 3
    clvd = M33 - exp
    ss1 = -M12
    ss2 = (M11 - M22) / 2
    ds1 = M13
    ds2 = -M23

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
    t = time_series[4] * m2[0] + time_series[7] * m2[1]
    r = (
        time_series[1] * m1[0]
        + time_series[3] * m1[1]
        + time_series[6] * m1[2]
        + time_series[9] * m1[3]
    )

    seismograms = [z, t, r]
    return seismograms


def read_spgrn(
    path_greenfunc,
    green_depth_in_km=None,
    az_in_deg=None,
    dist_in_km=None,
    focal_mechanism=None,
    srate=2,
    zero_phase=False,
    rotate=True,
    green_info=None,
    only_seismograms=True,
):
    """

    :param path_greenfunc:
    :param green_depth_in_km:
    :param az_in_deg:
    :param dist_in_km:
    :param focal_mechanism:
    :param srate:
    :param zero_phase:
    :param rotate:
    :param green_info: dict, returned by read_green_info_spgrn2020
    :param only_seismograms
    :return:
    """
    if green_info is None:
        green_info = read_green_info_spgrn2020(
            path_greenfunc=path_greenfunc, green_depth=green_depth_in_km
        )
    path_grn_data = os.path.join(path_greenfunc, "grn_d%.1f" % green_depth_in_km)
    time_series, green_dist = read_time_series(
        path_grn_data=path_grn_data, dist_in_km=dist_in_km, green_info=green_info
    )
    # with open('time_series.npy', 'wb') as fw:
    #     np.save(fw, time_series)
    # z,t,r
    seismograms = synthesize(
        az_in_deg=az_in_deg, time_series=time_series, focal_mechanism=focal_mechanism
    )
    if rotate:
        seismograms = rotate_rtz_to_enz(
            az_in_deg=az_in_deg, r=seismograms[2], t=seismograms[1], z=seismograms[0]
        )[:]
    for i in range(3):
        srate_old = 1 / green_info["sampling_interval"]
        seismograms[i] = resample(
            seismograms[i], srate_old=srate_old, srate_new=srate, zero_phase=zero_phase
        )
    if only_seismograms:
        return seismograms
    else:
        return seismograms, green_dist


def seek_spgrn(
    path_green_lib,
    event_depth_in_km=None,
    az_in_deg=None,
    dist_in_km=None,
    focal_mechanism=None,
    srate=2,
    zero_phase=False,
    rotate=True,
    green_before_p=None,
    only_seismograms=True,
):
    """

    :param path_green_lib:
    :param event_depth_in_km:
    :param az_in_deg:
    :param dist_in_km:
    :param focal_mechanism:
    :param srate:
    :param zero_phase:
    :param rotate:
    :param only_seismograms: only return seismograms
    :return:
    """
    grn_dep_list = find_green_depth_list(
        path_greenfunc=os.path.join(path_green_lib, "GreenFunc")
    )
    grn_dep = find_nearest_dichotomy(value=event_depth_in_km, value_list=grn_dep_list)[
        0
    ]
    path_greenfunc = str(os.path.join(path_green_lib, "GreenFunc", "%.1f" % grn_dep))

    green_info = read_green_info_spgrn2020(
        path_greenfunc=path_greenfunc, green_depth=grn_dep
    )
    tpts_table = read_tpts_table(
        path_greenfunc=path_greenfunc, dist_in_km=dist_in_km, green_info=green_info
    )
    seismograms, green_dist = read_spgrn(
        path_greenfunc=path_greenfunc,
        green_depth_in_km=grn_dep,
        az_in_deg=az_in_deg,
        dist_in_km=dist_in_km,
        focal_mechanism=focal_mechanism,
        srate=srate,
        zero_phase=zero_phase,
        rotate=rotate,
        green_info=green_info,
        only_seismograms=False,
    )
    if only_seismograms:
        return seismograms
    else:
        return seismograms, tpts_table, grn_dep, green_dist


if __name__ == "__main__":
    pass
