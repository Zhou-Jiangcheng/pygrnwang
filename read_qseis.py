import os
import numpy as np

from pygrnwang.read_green_info import read_green_info_qseis06
from pygrnwang.read_spgrn import cal_first_p_s, shift_green2real_tpts
from pygrnwang.utils import find_nearest_dichotomy, check_convert_fm
from pygrnwang.signal_process import linear_interp, resample
d2km = 111.19492664455873


def find_green_depth_list(path_green_lib):
    fname_list: list = os.listdir(path_green_lib)
    depth_list = []
    for fname in fname_list:
        try:
            dep = float(fname)
            depth_list.append(dep)
        except ValueError:
            pass
    depth_list.sort()
    return depth_list


def find_ind(dist_in_deg, greeninfo, num_each_group=100):
    dist_range = greeninfo["dist_range"]
    delta_dist = greeninfo["delta_dist"]
    ind = round((dist_in_deg - dist_range[0]) / delta_dist)
    ind_group = ind // num_each_group
    green_dist = dist_range[0] + ind * delta_dist
    sampling_num = int(greeninfo["sampling_num"])
    start_count = (
        (ind + 1 - ind_group * num_each_group) * sampling_num * 10
    )  # +1 because of T_sec line
    return ind, ind_group, green_dist, start_count, sampling_num


def read_time_series(path_greenfunc, start_count, sampling_num):
    fr = open(os.path.join(path_greenfunc, "grn.dat"), "rb")
    time_series = np.fromfile(
        file=fr, dtype=np.float32, count=sampling_num * 10, offset=start_count * 4
    )
    time_series = time_series.reshape(10, sampling_num)
    fr.close()
    return time_series


def synthesize(az_in_deg, time_series, focal_mechanism):
    """
    z upward, t north-east, r source-station
    return [z,t,r]
    """
    [M11, M12, M13, M22, M23, M33] = check_convert_fm(
        focal_mechanism=focal_mechanism)
    exp = (M11 + M22 + M33) / 3
    clvd = (-0.5 * M11 - 0.5 * M22 + M33) / 2
    ss1 = M12
    ss2 = (M11 - M22) / 2
    ds1 = M13
    ds2 = M23

    az = np.deg2rad(az_in_deg)
    sin_az, cos_az = np.sin(az), np.cos(az)
    sin_2az, cos_2az = np.sin(2 * az), np.cos(2 * az)
    m1 = [exp, ss1 * sin_2az + ss2 * cos_2az,
          ds1 * cos_az + ds2 * sin_az, clvd]
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

    seismograms = [-z, t, r]
    return seismograms


def rotate_rtz_to_enz(az_in_deg, r, t, z):
    az = np.deg2rad(az_in_deg)
    e = r * np.sin(az) + t * np.cos(az)
    n = r * np.cos(az) - t * np.sin(az)
    return [e, n, z]


def seek_qseis(
    path_green_lib,
    event_depth_in_km=None,
    az_in_deg=None,
    dist_in_km=None,
    focal_mechanism=None,
    srate=1.0,
    zero_phase=False,
    rotate=True,
    time_reduction_slowness=8,
    before_p=None,
    pad_zeros=False,
    shift=False,
    only_seismograms=True,
    model_name="ak135fc",
):
    grn_dep_list = find_green_depth_list(path_green_lib)
    grn_dep = find_nearest_dichotomy(event_depth_in_km, grn_dep_list)[0]
    path_greenfunc = os.path.join(path_green_lib, "%.1f" % grn_dep)
    green_info = read_green_info_qseis06(
        path_greenfunc=path_greenfunc, green_depth=grn_dep
    )
    ind, ind_group, green_dist, start_count, sampling_num = find_ind(
        dist_in_deg=dist_in_km/d2km,
        greeninfo=green_info,
        num_each_group=green_info["N_each_group"],
    )
    path_greenfunc = os.path.join(path_greenfunc, "%d" % ind_group)
    # print(path_greenfunc, sampling_num)
    time_series = read_time_series(
        path_greenfunc=path_greenfunc,
        start_count=start_count,
        sampling_num=sampling_num,
    )
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
    first_p = None
    tpts_table = None
    if (before_p is not None) or shift or pad_zeros:
        first_p, first_s = cal_first_p_s(
            event_depth=event_depth_in_km,
            dist_in_km=green_dist*d2km,
            model_name=model_name,
        )
        tpts_table = {"p_onset": first_p, "s_onset": first_s}

    if before_p is None:
        pass
    else:
        tp = first_p - green_dist * time_reduction_slowness
        ts_count = round((tp - before_p) * srate)
        # print(first_p, tp, ts_count)
        if ts_count >= 0:
            for i in range(3):
                seismograms[i] = seismograms[i][ts_count:]
        else:
            for i in range(3):
                seismograms[i] = np.concatenate(
                    [np.zeros(-ts_count), seismograms[i]])

    if shift:
        seismograms, first_p, first_s = shift_green2real_tpts(
            seismograms=seismograms,
            tpts_table=tpts_table,
            srate=srate,
            before_p=before_p,
            event_depth_in_km=event_depth_in_km,
            dist_in_km=dist_in_km,
            model_name=model_name,
        )

    if pad_zeros:
        if before_p:
            raise ValueError("can not set before_p and pad_zeros together")
        for i in range(3):
            seismograms[i] = np.concatenate(
                [np.zeros(round((first_p - before_p) * srate)), seismograms[i]]
            )

    if only_seismograms:
        return seismograms
    else:
        return seismograms, green_dist


if __name__ == "__main__":
    pass
