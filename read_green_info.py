import os
import numpy as np

from pygrnwang.utils import get_number_in_line


def read_green_info_spgrn2020(path_greenfunc: str, green_depth: float) -> dict:
    """
    read GreenInfo.dat
    :param path_greenfunc:
    :param green_depth: unit km
    :return: {
             "time_window": time_window, type float, unit s
             "sampling_interval": sampling_interval, type float, unit s
             "samples_num": int(samples_num), type int
             "dist_list": dist_list, type list, unit km
         }
    """
    with open(
        os.path.join(path_greenfunc, "GreenInfo%.1f.dat" % green_depth), "r"
    ) as fr:
        lines = fr.readlines()
    [time_window, sampling_interval,
        samples_num] = get_number_in_line(lines[6])
    number_of_distance = get_number_in_line(lines[9])[0]
    dist_list = []
    for i in range(int(np.ceil(number_of_distance / 5))):
        temp = get_number_in_line(lines[10 + i])
        for item in temp:
            dist_list.append(item)
    return {
        "time_window": time_window,
        "sampling_interval": sampling_interval,
        "samples_num": int(samples_num),
        "dist_list": dist_list,
    }


def read_green_info_qseis06(path_greenfunc: str, green_depth: float) -> dict:
    """
    read GreenInfo.dat
    :param path_greenfunc:
    :param green_depth: unit km
    :return: {
             "time_window": time_window, type float, unit s
             "sampling_interval": sampling_interval, type float, unit s
             "samples_num": int(samples_num), type int
             "dist_range": [dist_min, dist_max], type list, unit deg
             "delta_dist": delta_dist, type float, unit deg
             "N_dist": len(dist_list), type int
             "N_dist_group": type int
         }
    """
    with open(
        os.path.join(path_greenfunc, "GreenInfo%.1f.dat" % green_depth), "r"
    ) as fr:
        lines = fr.readlines()
    time_window = get_number_in_line(lines[1])[0]
    sampling_num = int(get_number_in_line(lines[2])[0])
    sampling_interval = get_number_in_line(lines[3])[0]
    dist_range = get_number_in_line(lines[4])
    delta_dist = get_number_in_line(lines[5])[0]
    N_dist = int(get_number_in_line(lines[6])[0])
    N_dist_group = int(get_number_in_line(lines[7])[0])
    N_each_group = int(get_number_in_line(lines[8])[0])
    return {
        "time_window": time_window,
        "sampling_num": sampling_num,
        "sampling_interval": sampling_interval,
        "dist_range": dist_range,
        "delta_dist": delta_dist,
        "N_dist": N_dist,
        "N_dist_group": N_dist_group,
        "N_each_group": N_each_group,
    }
