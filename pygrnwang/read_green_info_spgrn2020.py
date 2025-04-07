import os
import ast

import numpy as np


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
        os.path.join(path_greenfunc, "GreenInfo%.2f.dat" % green_depth), "r"
    ) as fr:
        lines = fr.readlines()
    [time_window, sampling_interval, samples_num] = ast.literal_eval(lines[6].strip())
    number_of_distance = ast.literal_eval(lines[9].strip())[0]
    dist_list = []
    for i in range(round(np.ceil(number_of_distance / 5))):
        temp = ast.literal_eval(lines[10 + i].strip())
        for item in temp:
            dist_list.append(item)
    return {
        "time_window": time_window,
        "sampling_interval": sampling_interval,
        "samples_num": int(samples_num),
        "dist_list": dist_list,
    }
