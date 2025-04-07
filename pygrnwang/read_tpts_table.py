import os
import struct

import numpy as np


def read_tpts_table(path_greenfunc: str, dist_in_km: float, green_info: dict) -> dict:
    """

    :param path_greenfunc:
    :param dist_in_km:
    :param green_info:
    :return:  {
        "p_onset": tp_onset,
        "p_takeoff": tp_takeoff,
        "p_slowness": tp_slowness,
        "s_onset": ts_onset,
        "s_takeoff": ts_takeoff,
        "s_slowness": ts_slowness
    }
    """
    dist_green_num = np.argmin(np.abs(np.array(green_info["dist_list"]) - dist_in_km))
    # print(dist_green_num)
    # print(green_info['dist_list'][dist_green_num])
    start_count = 4 + dist_green_num * 12

    def seek_time_table(start_count_in, phase):
        if phase == "P":
            fname = "tptable.dat"
        elif phase == "S":
            fname = "tstable.dat"
        else:
            raise ValueError("phase must be P or S")
        fr = open(os.path.join(path_greenfunc, fname), "rb")
        fr.seek(start_count_in)
        t_onset = struct.unpack("f", fr.read(4))[0]
        t_takeoff = struct.unpack("f", fr.read(4))[0]
        t_slowness = struct.unpack("f", fr.read(4))[0]
        fr.close()
        return [t_onset, t_takeoff, t_slowness]

    [tp_onset, tp_takeoff, tp_slowness] = seek_time_table(start_count, "P")
    [ts_onset, ts_takeoff, ts_slowness] = seek_time_table(start_count, "S")

    return {
        "p_onset": tp_onset,
        "p_takeoff": tp_takeoff,
        "p_slowness": tp_slowness,
        "s_onset": ts_onset,
        "s_takeoff": ts_takeoff,
        "s_slowness": ts_slowness,
    }
