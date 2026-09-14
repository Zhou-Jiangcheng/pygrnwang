import os

import numpy as np


def read_green_info_spgrn(path_greenfunc: str, green_depth: float) -> dict:
    """Read native SPGRN time sampling and distance nodes.

    Parameters
    ----------
    path_greenfunc : str
        SPGRN source/receiver depth folder containing GreenInfo and P/S table files.
    green_depth : float
        Exact source depth in km used in the GreenInfo filename.

    Returns
    -------
    info : dict
        time_window and sampling_interval in seconds, integer samples_num,
        and dist_list in km.

    Raises
    ------
    OSError
        The GreenInfo file cannot be read.
    ValueError
        The file does not follow the expected backend format.
    """
    with open(
        os.path.join(path_greenfunc, "GreenInfo%.2f.dat" % green_depth), "r"
    ) as fr:
        lines = fr.readlines()
    [time_window, sampling_interval, samples_num] = lines[6].strip().split()
    time_window = float(time_window)
    sampling_interval = float(sampling_interval)
    samples_num = int(samples_num)
    number_of_distance = float(lines[9].strip())
    dist_list = []
    for i in range(round(np.ceil(number_of_distance / 5))):
        temp = lines[10 + i].strip().split()
        for item in temp:
            dist_list.append(float(item))
    return {
        "time_window": time_window,
        "sampling_interval": sampling_interval,
        "samples_num": int(samples_num),
        "dist_list": dist_list,
    }
