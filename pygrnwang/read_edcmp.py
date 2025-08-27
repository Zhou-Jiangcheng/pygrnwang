import os
import json
import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

from .geo import d2km


def read_edcmp_raw(path_green, output_type, obs_depth, obs_depth_list):
    grn_obs_depth = obs_depth_list[
        np.argmin(np.abs(obs_depth - np.array(obs_depth_list)))
    ]
    df = pd.read_csv(
        str(
            os.path.join(
                path_green, "edcmp2", "%.2f" % grn_obs_depth, "hs.%s" % output_type
            )
        ),
        skiprows=3,
        sep="\\s+",
        header=None,
    )
    # Ux_m, Uy_m, Uz_m
    # Sxx_Pa Syy_Pa Szz_Pa Sxy_Pa Syz_Pa Szx_Pa
    values_raw = df.to_numpy()[:, 2:]
    return values_raw


def interpolate_values(xmin, xmax, nx, ymin, ymax, ny, v_array, obs_array):
    x = np.linspace(xmin, xmax, nx)
    y = np.linspace(ymin, ymax, ny)
    v_grid = v_array.reshape(ny, nx, v_array.shape[1])
    interpolator = RegularGridInterpolator((y, x), v_grid)
    pts = obs_array[:, [1, 0]]
    sigma_interp = interpolator(pts)
    return sigma_interp


def seek_edcmp2(
    path_green: str,
    output_type: str,
    obs_array: np.ndarray,
    geo_coordinate=True,
):
    """
    :param path_green: the root dir of Green's function lib
    :param output_type: 'disp','strain','stress','tilt'
    :param obs_array:
      Cartesian coordinate, [[x_1, y_1, z_1]], [x_2, y_2, z_2],..., [x_n, y_n, z_n]], or
      Geology coordinate, [[lat_1, lon_1, dep_1], [lat_2, lon_2, dep_2],
                          ..., [lat_n, lon_n, dep_n]]
      shape is (n, 3)
    :param geo_coordinate: True/False, use Geology/Cartesian coordinate

    :return values at obs_array
    """
    obs_array = obs_array.copy()
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)
    obs_x_range = green_info["obs_x_range"]
    obs_y_range = green_info["obs_y_range"]
    obs_delta_x = green_info["obs_delta_x"]
    obs_delta_y = green_info["obs_delta_y"]
    obs_depth_list = green_info["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]
    nx = math.ceil((obs_x_range[1] - obs_x_range[0]) / obs_delta_x) + 1
    ny = math.ceil((obs_y_range[1] - obs_y_range[0]) / obs_delta_y) + 1
    x_min, x_max = obs_x_range
    y_min, y_max = obs_y_range

    if output_type == "disp":
        cha_num = 3
    elif output_type == "strain":
        cha_num = 6
    elif output_type == "stress":
        cha_num = 6
    elif output_type == "tilt":
        cha_num = 2
    else:
        raise ValueError("output_type must in disp,strain,stress,tilt")
    values_output = np.zeros((obs_array.shape[0], cha_num)) + np.nan
    if geo_coordinate:
        obs_ref = green_info["obs_ref"]
        if obs_ref is None:
            raise ValueError("The Green's function lib dose not use obs_ref!")
        elif (
            np.min(obs_array[:, 0]) - obs_ref[0] < x_min / d2km
            or np.max(obs_array[:, 0]) - obs_ref[0] > x_max / d2km
            or np.min(obs_array[:, 1]) - obs_ref[1] < y_min / d2km
            or np.max(obs_array[:, 1]) - obs_ref[1] > y_max / d2km
            # or np.min(obs_array[:, 2]) < obs_depth_list[0]
            # or np.max(obs_array[:, 2]) > obs_depth_list[-1]
        ):
            print(
                np.min(obs_array[:, 0]) - obs_ref[0] < x_min / d2km,
                np.max(obs_array[:, 0]) - obs_ref[0] > x_max / d2km,
                np.min(obs_array[:, 1]) - obs_ref[1] < y_min / d2km,
                np.max(obs_array[:, 1]) - obs_ref[1] > y_max / d2km,
                # np.min(obs_array[:, 2]) < obs_depth_list[0],
                # np.max(obs_array[:, 2]) > obs_depth_list[-1]
            )
            raise ValueError("obs_array exceeds the range of Green's function lib!")
        else:
            obs_array[:, 0] = (obs_array[:, 0] - obs_ref[0]) * d2km
            obs_array[:, 1] = (obs_array[:, 1] - obs_ref[1]) * d2km

    unique_depths = np.unique(obs_array[:, 2])
    for i in range(len(unique_depths)):
        inds_i = np.where(obs_array[:, 2] == unique_depths[i])
        obs_array_dep_i = obs_array[inds_i]
        v_dep_i = read_edcmp_raw(
            path_green, output_type, unique_depths[i], obs_depth_list
        )
        v_interp_dep_i = interpolate_values(
            x_min, x_max, nx, y_min, y_max, ny, v_dep_i, obs_array_dep_i
        )
        values_output[inds_i] = v_interp_dep_i
    return values_output
