import os
import json
import math

import numpy as np
from scipy.interpolate import RegularGridInterpolator
import pandas as pd

from .focal_mechanism import (
    check_convert_fm,
    tensor2full_tensor_matrix,
    plane2mt,
    mt2plane,
)
from .utils import create_rotate_z_mat, read_material_nd
from .geo import rotate_rtz_to_enz


def read_edcmp_raw(path_green, output_type, grn_event_depth, grn_obs_depth, mt_ind):
    path_npy = str(
        os.path.join(
            path_green,
            "edcmp2",
            "%.2f" % grn_event_depth,
            "%.2f" % grn_obs_depth,
            "%d" % mt_ind,
            "%s.npy" % output_type,
        )
    )
    if os.path.exists(path_npy):
        values_raw = np.load(path_npy)
    else:
        df = pd.read_csv(
            str(
                os.path.join(
                    path_green,
                    "edcmp2",
                    "%.2f" % grn_event_depth,
                    "%.2f" % grn_obs_depth,
                    "%d" % mt_ind,
                    "hs.%s" % output_type,
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
    event_depth_km,
    receiver_depth_km,
    az_deg,
    dist_km,
    focal_mechanism,
    rotate=True,
    check_convert_pure_dp=True,
    output_type: str = "disp",
    times_mu: bool = False,
    model_name="ak135",
    green_info=None,
):
    """
    ee en ez nn nz zz
    :param path_green: the root dir of Green's function lib
    :param output_type: 'disp','strain','stress','tilt','mu_disp'
    :param times_mu: Whether to times mu (=rho*beta**2) in the result
    """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    grn_source_depth_range = green_info["grn_source_depth_range"]
    grn_source_delta_depth = green_info["grn_source_delta_depth"]
    event_depth_list = np.arange(
        grn_source_depth_range[0],
        grn_source_depth_range[1] + grn_source_delta_depth,
        grn_source_delta_depth,
    )
    obs_depth_list = green_info["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]
    obs_depth_list = np.array(obs_depth_list)

    grn_dist_range = green_info["grn_dist_range"]
    grn_dist_delta = green_info["grn_delta_dist"]
    dist_list = np.arange(
        grn_dist_range[0], grn_dist_range[1] + grn_dist_delta, grn_dist_delta
    )

    grn_event_depth = event_depth_list[
        np.argmin(np.abs(event_depth_list - event_depth_km))
    ]
    grn_obs_depth = obs_depth_list[
        np.argmin(np.abs(obs_depth_list - receiver_depth_km))
    ]
    grn_dist_ind = np.argmin(np.abs(dist_list - dist_km))

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

    A_rotate = create_rotate_z_mat(gamma=np.deg2rad(az_deg))
    focal_mechanism = check_convert_fm(focal_mechanism)
    mt_ned_full = tensor2full_tensor_matrix(mt=focal_mechanism, flag="ned")
    mt_rotate = A_rotate.T @ mt_ned_full @ A_rotate
    mt = np.array(
        [
            mt_rotate[0, 0],
            mt_rotate[0, 1],
            mt_rotate[0, 2],
            mt_rotate[1, 1],
            mt_rotate[1, 2],
            mt_rotate[2, 2],
        ]
    )
    if check_convert_pure_dp:
        mt_dp = plane2mt(1, *mt2plane(mt)[0])
    else:
        mt_dp = mt
    v_ned_green_north = np.zeros(cha_num, dtype=float)

    for i in range(5):
        v_raw = read_edcmp_raw(
            path_green=path_green,
            output_type=output_type,
            grn_event_depth=grn_event_depth,
            grn_obs_depth=grn_obs_depth,
            mt_ind=i,
        )
        v_ned_green_north = v_ned_green_north + v_raw[grn_dist_ind] * mt_dp[i]
    v_raw = read_edcmp_raw(
        path_green=path_green,
        output_type=output_type,
        grn_event_depth=grn_event_depth,
        grn_obs_depth=grn_obs_depth,
        mt_ind=3,
    )
    v_ned_green_north = v_ned_green_north + v_raw[grn_dist_ind] * mt_dp[0]

    v_rtz = np.zeros_like(v_ned_green_north)
    v_rtz[0] = v_ned_green_north[0]
    v_rtz[1] = -v_ned_green_north[1]
    v_rtz[2] = -v_ned_green_north[2]

    if rotate:
        v = rotate_rtz_to_enz(az_in_deg=az_deg, r=v_rtz[0], t=v_rtz[1], z=v_rtz[2])
    else:
        v = v_rtz

    if not times_mu:
        material_row = read_material_nd(model_name=model_name, depth=grn_event_depth)
        beta = material_row[2]
        rho = material_row[3]
        mu_pa = rho * beta**2 * 1e9
        v = v / mu_pa
    return v
