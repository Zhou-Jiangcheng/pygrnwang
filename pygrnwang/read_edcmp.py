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
from .utils import create_rotate_z_mat, read_material_nd, read_nd
from .geo import rotate_rtz_to_enz, rotate_symmetric_tensor_series


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
    area_km_sq: float = None,
    model_name="ak135",
    green_info=None,
):
    """
    Read static deformation from edcmp2 Green's function library for a single query point.

    Uses nearest-neighbor lookup in source depth, receiver depth, and epicentral distance.
    The result is synthesized from 5 elementary MT components stored per grid node.

    Output components for each output_type:
      disp   (3): east, north, up  [m / (N·m)]
      strain (6): ee, en, ez, nn, nz, zz  [1 / (N·m)]
      stress (6): ee, en, ez, nn, nz, zz  [Pa / (N·m)]
      tilt   (2): east tilt, north tilt  [rad / (N·m)]
    (all divided by mu = rho*beta^2 unless times_mu=True)

    :param path_green: Root directory of the Green's function library.
    :param event_depth_km: Event (source) depth in km.
    :param receiver_depth_km: Receiver depth in km.
    :param az_deg: Azimuth from source to receiver in degrees (measured from north).
    :param dist_km: Epicentral distance in km.
    :param focal_mechanism: [strike, dip, rake] in degrees, or
                            [M11, M12, M13, M22, M23, M33] moment tensor in N·m.
    :param rotate: If True, rotate output from rtz to enz coordinate system.
    :param check_convert_pure_dp: If True, project focal mechanism onto the nearest
                                  pure double-couple before synthesizing.
    :param output_type: One of 'disp', 'strain', 'stress', 'tilt'.
    :param times_mu: If True, return raw Green's function values (already multiplied by mu).
                     If False (default), divide by mu so units are per N·m.
    :param area_km_sq: Optional subfault area in km². When provided, the result is
                       multiplied by area in m² (area_km² * 1e6), scaling the output
                       by the subfault area contribution to M0 = mu * A [m²] * slip [m].
    :param model_name: Earth model name used to look up mu (e.g. 'ak135').
    :param green_info: Pre-loaded green_lib_info dict. If None, reads green_lib_info.json
                       from path_green automatically.
    :return: 1-D numpy array of length cha_num (3 for disp, 6 for strain/stress, 2 for tilt).
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
    if area_km_sq is not None:
        v = v * area_km_sq * 1e6
    return v


def seek_edcmp2_bulk(
    path_green: str,
    event_depth_km_arr: np.ndarray,
    receiver_depth_km_arr: np.ndarray,
    az_deg_arr: np.ndarray,
    dist_km_arr: np.ndarray,
    focal_mechanism_arr: np.ndarray,
    rotate: bool = True,
    check_convert_pure_dp: bool = True,
    output_type: str = "disp",
    times_mu: bool = False,
    area_km_sq_arr: np.ndarray = None,
    slip_m_arr: np.ndarray = None,
    model_name: str = "ak135",
    green_info=None,
):
    """
    Bulk version of seek_edcmp2. Loads edcmp2_{output_type}.npy once into memory,
    then slices and synthesizes for all query points without repeated file I/O.

    Requires convert_pd2np_edcmp2_all to have been called first to generate
    edcmp2_{output_type}.npy files in path_green.

    :param path_green: Root directory of the Green's function library.
    :param event_depth_km_arr: Array of event depths in km, shape (N,).
    :param receiver_depth_km_arr: Array of receiver depths in km, shape (N,).
    :param az_deg_arr: Array of azimuths in degrees, shape (N,).
    :param dist_km_arr: Array of epicentral distances in km, shape (N,).
    :param focal_mechanism_arr: Array of focal mechanisms, shape (N, 3) for [strike,dip,rake]
                                 or (N, 6) for MT in NED.
    :param rotate: Rotate from rtz to enz.
    :param check_convert_pure_dp: Convert to pure double-couple if True.
    :param output_type: 'disp', 'strain', 'stress', or 'tilt'.
    :param times_mu: If False, divide by mu (rho*beta^2).
    :param area_km_sq_arr: Optional array of subfault areas in km², shape (N,).
                        When provided, each result row is multiplied by the
                        corresponding area in m² (area_km² * 1e6), scaling the
                        Green's function output by the subfault area contribution
                        to the seismic moment M0 = mu * A [m²] * slip [m].
    :param slip_m_arr: Optional array of subfault slips in m, shape (N,).
                       When provided, each result row is multiplied by the
                       corresponding slip, scaling the Green's function output
                       by the slip contribution to M0 = mu * A [m²] * slip [m].
    :param model_name: Earth model name for mu lookup.
    :param green_info: Pre-loaded green_lib_info dict (avoids re-reading JSON).

    :return: numpy array of shape (N, cha_num), one row per query point.
    """
    event_depth_km_arr = np.asarray(event_depth_km_arr)
    receiver_depth_km_arr = np.asarray(receiver_depth_km_arr)
    az_deg_arr = np.asarray(az_deg_arr)
    dist_km_arr = np.asarray(dist_km_arr)
    focal_mechanism_arr = np.asarray(focal_mechanism_arr)

    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)

    grn_source_depth_range = green_info["grn_source_depth_range"]
    grn_source_delta_depth = green_info["grn_source_delta_depth"]
    event_depth_arr = np.arange(
        grn_source_depth_range[0],
        grn_source_depth_range[1] + grn_source_delta_depth,
        grn_source_delta_depth,
    )
    obs_depth_list = green_info["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]
    obs_depth_arr = np.array(obs_depth_list)

    grn_dist_range = green_info["grn_dist_range"]
    grn_dist_delta = green_info["grn_delta_dist"]
    dist_arr = np.arange(
        grn_dist_range[0], grn_dist_range[1] + grn_dist_delta, grn_dist_delta
    )

    # Load bulk array once: shape (n_dep, n_obs, 5, n_dist, cha_num)
    bulk = np.load(os.path.join(path_green, "edcmp2_%s.npy" % output_type))

    N = len(event_depth_km_arr)

    # Batch nearest-neighbor index lookup for all queries at once
    dep_idx = np.argmin(
        np.abs(event_depth_arr[:, None] - event_depth_km_arr[None, :]), axis=0
    )
    obs_idx = np.argmin(
        np.abs(obs_depth_arr[:, None] - receiver_depth_km_arr[None, :]), axis=0
    )
    dist_idx = np.argmin(
        np.abs(dist_arr[:, None] - dist_km_arr[None, :]), axis=0
    )

    # Batch extract Green's function data for all queries: shape (N, 5, cha_num)
    v_raws = bulk[
        dep_idx[:, None],        # (N, 1)
        obs_idx[:, None],        # (N, 1)
        np.arange(5)[None, :],   # (1, 5)
        dist_idx[:, None],       # (N, 1)
        :,
    ]  # -> (N, 5, cha_num)

    if output_type == "disp":
        cha_num = 3
    elif output_type in ("strain", "stress"):
        cha_num = 6
    elif output_type == "tilt":
        cha_num = 2
    else:
        raise ValueError("output_type must be one of disp, strain, stress, tilt")

    # Pre-compute mu_pa for every depth in event_depth_arr (nd file read once)
    if not times_mu:
        if model_name == "ak135fc":
            from .ak135fc import s as str_nd
            lines = str_nd.split("\n")
            lines_new = []
            for line in lines:
                temp = line.split()
                if len(temp) > 1:
                    lines_new.extend(float(x) for x in temp)
            nd_model = np.array(lines_new).reshape(-1, 4)
        else:
            nd_model = read_nd(model_name)
        nd_depths = nd_model[:, 0]
        nd_mu = nd_model[:, 3] * nd_model[:, 2] ** 2 * 1e9  # rho * vs^2 (Pa)
        ind_per_dep = np.searchsorted(nd_depths, event_depth_arr, side="left")
        ind_per_dep = np.clip(ind_per_dep, 0, len(nd_depths) - 1)
        mu_per_dep = nd_mu[ind_per_dep]   # shape: (n_event_depths,)

    results = np.zeros((N, cha_num), dtype=float)
    for n in range(N):
        az_deg = float(az_deg_arr[n])
        A_rotate = create_rotate_z_mat(gamma=np.deg2rad(az_deg))
        focal_mechanism_conv = check_convert_fm(focal_mechanism_arr[n])
        mt_ned_full = tensor2full_tensor_matrix(mt=focal_mechanism_conv, flag="ned")
        mt_rotate = A_rotate.T @ mt_ned_full @ A_rotate
        mt = np.array([
            mt_rotate[0, 0], mt_rotate[0, 1], mt_rotate[0, 2],
            mt_rotate[1, 1], mt_rotate[1, 2], mt_rotate[2, 2],
        ])
        if check_convert_pure_dp:
            mt_dp = plane2mt(1, *mt2plane(mt)[0])
        else:
            mt_dp = mt

        # mt_ind=3 contributes twice: once for mt_dp[3], once extra for mt_dp[0]
        weights = mt_dp[:5].copy()
        weights[3] += mt_dp[0]

        # Synthesize: (cha_num, 5) @ (5,) = (cha_num,)
        v_ned_green_north = np.asarray(v_raws[n], dtype=float).T @ weights

        v_rtz = np.zeros(cha_num)
        if cha_num <= 3:
            # displacement / tilt: [Ux(R), Uy(T_edcmp), Uz(Up)] → t_code=-y, z_code=-z
            v_rtz[0] = v_ned_green_north[0]
            if cha_num > 1:
                v_rtz[1] = -v_ned_green_north[1]
            if cha_num > 2:
                v_rtz[2] = -v_ned_green_north[2]
            if rotate:
                v = rotate_rtz_to_enz(az_in_deg=az_deg, r=v_rtz[0], t=v_rtz[1], z=v_rtz[2])
            else:
                v = v_rtz
        else:
            # stress / strain (cha_num=6)
            # edcmp columns: [Sxx, Syy, Szz, Sxy, Syz, Szx]  (x=R, y=T_edcmp, z=Up)
            # sign convention t_code=-y, z_code=-z → Q=diag(1,-1,-1):
            #   σ_ij_code = Q_ia Q_jb σ_ab_edcmp
            # reorder to [rr, rt, rz, tt, tz, zz] = [xx, xy, xz, yy, yz, zz] for rotate fn
            v_rtz[0] =  v_ned_green_north[0]   # σ_rr = +Sxx
            v_rtz[1] = -v_ned_green_north[3]   # σ_rt = -Sxy  (one sign flip)
            v_rtz[2] = -v_ned_green_north[5]   # σ_rz = -Szx  (one sign flip)
            v_rtz[3] =  v_ned_green_north[1]   # σ_tt = +Syy  (two flips cancel)
            v_rtz[4] =  v_ned_green_north[4]   # σ_tz = +Syz  (two flips cancel)
            v_rtz[5] =  v_ned_green_north[2]   # σ_zz = +Szz  (two flips cancel)
            if rotate:
                # rotate_symmetric_tensor_series with γ = az - π/2 gives:
                #   R.T @ σ_rtz @ R = A_enz @ σ_rtz @ A_enz.T
                # where A_enz = [[sin(az),-cos(az),0],[cos(az),sin(az),0],[0,0,1]]
                # output order: [ee, en, ez, nn, nz, zz]
                gamma = np.deg2rad(az_deg) - np.pi / 2
                v = rotate_symmetric_tensor_series(v_rtz.reshape(1, 6), gamma)[0]
            else:
                v = v_rtz

        if not times_mu:
            v = v / mu_per_dep[dep_idx[n]]

        if area_km_sq_arr is not None:
            v = v * area_km_sq_arr[n] * 1e6  # km² -> m²

        if slip_m_arr is not None:
            v = v * slip_m_arr[n]

        results[n] = v

    return results
