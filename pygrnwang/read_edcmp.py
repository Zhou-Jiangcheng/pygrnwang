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

_EDCMP_CHA_NUM = {"disp": 3, "strain": 6, "stress": 6, "tilt": 2}


def _get_edcmp_cha_num(output_type):
    try:
        return _EDCMP_CHA_NUM[output_type]
    except KeyError:
        raise ValueError("output_type must be one of disp, strain, stress, tilt")


def _nearest_grid_indices(values, grid):
    values = np.asarray(values)
    grid = np.asarray(grid)
    if grid.ndim != 1 or len(grid) == 0:
        raise ValueError("grid must be a non-empty 1-D array")
    if len(grid) == 1:
        return np.zeros(values.shape, dtype=np.intp)

    grid_order = np.argsort(grid, kind="stable")
    grid_sorted = grid[grid_order]
    right_idx = np.searchsorted(grid_sorted, values, side="left")
    right_idx = np.clip(right_idx, 1, len(grid_sorted) - 1)
    left_idx = right_idx - 1

    left_dist = np.abs(values - grid_sorted[left_idx])
    right_dist = np.abs(grid_sorted[right_idx] - values)
    left_orig_idx = grid_order[left_idx]
    right_orig_idx = grid_order[right_idx]
    choose_left = (left_dist < right_dist) | (
        (left_dist == right_dist) & (left_orig_idx <= right_orig_idx)
    )
    sorted_indices = np.where(choose_left, left_idx, right_idx)
    return grid_order[sorted_indices].astype(np.intp)


def _nearest_regular_grid_indices(values, grid_start, grid_delta, grid_size):
    values = np.asarray(values)
    indices = np.ceil((values - grid_start) / grid_delta - 0.5).astype(np.intp)
    return np.clip(indices, 0, grid_size - 1)


def read_edcmp_raw(path_green, output_type, grn_event_depth, grn_obs_depth, mt_ind):
    cha_num = _get_edcmp_cha_num(output_type)
    path_bin = str(
        os.path.join(
            path_green,
            "edcmp2",
            "%.2f" % grn_event_depth,
            "%.2f" % grn_obs_depth,
            "%d" % mt_ind,
            "%s.bin" % output_type,
        )
    )
    if os.path.exists(path_bin):
        values_raw = np.fromfile(path_bin, dtype=np.float32).reshape(-1, cha_num)
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
        values_raw = np.asarray(df.to_numpy()[:, 2:], dtype=np.float32)
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

    cha_num = _get_edcmp_cha_num(output_type)

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

    if cha_num == 3:
        v_rtz = np.zeros_like(v_ned_green_north)
        v_rtz[0] = v_ned_green_north[0]
        v_rtz[1] = -v_ned_green_north[1]
        v_rtz[2] = -v_ned_green_north[2]
        if rotate:
            v = rotate_rtz_to_enz(az_in_deg=az_deg, r=v_rtz[0], t=v_rtz[1], z=v_rtz[2])
        else:
            v = v_rtz
    elif cha_num == 2:
        v_rt = np.zeros_like(v_ned_green_north)
        v_rt[0] = v_ned_green_north[0]
        v_rt[1] = -v_ned_green_north[1]
        if rotate:
            v = rotate_rtz_to_enz(az_in_deg=az_deg, r=v_rt[0], t=v_rt[1], z=0.0)[:2]
        else:
            v = v_rt
    else:
        v_rtz = np.zeros_like(v_ned_green_north)
        v_rtz[0] = v_ned_green_north[0]
        v_rtz[1] = -v_ned_green_north[3]
        v_rtz[2] = -v_ned_green_north[5]
        v_rtz[3] = v_ned_green_north[1]
        v_rtz[4] = v_ned_green_north[4]
        v_rtz[5] = v_ned_green_north[2]
        if rotate:
            gamma = np.deg2rad(az_deg) - np.pi / 2
            v = rotate_symmetric_tensor_series(v_rtz.reshape(1, 6), gamma)[0]
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
    Bulk version of seek_edcmp2. Synthesizes all N query points in a single
    vectorized pass (no Python loop).

    Requires convert_pd2bin_edcmp2_all to have been called first to generate
    edcmp2_{output_type}.bin files in path_green.

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
    event_depth_km_arr = np.asarray(event_depth_km_arr, dtype=float)
    receiver_depth_km_arr = np.asarray(receiver_depth_km_arr, dtype=float)
    az_deg_arr = np.asarray(az_deg_arr, dtype=float)
    dist_km_arr = np.asarray(dist_km_arr, dtype=float)
    focal_mechanism_arr = np.asarray(focal_mechanism_arr, dtype=float)

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

    cha_num = _get_edcmp_cha_num(output_type)

    # Load bulk array once: shape (n_dep, n_obs, 5, n_dist, cha_num)
    bulk_shape = (len(event_depth_arr), len(obs_depth_arr), 5, len(dist_arr), cha_num)
    bulk = np.fromfile(
        os.path.join(path_green, "edcmp2_%s.bin" % output_type), dtype=np.float32
    ).reshape(bulk_shape)

    N = len(event_depth_km_arr)

    dep_idx = _nearest_regular_grid_indices(
        event_depth_km_arr,
        event_depth_arr[0],
        grn_source_delta_depth,
        len(event_depth_arr),
    )
    obs_idx = _nearest_grid_indices(receiver_depth_km_arr, obs_depth_arr)
    dist_idx = _nearest_regular_grid_indices(
        dist_km_arr,
        dist_arr[0],
        grn_dist_delta,
        len(dist_arr),
    )

    # Batch extract Green's function data for all queries: shape (N, 5, cha_num)
    v_raws = np.array(
        bulk[
            dep_idx[:, None],  # (N, 1)
            obs_idx[:, None],  # (N, 1)
            np.arange(5)[None, :],  # (1, 5)
            dist_idx[:, None],  # (N, 1)
            :,
        ],
        dtype=float,
    )  # (N, 5, cha_num)

    # Pre-compute mu for all grid depths (nd file read once)
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
        mu_per_dep = nd_mu[ind_per_dep]  # shape: (n_event_depths,)

    # --- vectorized focal mechanism → per-point weights (N, 5) ---

    fm = focal_mechanism_arr
    if fm.ndim == 1:
        fm = fm.reshape(1, -1)

    if fm.shape[1] == 6:
        mt_ned = fm  # already [M11, M12, M13, M22, M23, M33] in NED
    else:
        mt_ned = np.array([check_convert_fm(fm[n]) for n in range(N)], dtype=float)

    # Rotate MT by azimuth: mt_r[n] = A(az[n]).T @ MT[n] @ A(az[n])
    # For a z-rotation matrix A with c=cos(az), s=sin(az), the result is:
    #   M11' = c²M11 + 2cs·M12 + s²M22
    #   M12' = -cs·M11 + (c²-s²)·M12 + cs·M22
    #   M13' = c·M13 + s·M23
    #   M22' = s²M11 - 2cs·M12 + c²M22
    #   M23' = -s·M13 + c·M23
    #   M33' = M33
    gamma_arr = np.deg2rad(az_deg_arr)
    c, s = np.cos(gamma_arr), np.sin(gamma_arr)
    c2, s2, cs = c * c, s * s, c * s
    M11, M12, M13 = mt_ned[:, 0], mt_ned[:, 1], mt_ned[:, 2]
    M22, M23, M33 = mt_ned[:, 3], mt_ned[:, 4], mt_ned[:, 5]
    mt_r = np.stack(
        [
            c2 * M11 + 2 * cs * M12 + s2 * M22,
            -cs * M11 + (c2 - s2) * M12 + cs * M22,
            c * M13 + s * M23,
            s2 * M11 - 2 * cs * M12 + c2 * M22,
            -s * M13 + c * M23,
            M33,
        ],
        axis=1,
    )  # (N, 6)

    if check_convert_pure_dp:
        for n in range(N):
            mt_r[n] = plane2mt(1, *mt2plane(mt_r[n])[0])

    # mt_ind=3 contributes twice: once for mt_r[3], once extra for mt_r[0]
    weights = mt_r[:, :5].copy()  # (N, 5)
    weights[:, 3] += mt_r[:, 0]

    # Synthesize all N points at once: v_ned[n] = v_raws[n].T @ weights[n]
    # v_raws: (N, 5, cha_num), weights: (N, 5) → v_ned: (N, cha_num)
    v_ned = np.einsum("nkc,nk->nc", v_raws, weights)

    # --- coordinate transforms (fully vectorised) ---

    if cha_num == 3:
        # displacement: [Ux(R), Uy(T_edcmp), Uz(Up)] → t_code=-y, z_code=-z
        v_rtz = np.empty((N, 3))
        v_rtz[:, 0] = v_ned[:, 0]
        v_rtz[:, 1] = -v_ned[:, 1]
        v_rtz[:, 2] = -v_ned[:, 2]
        if rotate:
            az = gamma_arr
            results = np.stack(
                [
                    v_rtz[:, 0] * np.sin(az) - v_rtz[:, 1] * np.cos(az),
                    v_rtz[:, 0] * np.cos(az) + v_rtz[:, 1] * np.sin(az),
                    v_rtz[:, 2],
                ],
                axis=1,
            )
        else:
            results = v_rtz

    elif cha_num == 2:
        # tilt: horizontal components only
        v_rt = np.empty((N, 2))
        v_rt[:, 0] = v_ned[:, 0]
        v_rt[:, 1] = -v_ned[:, 1]
        if rotate:
            az = gamma_arr
            results = np.stack(
                [
                    v_rt[:, 0] * np.sin(az) - v_rt[:, 1] * np.cos(az),
                    v_rt[:, 0] * np.cos(az) + v_rt[:, 1] * np.sin(az),
                ],
                axis=1,
            )
        else:
            results = v_rt

    else:
        # stress / strain (cha_num=6)
        # edcmp columns: [Sxx, Syy, Szz, Sxy, Syz, Szx]  (x=R, y=T_edcmp, z=Up)
        # sign convention t_code=-y, z_code=-z → Q=diag(1,-1,-1):
        #   σ_ij_code = Q_ia Q_jb σ_ab_edcmp
        # reorder to [rr, rt, rz, tt, tz, zz] = [xx, xy, xz, yy, yz, zz] for rotate fn
        v_rtz = np.empty((N, 6))
        v_rtz[:, 0] = v_ned[:, 0]  # σ_rr = +Sxx
        v_rtz[:, 1] = -v_ned[:, 3]  # σ_rt = -Sxy  (one sign flip)
        v_rtz[:, 2] = -v_ned[:, 5]  # σ_rz = -Szx  (one sign flip)
        v_rtz[:, 3] = v_ned[:, 1]  # σ_tt = +Syy  (two flips cancel)
        v_rtz[:, 4] = v_ned[:, 4]  # σ_tz = +Syz  (two flips cancel)
        v_rtz[:, 5] = v_ned[:, 2]  # σ_zz = +Szz  (two flips cancel)
        if rotate:
            # Analytical z-rotation by gamma2 = az - π/2 for each point.
            # For symmetric tensor [xx,xy,xz,yy,yz,zz] and R=R_z(gamma2):
            #   xx' = c²xx + 2cs·xy + s²yy
            #   xy' = -cs·xx + (c²-s²)·xy + cs·yy
            #   xz' = c·xz + s·yz
            #   yy' = s²xx - 2cs·xy + c²yy
            #   yz' = -s·xz + c·yz
            #   zz' = zz
            # output order matches rotate_symmetric_tensor_series: [ee,en,ez,nn,nz,zz]
            gamma2 = gamma_arr - np.pi / 2
            cr, sr = np.cos(gamma2), np.sin(gamma2)
            cr2, sr2, csr = cr * cr, sr * sr, cr * sr
            xx, xy, xz = v_rtz[:, 0], v_rtz[:, 1], v_rtz[:, 2]
            yy, yz, zz = v_rtz[:, 3], v_rtz[:, 4], v_rtz[:, 5]
            results = np.stack(
                [
                    cr2 * xx + 2 * csr * xy + sr2 * yy,
                    -csr * xx + (cr2 - sr2) * xy + csr * yy,
                    cr * xz + sr * yz,
                    sr2 * xx - 2 * csr * xy + cr2 * yy,
                    -sr * xz + cr * yz,
                    zz,
                ],
                axis=1,
            )
        else:
            results = v_rtz

    # --- scaling ---

    if not times_mu:
        results = results / mu_per_dep[dep_idx][:, np.newaxis]

    if area_km_sq_arr is not None:
        results = results * (area_km_sq_arr[:, np.newaxis] * 1e6)  # km² -> m²

    if slip_m_arr is not None:
        results = results * slip_m_arr[:, np.newaxis]

    return results.astype(np.float32)
