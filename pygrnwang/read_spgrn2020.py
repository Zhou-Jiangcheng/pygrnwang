import os
import json
from collections import OrderedDict
from typing import Union

import numpy as np
from scipy import signal

from .read_tpts_table import read_tpts_table
from .utils import shift_green2real_tpts
from .focal_mechanism import check_convert_fm
from .geo import rotate_rtz_to_enz
from .signal_process import resample, filter_butter


def synthesize_spgrn(az_in_deg, time_series, focal_mechanism):
    [M11, M12, M13, M22, M23, M33] = check_convert_fm(focal_mechanism=focal_mechanism)
    # expl.,strike-slip,dip-slip,clvd
    exp = (M11 + M22 + M33) / 3
    ss1 = -M12
    ss2 = (M11 - M22) / 2
    ds1 = M13
    ds2 = -M23
    clvd = M33 - exp

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
    r = (
        time_series[1] * m1[0]
        + time_series[3] * m1[1]
        + time_series[6] * m1[2]
        + time_series[9] * m1[3]
    )
    t = time_series[4] * m2[0] + time_series[7] * m2[1]

    seismograms = np.array([r, t, z])
    return seismograms


def read_spgrn_data_by_index(path_grn_data, dist_index, green_info):
    """Helper to read raw data at a specific distance index."""
    N_T = round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    # Structure based on original code: 3 header floats + (2 + N_T) floats * 10 components
    length_each = 3 + (2 + N_T) * 10
    start_count = dist_index * length_each

    with open(path_grn_data, "rb") as fr:
        time_series = np.fromfile(
            file=fr, dtype=np.float32, count=length_each, offset=start_count * 4
        )

    # Reshape logic from original code
    time_series = time_series[3:].reshape(10, 2 + N_T)
    time_series = time_series[:, 1:-1]
    return time_series


def read_spgrn_data_two_indices(path_grn_data, dist_idx_low, dist_idx_high, green_info):
    """Read two distance indices in a single file open for efficiency."""
    N_T = round(green_info["time_window"] / green_info["sampling_interval"]) + 1
    length_each = 3 + (2 + N_T) * 10
    bytes_each = length_each * 4

    with open(path_grn_data, "rb") as fr:
        fr.seek(dist_idx_low * bytes_each)
        raw_low = np.frombuffer(fr.read(bytes_each), dtype=np.float32).copy()
        fr.seek(dist_idx_high * bytes_each)
        raw_high = np.frombuffer(fr.read(bytes_each), dtype=np.float32).copy()

    shape = (10, 2 + N_T)
    data_low = raw_low[3:].reshape(shape)[:, 1:-1]
    data_high = raw_high[3:].reshape(shape)[:, 1:-1]
    return data_low, data_high


def get_sorted_grid_params(target, grid_list):
    """
    Helper to find neighbors and weight for 1D interpolation on a sorted list.
    Returns: val_low, val_high, weight_high, idx_low
    """
    arr = np.array(grid_list)
    if len(arr) == 0:
        raise ValueError("Grid list is empty")
    if len(arr) == 1:
        return arr[0], arr[0], 0.0, 0

    # Handle out of bounds by clamping
    if target <= arr[0]:
        return arr[0], arr[0], 0.0, 0
    if target >= arr[-1]:
        return arr[-1], arr[-1], 0.0, len(arr) - 1

    idx = np.searchsorted(arr, target)
    # arr[idx-1] <= target <= arr[idx]
    val_low = arr[idx - 1]
    val_high = arr[idx]

    if val_high == val_low:
        weight = 0.0
    else:
        weight = (target - val_low) / (val_high - val_low)

    return val_low, val_high, weight, idx - 1


def seek_spgrn2020(
    path_green: str,
    event_depth_km: float,
    receiver_depth_km: float,
    az_deg: float,
    dist_km: float,
    focal_mechanism: Union[np.ndarray, list],
    srate: float,
    output_type: str = "disp",
    rotate: bool = True,
    before_p: Union[float, None] = None,
    pad_zeros: bool = False,
    shift: bool = False,
    only_seismograms: bool = True,
    model_name: str = "ak135fc",
    green_info: Union[dict, None] = None,
    interpolate_type: int = 0,
    freq_band=None,
    butter_order: int = 4,
    zero_phase: bool = False,
):
    """
    Read synthetic seismograms.

    :param path_green: Root directory of the data.
    :param event_depth_km: Event depth in km.
    :param receiver_depth_km: Receiver depth in km.
    :param az_deg: Azimuth in degrees.
    :param dist_km: Epicentral distance in km.
    :param focal_mechanism: [strike, dip, rake] or [M11, M12, M13, M22, M23, M33].
    :param output_type: 'disp', 'velo', 'acce'.
    :param srate: Sampling rate in Hz.
    :param before_p: Time before P-wave.
    :param pad_zeros: Pad with zeros.
    :param shift: Shift seismograms based on tpts.
    :param rotate: Rotate rtz2ned.
    :param only_seismograms: Return only seismograms.
    :param model_name: Model name.
    :param green_info: Green's function library info.
    :param interpolate_type:
            0 for nearest neighbor,
            1 for trilinear interpolation (Source Depth, Receiver Depth, Distance).
    :param freq_band: Frequency band for bandpass filter [low_freq, high_freq] in Hz.
            Use None or [None, None] for no filtering (default).
            Use [low_freq, None] for highpass, [None, high_freq] for lowpass.
    :param butter_order: Order of Butterworth filter (default: 4).
    :param zero_phase: Whether to use zero-phase filtering (default: False).
    :return: (
            seismograms_resample,
            tpts_table,
            first_p,
            first_s,
            grn_dep_source,
            grn_dep_receiver,
            grn_dist,
        )
    """
    if green_info is None:
        with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
            green_info = json.load(fr)
    srate_grn = 1 / green_info["sampling_interval"]
    sampling_num = green_info["samples_num"]
    green_before_p = green_info["green_before_p"]
    grn_dep_list = green_info["event_depth_list"]
    grn_receiver_list = green_info["receiver_depth_list"]
    dist_list = green_info["dist_list"]

    # --- 1. Identify Nearest Neighbors (Used for Metadata) ---
    # These are still needed for reading tpts tables which typically don't support interpolation well
    if not isinstance(grn_dep_list, list):
        grn_dep_source = grn_dep_list
    else:
        grn_dep_source = grn_dep_list[
            np.argmin(np.abs(event_depth_km - np.array(grn_dep_list)))
        ]
    if not isinstance(grn_receiver_list, list):
        grn_dep_receiver = grn_receiver_list
    else:
        grn_dep_receiver = grn_receiver_list[
            np.argmin(np.abs(receiver_depth_km - np.array(grn_receiver_list)))
        ]

    # Nearest distance logic for metadata
    nearest_dist_idx = np.argmin(np.abs(np.array(dist_list) - dist_km))
    grn_dist = dist_list[nearest_dist_idx]

    # --- 2. Retrieve Data based on Interpolation Type ---

    if interpolate_type == 0:
        # === Type 0: Nearest Neighbor ===
        path_greenfunc = str(
            os.path.join(
                path_green,
                "GreenFunc",
                "%.2f" % grn_dep_source,
                "%.2f" % grn_dep_receiver,
            )
        )
        path_grn_data = os.path.join(path_greenfunc, "grn_d%.2f" % grn_dep_source)
        time_series = read_spgrn_data_by_index(
            path_grn_data, nearest_dist_idx, green_info
        )

    else:
        # === Type 1: Trilinear Interpolation (Source Depth & Receiver Depth & Distance) ===

        # A. Source Depth Interpolation Parameters
        if not isinstance(grn_dep_list, list):
            d_src_low, d_src_high, w_src, _ = grn_dep_list, grn_dep_list, 0.0, 0
        else:
            d_src_low, d_src_high, w_src, _ = get_sorted_grid_params(
                event_depth_km, grn_dep_list
            )

        # B. Receiver Depth Interpolation Parameters (NEW)
        if not isinstance(grn_receiver_list, list):
            d_rec_low, d_rec_high, w_rec, _ = (
                grn_receiver_list,
                grn_receiver_list,
                0.0,
                0,
            )
        else:
            d_rec_low, d_rec_high, w_rec, _ = get_sorted_grid_params(
                receiver_depth_km, grn_receiver_list
            )

        # C. Distance Interpolation Parameters
        _, _, w_dist, dist_idx_low = get_sorted_grid_params(dist_km, dist_list)
        dist_idx_high = min(dist_idx_low + 1, len(dist_list) - 1)

        if dist_idx_low == dist_idx_high:
            w_dist = 0.0

        # D. Helper to fetch and interpolate distance for a specific (Source, Receiver) pair
        def fetch_distance_layer(src_depth, rec_depth):
            path_gf_layer = str(
                os.path.join(
                    path_green,
                    "GreenFunc",
                    "%.2f" % src_depth,
                    "%.2f" % rec_depth,
                )
            )
            path_data_layer = os.path.join(path_gf_layer, "grn_d%.2f" % src_depth)

            if w_dist > 1e-4:
                # Read both indices in one file open
                data_d0, data_d1 = read_spgrn_data_two_indices(
                    path_data_layer, dist_idx_low, dist_idx_high, green_info
                )
                return (1 - w_dist) * data_d0 + w_dist * data_d1
            else:
                return read_spgrn_data_by_index(
                    path_data_layer, dist_idx_low, green_info
                )

        # E. Helper to fetch and interpolate Receiver Depth for a specific Source Depth
        def fetch_receiver_layer(src_depth):
            # 1. Fetch Low Receiver Depth
            data_r0 = fetch_distance_layer(src_depth, d_rec_low)

            # 2. Fetch High Receiver Depth (if needed)
            if w_rec > 1e-4 and d_rec_high != d_rec_low:
                data_r1 = fetch_distance_layer(src_depth, d_rec_high)
                return (1 - w_rec) * data_r0 + w_rec * data_r1
            else:
                return data_r0

        # F. Perform Source Depth Interpolation (Top Level)
        # 1. Low Source Depth
        ts_low_src = fetch_receiver_layer(d_src_low)

        # 2. High Source Depth (if needed)
        if w_src > 1e-4 and d_src_high != d_src_low:
            ts_high_src = fetch_receiver_layer(d_src_high)
            time_series = (1 - w_src) * ts_low_src + w_src * ts_high_src
        else:
            time_series = ts_low_src

    # --- 3. Synthesize Seismograms ---
    # r,t,z
    seismograms = synthesize_spgrn(
        az_in_deg=az_deg, time_series=time_series, focal_mechanism=focal_mechanism
    )
    if rotate:
        seismograms = rotate_rtz_to_enz(
            az_in_deg=az_deg, r=seismograms[0], t=seismograms[1], z=seismograms[2]
        )[:]

    # --- 4. Post-processing (Time shifting, resampling) ---
    # Get metadata from the nearest neighbor path (Metadata uses Nearest Neighbor)
    path_greenfunc_meta = str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % grn_dep_source,
            "%.2f" % grn_dep_receiver,
        )
    )
    tpts_table = read_tpts_table(
        path_greenfunc=path_greenfunc_meta, dist_in_km=dist_km, green_info=green_info
    )

    # Apply bandpass filter (vectorized over all 3 components at once)
    if freq_band is not None and (freq_band[0] is not None or freq_band[1] is not None):
        seismograms = filter_butter(
            seismograms, srate_grn, freq_band, butter_order, zero_phase
        )

    ts_count = 0
    if before_p is not None:
        ts_count = round((green_before_p - before_p) * srate_grn)
    if pad_zeros:
        if before_p is not None:
            raise ValueError("can not set before_p and pad_zeros together")
        ts_count = round((green_before_p - tpts_table["p_onset"]) * srate_grn)
        before_p = tpts_table["p_onset"]
    seismograms = np.roll(seismograms, -ts_count)
    if ts_count > 0:
        seismograms[:, -ts_count:] = 0
    elif ts_count < 0:
        seismograms[:, :-ts_count] = 0

    first_p = None
    first_s = None
    if shift:
        seismograms, first_p, first_s = shift_green2real_tpts(
            seismograms=seismograms,
            tpts_table=tpts_table,
            srate=srate_grn,
            green_before_p=before_p,
            event_depth_km=event_depth_km,
            dist_in_km=dist_km,
            model_name=model_name,
        )

    len_after_resample = round(sampling_num * srate / srate_grn)
    # Vectorized resample: resample_poly supports 2D arrays via axis parameter
    if float(srate_grn).is_integer() and float(srate).is_integer():
        gcd = np.gcd(int(srate), int(srate_grn))
        p = int(srate) // gcd
        q = int(srate_grn) // gcd
        seismograms_resample = signal.resample_poly(seismograms, p, q, axis=1)[
            :, :len_after_resample
        ]
    else:
        seismograms_resample = np.zeros((3, len_after_resample))
        for i in range(3):
            seismograms_resample[i] = resample(
                seismograms[i], srate_old=srate_grn, srate_new=srate, zero_phase=True
            )[:len_after_resample]

    if output_type == "disp":
        seismograms_resample = np.cumsum(seismograms_resample, axis=1) / srate
    elif output_type == "acce":
        seismograms_resample = (
            signal.convolve(
                seismograms_resample.T,
                np.array([1, -1])[:, None],
                mode="same",
                method="auto",
            ).T
            / srate
        )

    if only_seismograms:
        return seismograms_resample
    else:
        return (
            seismograms_resample,
            tpts_table,
            first_p,
            first_s,
            grn_dep_source,
            grn_dep_receiver,
            grn_dist,
        )


# =============================================================================
# Cached / deduplicated reading backend for dynamic FFI
# =============================================================================
# seek_spgrn2020 above re-reads the (<=4) interpolation-corner library blocks and
# re-runs the whole filter/resample/integration chain for every (station,
# subfault, component). For a finite-fault inversion that is hundreds of
# thousands of redundant reads/decodes. The cache below reads each distinct
# (depth_node, receiver_node, dist_index) block once; synthesize_from_cache then
# (a) bilinearly interpolates from the cache with the *same* weights as
# seek_spgrn2020 (bit-identical time_series), (b) synthesizes a *list* of focal
# mechanisms in one pass (dedup of the per-component duplication), and (c)
# optionally truncates the processing window to the only samples that survive
# (exact for causal filtering + integer-ratio resampling). seek_spgrn2020 is left
# untouched so other callers are unaffected.


class GridGFCache:
    """In-process cache of raw spgrn2020 library grid blocks and tpts tables.

    A *block* is the decoded ``(10, N_T)`` array at one
    ``(event_depth_node, receiver_depth_node, dist_index)`` library node -- the
    exact array :func:`read_spgrn_data_by_index` returns. Reading each distinct
    block once removes the cross-pair redundancy (~10x for teleseismic FFI).

    If ``max_blocks`` is set, the cache is an LRU over individual blocks: once
    more than ``max_blocks`` are held the least-recently-used block is evicted.
    Each block is ``10 * N_T`` float32 (~0.3 MB for N_T=8192), so this bounds the
    cache memory precisely (low-memory mode). ``None`` keeps every distinct block
    (fastest; ~5 GB for a full teleseismic FFI). The tpts tables (~24 KB each,
    one per depth folder) are always kept -- they are negligible.
    """

    def __init__(self, path_green, green_info, max_blocks=None):
        self.path_green = path_green
        self.green_info = green_info
        self.grn_dep_list = green_info["event_depth_list"]
        self.grn_receiver_list = green_info["receiver_depth_list"]
        self.dist_list = green_info["dist_list"]
        self.dist_arr = np.asarray(self.dist_list, dtype=float)
        self.max_blocks = max_blocks
        # _blocks[(depth_node, rec_node, dist_idx)] -> (10, N_T) float32 (LRU order)
        self._blocks = OrderedDict()
        # _tpts[(depth_node, rec_node)] -> ((n_dist, 3) P, (n_dist, 3) S)
        self._tpts = {}

    def get_block(self, depth_node, rec_node, dist_idx):
        key = (depth_node, rec_node, dist_idx)
        b = self._blocks.get(key)
        if b is not None:
            self._blocks.move_to_end(key)
            return b
        path_grn_data = os.path.join(
            self.path_green,
            "GreenFunc",
            "%.2f" % depth_node,
            "%.2f" % rec_node,
            "grn_d%.2f" % depth_node,
        )
        b = read_spgrn_data_by_index(path_grn_data, dist_idx, self.green_info)
        self._blocks[key] = b
        if self.max_blocks is not None and len(self._blocks) > self.max_blocks:
            self._blocks.popitem(last=False)
        return b

    def _get_tpts_tables(self, depth_node, rec_node):
        key = (depth_node, rec_node)
        tab = self._tpts.get(key)
        if tab is None:
            folder = os.path.join(
                self.path_green, "GreenFunc", "%.2f" % depth_node, "%.2f" % rec_node
            )
            n = len(self.dist_list)
            # tptable.dat / tstable.dat: 4-byte Fortran record marker, then n
            # triples (onset, takeoff, slowness) of float32, then trailer.
            p_tab = np.fromfile(
                os.path.join(folder, "tptable.dat"),
                dtype=np.float32,
                count=n * 3,
                offset=4,
            ).reshape(n, 3)
            s_tab = np.fromfile(
                os.path.join(folder, "tstable.dat"),
                dtype=np.float32,
                count=n * 3,
                offset=4,
            ).reshape(n, 3)
            tab = (p_tab, s_tab)
            self._tpts[key] = tab
        return tab

    def _interp_params(self, event_depth_km, receiver_depth_km, dist_km):
        if not isinstance(self.grn_dep_list, list):
            d_src_low, d_src_high, w_src = self.grn_dep_list, self.grn_dep_list, 0.0
        else:
            d_src_low, d_src_high, w_src, _ = get_sorted_grid_params(
                event_depth_km, self.grn_dep_list
            )
        if not isinstance(self.grn_receiver_list, list):
            d_rec_low, d_rec_high, w_rec = (
                self.grn_receiver_list,
                self.grn_receiver_list,
                0.0,
            )
        else:
            d_rec_low, d_rec_high, w_rec, _ = get_sorted_grid_params(
                receiver_depth_km, self.grn_receiver_list
            )
        _, _, w_dist, dist_idx_low = get_sorted_grid_params(dist_km, self.dist_list)
        dist_idx_high = min(dist_idx_low + 1, len(self.dist_list) - 1)
        if dist_idx_low == dist_idx_high:
            w_dist = 0.0
        return (
            d_src_low,
            d_src_high,
            w_src,
            d_rec_low,
            d_rec_high,
            w_rec,
            w_dist,
            dist_idx_low,
            dist_idx_high,
        )

    def time_series(self, event_depth_km, receiver_depth_km, dist_km, n_cols=None):
        """Interpolated ``(10, N)`` time series, identical to seek_spgrn2020.

        When ``n_cols`` is given, each cached block is sliced to its first
        ``n_cols`` columns before combining -- the first ``n_cols`` columns of the
        full interpolation, used by the truncated processing window.
        """
        (
            d_src_low,
            d_src_high,
            w_src,
            d_rec_low,
            d_rec_high,
            w_rec,
            w_dist,
            dist_idx_low,
            dist_idx_high,
        ) = self._interp_params(event_depth_km, receiver_depth_km, dist_km)

        def fetch_distance_layer(src_depth, rec_depth):
            if w_dist > 1e-4:
                d0 = self.get_block(src_depth, rec_depth, dist_idx_low)
                d1 = self.get_block(src_depth, rec_depth, dist_idx_high)
                if n_cols is not None:
                    d0 = d0[:, :n_cols]
                    d1 = d1[:, :n_cols]
                return (1 - w_dist) * d0 + w_dist * d1
            d0 = self.get_block(src_depth, rec_depth, dist_idx_low)
            if n_cols is not None:
                d0 = d0[:, :n_cols]
            return d0

        def fetch_receiver_layer(src_depth):
            data_r0 = fetch_distance_layer(src_depth, d_rec_low)
            if w_rec > 1e-4 and d_rec_high != d_rec_low:
                data_r1 = fetch_distance_layer(src_depth, d_rec_high)
                return (1 - w_rec) * data_r0 + w_rec * data_r1
            return data_r0

        ts_low = fetch_receiver_layer(d_src_low)
        if w_src > 1e-4 and d_src_high != d_src_low:
            ts_high = fetch_receiver_layer(d_src_high)
            return (1 - w_src) * ts_low + w_src * ts_high
        return ts_low

    def tpts_table(self, event_depth_km, receiver_depth_km, dist_km):
        """Nearest-neighbour tpts dict, identical to read_tpts_table."""
        if not isinstance(self.grn_dep_list, list):
            dep_node = self.grn_dep_list
        else:
            arr = np.asarray(self.grn_dep_list, dtype=float)
            dep_node = self.grn_dep_list[int(np.argmin(np.abs(event_depth_km - arr)))]
        if not isinstance(self.grn_receiver_list, list):
            rec_node = self.grn_receiver_list
        else:
            arr = np.asarray(self.grn_receiver_list, dtype=float)
            rec_node = self.grn_receiver_list[
                int(np.argmin(np.abs(receiver_depth_km - arr)))
            ]
        p_tab, s_tab = self._get_tpts_tables(dep_node, rec_node)
        idx = int(np.argmin(np.abs(self.dist_arr - dist_km)))
        p = p_tab[idx]
        s = s_tab[idx]
        return {
            "p_onset": float(p[0]),
            "p_takeoff": float(p[1]),
            "p_slowness": float(p[2]),
            "s_onset": float(s[0]),
            "s_takeoff": float(s[1]),
            "s_slowness": float(s[2]),
        }


def synthesize_from_cache(
    cache: GridGFCache,
    event_depth_km: float,
    receiver_depth_km: float,
    az_deg: float,
    dist_km: float,
    focal_mechanisms,
    srate: float,
    output_type: str = "disp",
    rotate: bool = True,
    before_p: Union[float, None] = None,
    pad_zeros: bool = False,
    shift: bool = False,
    freq_band=None,
    butter_order: int = 4,
    zero_phase: bool = False,
    n_keep: Union[int, None] = None,
):
    """Synthesize seismograms for a *list* of focal mechanisms from a cache.

    Equivalent to calling :func:`seek_spgrn2020` once per focal mechanism with
    ``shift=False`` and slicing the result to ``n_keep`` samples, but it reads
    each library block at most once (via ``cache``), runs the post-processing
    chain a single time over the stacked mechanisms, and -- when safe -- only
    processes the ``n_keep`` samples that survive.

    :return: ``(seis_stack, tpts_table)`` where ``seis_stack`` has shape
             ``(n_fm, 3, n_out)`` (ENZ if ``rotate``), ``n_out == n_keep`` when
             ``n_keep`` is given.
    """
    green_info = cache.green_info
    srate_grn = 1 / green_info["sampling_interval"]
    sampling_num = green_info["samples_num"]
    green_before_p = green_info["green_before_p"]

    tpts_table = cache.tpts_table(event_depth_km, receiver_depth_km, dist_km)

    # Shift amount (mirror seek_spgrn2020).
    ts_count = 0
    if before_p is not None:
        ts_count = round((green_before_p - before_p) * srate_grn)
    if pad_zeros:
        if before_p is not None:
            raise ValueError("can not set before_p and pad_zeros together")
        ts_count = round((green_before_p - tpts_table["p_onset"]) * srate_grn)

    # Truncation is exact only with causal filtering + integer-ratio resampling
    # (the FFT resample path has global support). Disabled otherwise.
    int_srate = float(srate_grn).is_integer() and float(srate).is_integer()
    truncate = (
        n_keep is not None
        and not shift
        and not zero_phase
        and int_srate
        and ts_count >= 0
        and before_p is not None
        and not pad_zeros
    )
    n_cols = None
    if truncate:
        gcd = np.gcd(int(srate), int(srate_grn))
        p = int(srate) // gcd
        q = int(srate_grn) // gcd
        n_keep_grn = int(np.ceil(n_keep * q / p))
        guard = 10 * max(p, q) + 64
        n_proc = ts_count + n_keep_grn + guard
        if n_proc >= sampling_num:
            truncate = False
        else:
            n_cols = n_proc

    time_series = cache.time_series(
        event_depth_km, receiver_depth_km, dist_km, n_cols=n_cols
    )

    seis_list = []
    for fm in focal_mechanisms:
        s = synthesize_spgrn(
            az_in_deg=az_deg, time_series=time_series, focal_mechanism=fm
        )
        if rotate:
            s = rotate_rtz_to_enz(az_in_deg=az_deg, r=s[0], t=s[1], z=s[2])[:]
        seis_list.append(s)
    seis = np.stack(seis_list, axis=0)  # (n_fm, 3, N)

    if freq_band is not None and (freq_band[0] is not None or freq_band[1] is not None):
        seis = filter_butter(seis, srate_grn, freq_band, butter_order, zero_phase)

    # Time-align by ts_count. seek_spgrn2020 does np.roll(-ts_count) over the
    # flattened (3, N) array then zeros the wrapped tail/head; the net effect is a
    # per-trace shift with zero padding, reproduced here on the (n_fm, 3, N) stack.
    N = seis.shape[-1]
    if ts_count > 0:
        shifted = np.zeros_like(seis)
        shifted[..., : N - ts_count] = seis[..., ts_count:]
        seis = shifted
    elif ts_count < 0:
        a = -ts_count
        shifted = np.zeros_like(seis)
        shifted[..., a:] = seis[..., : N - a]
        seis = shifted

    len_after_resample = round(sampling_num * srate / srate_grn)
    if int_srate:
        gcd = np.gcd(int(srate), int(srate_grn))
        p = int(srate) // gcd
        q = int(srate_grn) // gcd
        seis = signal.resample_poly(seis, p, q, axis=-1)
        if not truncate:
            seis = seis[..., :len_after_resample]
    else:
        out = np.zeros((seis.shape[0], 3, len_after_resample))
        for fi in range(seis.shape[0]):
            for ci in range(3):
                out[fi, ci] = resample(
                    seis[fi, ci], srate_old=srate_grn, srate_new=srate, zero_phase=True
                )[:len_after_resample]
        seis = out

    if output_type == "disp":
        seis = np.cumsum(seis, axis=-1) / srate
    elif output_type == "acce":
        kernel = np.array([1.0, -1.0])[:, None]
        out = np.empty_like(seis)
        for fi in range(seis.shape[0]):
            out[fi] = (
                signal.convolve(seis[fi].T, kernel, mode="same", method="auto").T
                / srate
            )
        seis = out

    if n_keep is not None:
        seis = seis[..., :n_keep]
    return seis, tpts_table


if __name__ == "__main__":
    pass
