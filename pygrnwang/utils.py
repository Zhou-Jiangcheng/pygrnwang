import os
import sys
import math
import platform
import subprocess

import numpy as np
import pandas as pd

from .signal_process import linear_interp
from .pytaup import cal_first_p_s
from . import geo


# These three live in geo; they used to be duplicated here verbatim.
# Re-exported so the historical pygrnwang.utils import paths keep working.
cal_max_dist_from_2d_points = geo.cal_max_dist_from_2d_points
create_rotate_z_mat = geo.create_rotate_z_mat
rotate_symmetric_tensor_series = geo.rotate_symmetric_tensor_series


def read_source_array(source_inds, path_input, shift2corner=False, source_shapes=None):
    source_array = None
    for ind_src in range(len(source_inds)):
        source_plane = pd.read_csv(
            str(os.path.join(path_input, "source_plane%d.csv" % source_inds[ind_src])),
            index_col=False,
            header=None,
        ).to_numpy()
        if shift2corner:
            mu_strike = (
                source_plane[source_shapes[ind_src][1], :3] - source_plane[0, :3]
            )
            mu_dip = source_plane[1, :3] - source_plane[0, :3]
            source_plane[:, :3] = source_plane[:, :3] - mu_strike / 2 - mu_dip / 2
        if ind_src == 0:
            source_array = source_plane.copy()
        else:
            source_array = np.concatenate([source_array, source_plane.copy()], axis=0)
    return source_array


def cal_grid(v_min, v_max, delta):
    """Construct the regular grid shared by writers and readers.

    Parameters
    ----------
    v_min : float
        First grid value in any consistent unit.
    v_max : float
        Minimum required terminal value, in the same unit as v_min.
    delta : float
        Positive grid spacing in the same unit as v_min.

    Returns
    -------
    grid : numpy.ndarray

    Raises
    ------
    ZeroDivisionError
        delta is zero.
    ValueError
        Non-finite grid parameters prevent calculation of the sample count.

    Notes
    -----
    """
    n = math.ceil((v_max - v_min) / delta) + 1
    return v_min + np.arange(n) * delta


def group(inp_list, num_in_each_group):
    group_list = []
    for i in range(len(inp_list) // num_in_each_group):
        group_list.append(inp_list[i * num_in_each_group : (i + 1) * num_in_each_group])
    rest = len(inp_list) % num_in_each_group
    if rest != 0:
        group_list.append(inp_list[-rest:])
    return group_list


def shift_green2real_tpts(
    seismograms,
    tpts_table,
    green_before_p,
    srate,
    event_depth_km,
    dist_in_km,
    receiver_depth_km=0,
    model_name="ak135",
):
    first_p, first_s = cal_first_p_s(
        event_depth_km=event_depth_km,
        dist_km=dist_in_km,
        receiver_depth_km=receiver_depth_km,
        model_name=model_name,
    )
    p_count = round(green_before_p * srate)
    s_count = round(
        (tpts_table["s_onset"] - tpts_table["p_onset"] + green_before_p) * srate
    )
    p_count_new = round((first_p - tpts_table["p_onset"] + green_before_p) * srate)
    s_count_new = min(
        len(seismograms[0]),
        round((first_s - tpts_table["p_onset"] + green_before_p) * srate),
    )
    if s_count == p_count or s_count_new == p_count_new:
        return seismograms, first_p, first_s

    n_samples = seismograms.shape[1]
    for i in range(seismograms.shape[0]):
        # own local name: green_before_p must stay the scalar parameter
        before_p_part = seismograms[i][:p_count]
        p_s = linear_interp(seismograms[i][p_count:s_count], s_count_new - p_count_new)
        after_s = seismograms[i][s_count:]
        if len(after_s) > 0:
            after_s = linear_interp(
                after_s, max(0, n_samples - len(before_p_part) - len(p_s))
            )
            row = np.concatenate([before_p_part, p_s, after_s])
        else:
            row = np.concatenate([before_p_part, p_s])
        # the pieces do not always add up to the original length, e.g. when after_s
        # is empty and the real S is earlier than the one in the library
        if len(row) < n_samples:
            row = np.concatenate([row, np.zeros(n_samples - len(row))])
        seismograms[i] = row[:n_samples]

    return seismograms, first_p, first_s


def convert_earth_model_nd2inp(path_nd, path_output):
    """Convert ND numeric rows to numbered backend model input lines.

    Parameters
    ----------
    path_nd : str
        Path to a named-discontinuity text model with either four or six numeric columns as required by the operation.
    path_output : str
        Destination path for the converted model.

    Returns
    -------
    lines : list of str
        Numeric rows prefixed by a one-based row number and terminated by newlines.

    Raises
    ------
    OSError
        The input model cannot be read.

    Notes
    -----
    path_output is retained for API compatibility but is not written by this function. Discontinuity labels are removed.
    """
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            lines_new.append(temp)
    for i in range(len(lines_new)):
        # print(lines_new[i])
        lines_new[i] = "  ".join([str(int(i + 1))] + lines_new[i]) + "\n"  # type: ignore
    # with open(path_output, "w") as fw:
    #     fw.writelines(lines_new)
    return lines_new


def convert_earth_model_nd2nd_without_Q(path_nd, path_output):
    """Write a four-column ND model by removing the final two Q columns.

    Parameters
    ----------
    path_nd : str
        Path to a named-discontinuity text model with either four or six numeric columns as required by the operation.
    path_output : str
        Destination path for the converted model.

    Returns
    -------
    lines : list of str
        The exact converted lines written to path_output, with labels retained.

    Raises
    ------
    OSError
        Input or output cannot be accessed.

    Notes
    -----
    The input must have six numeric columns. Passing a four-column file removes real physical columns and is invalid.
    """
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            lines_new.append([str(float(_)) for _ in temp[:-2]])
            lines_new[i] = "  ".join(lines_new[i]) + "\n"
        else:
            lines_new.append(lines[i].strip() + "\n")
    with open(path_output, "w") as fw:
        fw.writelines(lines_new)
    return lines_new


def read_nd(path_nd, with_Q=False):
    """Read numeric rows from a named-discontinuity model.

    Parameters
    ----------
    path_nd : str
        Path to a named-discontinuity text model with either four or six numeric columns as required by the operation.
    with_Q : bool, optional
        True expects exactly six numeric columns; False expects four. The flag describes the file, it does not remove Q columns. Default: False.

    Returns
    -------
    model : numpy.ndarray
        Shape (N, 4) or (N, 6): depth km, Vp/Vs km/s, density g/cm3, and
        optional dimensionless Qp/Qs. Single-token discontinuity labels are skipped.

    Raises
    ------
    OSError
        The model cannot be read.
    ValueError
        Numeric rows do not match the requested column layout.
    """
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            for j in range(len(temp)):
                lines_new.append(float(temp[j]))
    if with_Q:
        nd_model = np.array(lines_new).reshape(-1, 6)
    else:
        nd_model = np.array(lines_new).reshape(-1, 4)
    return nd_model


def read_material_nd(model_name, depth):
    """Select the first material row at or below the requested depth.

    Parameters
    ----------
    model_name : str
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library.
    depth : float
        Material lookup depth in km, positive down.

    Returns
    -------
    material : numpy.ndarray
        [depth km, Vp km/s, Vs km/s, density g/cm3]; below the model bottom,
        the final row is returned.

    Raises
    ------
    FileNotFoundError
        model_name is neither ak135fc nor an existing model file.

    Notes
    -----
    This is a row selection, not interpolation. model_name accepts only the built-in ak135fc or a four-column no-Q ND file path; ak135 is not a material-lookup alias.
    """
    if model_name == "ak135fc":
        from .ak135fc import s as str_nd

        lines = str_nd.split("\n")
        lines_new = []
        for i in range(len(lines)):
            temp = lines[i].split()
            if len(temp) > 1:
                for j in range(len(temp)):
                    lines_new.append(float(temp[j]))
        nd_model = np.array(lines_new).reshape(-1, 4)
    else:
        if not os.path.isfile(model_name):
            raise FileNotFoundError(
                "model_name must be the built-in 'ak135fc' or a path to an nd file, "
                "got %r" % (model_name,)
            )
        nd_model = read_nd(model_name)
    inds = np.argwhere((nd_model[:, 0] - depth) >= 0)
    # below the bottom of the model: use its deepest layer
    ind = inds[0][0] if len(inds) > 0 else len(nd_model) - 1
    return nd_model[ind]


def read_layerd_material(path_layerd_dat, depth_in_km):
    # thickness, rho, vp, vs, qp, qs
    """Read the material layer containing a depth from a thickness table.

    Parameters
    ----------
    path_layerd_dat : str
        Text table of thickness (m), density, Vp, Vs, Qp, Qs; returned material values retain the file units.
    depth_in_km : float
        Material lookup depth in km.

    Returns
    -------
    material : numpy.ndarray
        One row [thickness, density, Vp, Vs, Qp, Qs] in the file units;
        below the model bottom the final row is returned.

    Raises
    ------
    OSError
        The table cannot be read.
    """
    depth_in_m = depth_in_km * 1e3
    dat = np.loadtxt(path_layerd_dat)
    inds = np.argwhere((np.cumsum(dat[:, 0]) - depth_in_m) >= 0)
    # below the bottom of the model: use its deepest layer
    ind = inds[0][0] if len(inds) > 0 else len(dat) - 1
    return dat[ind]


def create_stf(tau, srate):
    """Sample a normalized squared half-sine source-rate function.

    Parameters
    ----------
    tau : float
        Positive source duration in seconds.
    srate : float
        Positive output sampling rate in Hz.

    Returns
    -------
    stf : numpy.ndarray
        Shape (round(tau*srate)+1,), samples of 2/tau * sin(pi*t/tau)^2,
        in 1/s. Its continuous integral over [0, tau] is one.

    Notes
    -----
    The discrete integral depends on sampling; normalize explicitly when exact discrete convolution normalization is required.
    """
    t = np.linspace(0, tau, round(tau * srate) + 1, endpoint=True)
    stf = (2 / tau) * (np.sin(np.pi * t / tau)) ** 2
    return stf


def group_planes(strike_array):
    """
    It is necessary to ensure that the sub faults on
    the same fault plane have the same strikes!!!
    :param strike_array: numpy array

    Returns:
    np.array: An array containing the lengths of each group.
    """
    # Find the indices where the value changes
    # a[1:] != a[:-1] produces a boolean array that's True
    # at positions where a value differs from its predecessor.
    change_indices = np.where(strike_array[1:] != strike_array[:-1])[0] + 1

    # Include the start and end indices to get boundaries for each group.
    boundaries = np.concatenate(([0], change_indices, [len(strike_array)]))

    # The difference between consecutive boundaries gives the group lengths.
    lengths = np.diff(boundaries)
    return lengths


def reshape_sub_faults(sub_faults, num_strike, num_dip):
    mu_strike = sub_faults[num_dip] - sub_faults[0]
    mu_dip = sub_faults[1] - sub_faults[0]
    sub_faults = sub_faults - mu_strike / 2 - mu_dip / 2
    X: np.ndarray = sub_faults[:, 0]
    Y: np.ndarray = sub_faults[:, 1]
    Z: np.ndarray = sub_faults[:, 2]

    X = X.reshape(num_strike, num_dip)
    Y = Y.reshape(num_strike, num_dip)
    Z = Z.reshape(num_strike, num_dip)

    X = np.concatenate([X, np.array([X[:, -1] + mu_dip[0]]).T], axis=1)
    Y = np.concatenate([Y, np.array([Y[:, -1] + mu_dip[1]]).T], axis=1)
    Z = np.concatenate([Z, np.array([Z[:, -1] + mu_dip[2]]).T], axis=1)

    X = np.concatenate([X, np.array([X[-1, :] + mu_strike[0]])], axis=0)
    Y = np.concatenate([Y, np.array([Y[-1, :] + mu_strike[1]])], axis=0)
    Z = np.concatenate([Z, np.array([Z[-1, :] + mu_strike[2]])], axis=0)
    return X, Y, Z


def call_exe(path_inp, path_finished, name):
    if platform.system() == "Windows":
        name_exe = "%s.exe" % name
        path_exe = os.path.join(sys.exec_prefix, "Scripts", name_exe)
    else:
        name_exe = "%s.bin" % name
        path_exe = os.path.join(sys.exec_prefix, "bin", name_exe)
    proc = subprocess.Popen(
        [path_exe],
        stdin=subprocess.PIPE,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    stdout_bytes, stderr_bytes = proc.communicate(str.encode(path_inp))
    stdout_text = stdout_bytes.decode(errors="ignore")
    stderr_text = stderr_bytes.decode(errors="ignore")
    output = stdout_text + stderr_text
    with open(path_finished, "w", encoding="utf-8") as fw:
        fw.writelines(output)
        return None


def read_tpts_table(path_green, event_depth_km, receiver_depth_km, ind):
    """Read one pair of flat float32 TauP arrival-table entries.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    event_depth_km : float
        Requested source depth in km, positive down.
    receiver_depth_km : float
        Requested receiver depth in km, positive down.
    ind : int
        Zero-based index into the stored distance sequence.

    Returns
    -------
    first_p, first_s : float

    Raises
    ------
    OSError
        Required inputs or outputs cannot be accessed.
    ValueError
        Parameters do not describe a supported grid or observable.

    Notes
    -----


    """
    fr_tp = open(
        os.path.join(
            path_green,
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "tp_table.bin",
        ),
        "rb",
    )
    tp = np.fromfile(file=fr_tp, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_tp.close()

    fr_ts = open(
        os.path.join(
            path_green,
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "ts_table.bin",
        ),
        "rb",
    )
    ts = np.fromfile(file=fr_ts, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_ts.close()
    return float(tp), float(ts)
