import warnings

import numpy as np

from .crust1 import CrustModel
from .utils import read_nd


def create_nd_by_crust1_ak135(
    lat: float,
    lon: float,
    path_crust1: str,
    path_ak135: str,
    path_output: str,
    no_low_velo_layer: bool = False,
    layered_crust: bool = True,
):
    """Join a location-specific CRUST1.0 crust to an AK135 mantle model.

    Parameters
    ----------
    lat : float
        Latitude in degrees.
    lon : float
        Longitude in degrees.
    path_crust1 : str
        Directory containing crust1.vp, crust1.vs, crust1.rho and crust1.bnds.
    path_ak135 : str
        Six-column no-water AK135 ND file with Qp and Qs in the last columns.
    path_output : str
        Destination path for the converted model.
    no_low_velo_layer : bool, optional
        Remove conflicting shallow mantle rows to avoid an artificial low-velocity join. Default: False.
    layered_crust : bool, optional
        True repeats interface depths to encode constant-property CRUST1 layers; False uses linear interpolation between layer tops. Default: True.

    Returns
    -------
    model : numpy.ndarray
        Shape (N, 6): depth km, Vp/Vs km/s, density g/cm3, Qp/Qs.
        The corresponding ND file is also written to path_output.

    Raises
    ------
    OSError
        A required model file is unavailable.
    ValueError
        The CRUST1 columns lack suitable crust/mantle rows or increasing layer depths.

    Notes
    -----
    Water and upper sediments are omitted. CRUST1 rows use Qp=927.34 and Qs=599.99. The mantle label is placed for TauP compatibility; verify that it represents the intended model discontinuity.
    """
    crust1 = CrustModel(path_crust1)
    crust1_dict = crust1.get_point(lat, lon)

    dep_list = []
    vp_list = []
    vs_list = []
    rho_list = []
    qp_list = []
    qs_list = []

    for key, val in crust1_dict.items():
        if key == "water" or key == "upper_sediments":
            continue
        vp_list.append(crust1_dict[key][0])
        vs_list.append(crust1_dict[key][1])
        rho_list.append(crust1_dict[key][2])
        dep_list.append(-crust1_dict[key][4])
        qp_list.append(927.34)
        qs_list.append(599.99)

    nd_crust1 = np.concatenate(
        [
            np.array([dep_list]),
            np.array([vp_list]),
            np.array([vs_list]),
            np.array([rho_list]),
            np.array([qp_list]),
            np.array([qs_list]),
        ]
    ).T
    if nd_crust1[0, 0] < 0:
        nd_crust1[:, 0] = nd_crust1[:, 0] - nd_crust1[0, 0]
    else:
        nd_crust1[0, 0] = 0

    nd_ak135 = read_nd(path_ak135, True)
    N_ak135 = len(nd_ak135)
    ind_cut = 0
    for i in range(N_ak135):
        if nd_ak135[i, 0] >= np.max(nd_crust1[:, 0]):
            ind_cut = i
            break
    if no_low_velo_layer:
        ind_cut_new = ind_cut
        for i in range(ind_cut_new, N_ak135):
            if nd_ak135[i, 1] >= np.max(nd_crust1[:, 1]):
                ind_cut_new = i
                break
        for i in range(ind_cut_new, N_ak135):
            if nd_ak135[i, 2] >= np.max(nd_crust1[:, 2]):
                ind_cut_new = i
                break
        for i in range(ind_cut_new, N_ak135):
            if nd_ak135[i, 3] >= np.max(nd_crust1[:, 3]):
                ind_cut_new = i
                break
        if dep_list[-1] > 660:
            warnings.warn(
                "The cutoff depth exceeds 660km, "
                "ignoring parameter no_low_velo_layer."
            )
        else:
            ind_cut = ind_cut_new

    if len(nd_crust1) < 2:
        raise ValueError("CRUST1.0 model must contain crust and mantle rows")

    # CRUST1.0 supplies one constant-property value per layer at the layer
    # top. An nd file, however, linearly interpolates between adjacent rows.
    # Repeat each interface depth to preserve CRUST1.0's stepwise layering.
    crust_rows = nd_crust1[:-1]
    crust_mantle_row = nd_crust1[-1].copy()
    moho_depth = float(crust_mantle_row[0])
    if np.any(np.diff(crust_rows[:, 0]) <= 0) or moho_depth <= crust_rows[-1, 0]:
        raise ValueError("CRUST1.0 layer-top depths must increase toward the mantle")

    if layered_crust:
        rows_above_mantle = []
        for i, top_row in enumerate(crust_rows):
            bottom_depth = (
                crust_rows[i + 1, 0] if i + 1 < len(crust_rows) else moho_depth
            )
            bottom_row = top_row.copy()
            bottom_row[0] = bottom_depth
            rows_above_mantle.extend((top_row.copy(), bottom_row))
        rows_above_mantle = np.asarray(rows_above_mantle)
    else:
        bottom_row = crust_rows[-1].copy()
        bottom_row[0] = moho_depth
        rows_above_mantle = np.vstack((crust_rows, bottom_row))

    mantle_rows = nd_ak135[ind_cut:, :].copy()
    if not len(mantle_rows):
        raise ValueError("AK135 model has no rows at or below the Moho")
    if not np.isclose(mantle_rows[0, 0], moho_depth):
        # When no_low_velo_layer skips deeper into AK135, retain the CRUST1.0
        # mantle value at the Moho so the named boundary remains well formed.
        crust_mantle_row[0] = moho_depth
        mantle_rows = np.vstack((crust_mantle_row, mantle_rows))

    mantle_index = len(rows_above_mantle)
    nd_new = np.vstack((rows_above_mantle, mantle_rows))

    fluid_indices = np.flatnonzero(nd_new[:, 2] == 0)
    if not len(fluid_indices) or fluid_indices[-1] + 1 >= len(nd_new):
        raise ValueError("AK135 model must contain outer- and inner-core rows")
    boundary_labels = {
        mantle_index: "mantle",
        int(fluid_indices[0]): "outer-core",
        int(fluid_indices[-1] + 1): "inner-core",
    }

    lines = []
    for i, row in enumerate(nd_new):
        if i in boundary_labels:
            lines.append(boundary_labels[i] + "\n")
        lines.append(" ".join("%12.5f" % float(value) for value in row[:6]) + "\n")
    lines.append('\n')
    with open(path_output, "w") as fw:
        fw.writelines(lines)
    return nd_new


if __name__ == "__main__":
    pass
