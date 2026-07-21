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
    """
    Merge the CRUST1.0 model (without water and upper_sediments) of a given location
    with the AK135fc model.
    Note, the mantle flag in output nd file is fake, it is only for convenient to call taup.
    Args:
        lat: latitude, unit deg
        lon: longitude, unit deg
        path_crust1: Dir contains crust1.vp,crust1.vs,crust1.rho,crust1.bnds.
            (data in repo https://github.com/jrleeman/Crust1.0.git is recommended)
        path_ak135: Path to ak135fc.nd, without water, with Qp,Qs in the last two cols.
        path_output: Path to output nd file.
        no_low_velo_layer: Ensure that the combination of CRUST1.0 and AK135fc models
            produces no spurious low-velocity zones at the interface, which may remove
            several layers in the upper-mantle of the AK135fc model.
        layered_crust: If True (default), encode every CRUST1.0 layer above
            the mantle as a constant-property layer using repeated interface
            depths. If False, retain linear interpolation between CRUST1.0
            layer-top samples. In both modes, the rows immediately above and
            below the ``mantle`` marker have the same depth.
    Returns:
        nd_new: np.ndarray [[dep, vp, vs, rho, qp, qs], ... ]
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
    with open(path_output, "w") as fw:
        fw.writelines(lines)
    return nd_new


if __name__ == "__main__":
    pass
