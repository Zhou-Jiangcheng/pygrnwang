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

    nd_new = np.concatenate([nd_crust1, nd_ak135[ind_cut:, :]])
    inds = np.argwhere(nd_new[:, 2] == 0)
    ind_cmb = inds[0][0]
    ind_icocb = inds[-1][0]
    ind_660 = np.argwhere(nd_new[:, 0] == 660)[0][0]

    lines = []
    for i in range(len(nd_new)):
        line = ""
        for j in range(6):
            line = line + "%12.5f " % float(nd_new[i, j])
        line = line[:-1] + "\n"
        lines.append(line)
    # fake crust-mantle boundary, only for convenient to taup
    lines.insert(ind_660 + 1, "mantle\n")
    lines.insert(ind_cmb + 1, "outer-core\n")
    lines.insert(ind_icocb + 3, "inner-core\n")
    with open(path_output, "w") as fw:
        fw.writelines(lines)
    return nd_new


if __name__ == "__main__":
    pass
