from typing import Tuple

import numpy as np


def convert_mt_axis(mt, convert_flag):
    """
    convert moment tensor from one axis to another axis.
    :param mt: moment tensor , if in ned axis, [M11, M12, M13, M22, M23, M33],
                               if in rtp axis, [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp].
    :param convert_flag: 'ned2rtp' or 'rtp2ned'.
    :return:
    """
    if convert_flag == "ned2rtp":
        Mtt = mt[0]
        Mtp = -mt[1]
        Mrt = mt[2]
        Mpp = mt[3]
        Mrp = -mt[4]
        Mrr = mt[5]
        mt = [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]
    elif convert_flag == "rtp2ned":
        M11 = mt[1]
        M12 = -mt[5]
        M13 = mt[3]
        M22 = mt[2]
        M23 = -mt[4]
        M33 = mt[0]
        mt = [M11, M12, M13, M22, M23, M33]
    return mt


def mt2full_mt_matrix(mt, flag="ned") -> np.ndarray:
    """
    create full moment tensor matrix from 6 components.
    :param mt: in NED axis, [M11, M12, M13, M22, M23, M33].
               in rtp axis, [mrr, mtt, mpp, mrt, mrp, mtp]
    :param flag: 'ned'/'rtp'
    :return: full moment tensor matrix,
    in NED axis,
    np.array([[M11, M12, M13],
              [M12, M22, M23],
              [M13, M23, M33]])
    in rtp axis,
    np.array([[mrr, mrt, mrp],
              [mrt, mtt, mtp],
              [mrp, mtp, mpp],
            ])
    """
    mpq = np.zeros((3, 3))
    if flag == "ned":
        mpq[0, 0] = mt[0]
        mpq[0, 1] = mt[1]
        mpq[0, 2] = mt[2]
        mpq[1, 0] = mpq[0, 1]
        mpq[1, 1] = mt[3]
        mpq[1, 2] = mt[4]
        mpq[2, 0] = mpq[0, 2]
        mpq[2, 1] = mpq[1, 2]
        mpq[2, 2] = mt[5]
    elif flag == "rtp":
        mpq[0, 0] = mt[0]
        mpq[0, 1] = mt[3]
        mpq[0, 2] = mt[4]
        mpq[1, 0] = mpq[0, 1]
        mpq[1, 1] = mt[1]
        mpq[1, 2] = mt[5]
        mpq[2, 0] = mpq[0, 2]
        mpq[2, 1] = mpq[1, 2]
        mpq[2, 2] = mt[2]
    else:
        raise ValueError("axis flag wrong")
    return mpq


def moment_from_moment_tensor(mt):
    m0 = np.sqrt(
        1
        / 2
        * (
            mt[0] ** 2
            + 2 * mt[1] ** 2
            + 2 * mt[2] ** 2
            + mt[3] ** 2
            + 2 * mt[4] ** 2
            + mt[5] ** 2
        )
    )
    return m0


def check_convert_fm(focal_mechanism):
    """

    :param focal_mechanism:
    :return: [M11, M12, M13, M22, M23, M33]
    """
    if len(focal_mechanism) == 3:
        mt = plane2mt(1, focal_mechanism[0], focal_mechanism[1], focal_mechanism[2])
        [M11, M12, M13, M22, M23, M33] = list(mt)
    elif len(focal_mechanism) == 4:
        mt = plane2mt(
            focal_mechanism[0],
            focal_mechanism[1],
            focal_mechanism[2],
            focal_mechanism[3],
        )
        [M11, M12, M13, M22, M23, M33] = list(mt)
    elif len(focal_mechanism) == 6:
        [M11, M12, M13, M22, M23, M33] = focal_mechanism
    elif len(focal_mechanism) == 7:
        M0 = focal_mechanism[0]
        temp = np.array(focal_mechanism[1:])
        M0_temp = moment_from_moment_tensor(temp)
        temp = temp / M0_temp
        M11 = M0 * temp[0]
        M12 = M0 * temp[1]
        M13 = M0 * temp[2]
        M22 = M0 * temp[3]
        M23 = M0 * temp[4]
        M33 = M0 * temp[5]
    else:
        raise ValueError("focal mechanism wrong")
    return [M11, M12, M13, M22, M23, M33]


def cal_m0_from_mt(mt):
    m0 = np.sqrt(
        1
        / 2
        * (
            mt[0] ** 2
            + 2 * mt[1] ** 2
            + 2 * mt[2] ** 2
            + mt[3] ** 2
            + 2 * mt[4] ** 2
            + mt[5] ** 2
        )
    )
    return m0


def plane2mt(M0, strike, dip, rake):
    """

    :param M0: scalar moment, unit: Nm
    :param strike: strike angle, unit: degree
    :param dip: dip angle, unit: degree
    :param rake: rake angle, unit: degree
    :return: mt : numpy array
        in NEZ(NED) axis, [M11, M12, M13, M22, M23, M33].
    """
    strike, dip, rake = strike * np.pi / 180, dip * np.pi / 180, rake * np.pi / 180

    sin_strike, cos_strike = np.sin(strike), np.cos(strike)
    sin_2strike, cos_2strike = np.sin(2 * strike), np.cos(2 * strike)
    sin_dip, cos_dip = np.sin(dip), np.cos(dip)
    sin_2dip, cos_2dip = np.sin(2 * dip), np.cos(2 * dip)
    sin_lambda, cos_lambda = np.sin(rake), np.cos(rake)

    M11 = -M0 * (
        sin_dip * cos_lambda * sin_2strike + sin_2dip * sin_lambda * sin_strike**2
    )  # Mtt
    M12 = M0 * (
        sin_dip * cos_lambda * cos_2strike + 1 / 2 * sin_2dip * sin_lambda * sin_2strike
    )  # -Mtp
    M13 = -M0 * (
        cos_dip * cos_lambda * cos_strike + cos_2dip * sin_lambda * sin_strike
    )  # Mrt
    M22 = M0 * (
        sin_dip * cos_lambda * sin_2strike - sin_2dip * sin_lambda * cos_strike**2
    )  # Mpp
    M23 = -M0 * (
        cos_dip * cos_lambda * sin_strike - cos_2dip * sin_lambda * cos_strike
    )  # -Mrp
    M33 = M0 * sin_2dip * sin_lambda  # Mrr

    mt = np.array([M11, M12, M13, M22, M23, M33])
    return mt


def plane2nd(strike, dip, rake) -> Tuple[np.ndarray, np.ndarray]:
    """

    :param strike: unit: degree
    :param dip: unit: degree
    :param rake: unit: degree
    :return: n, np.ndarray
             normal vector of the fault plane, in NED axis.
             d, np.ndarray
             rupture vector on the fault plane, in NED axis.
    """
    strike, dip, rake = strike * np.pi / 180, dip * np.pi / 180, rake * np.pi / 180
    sin_strike, cos_strike = np.sin(strike), np.cos(strike)
    sin_dip, cos_dip = np.sin(dip), np.cos(dip)
    sin_rake, cos_rake = np.sin(rake), np.cos(rake)

    n_nwu = np.array([-sin_dip * sin_strike, -sin_dip * cos_strike, cos_dip])
    n = np.array([n_nwu[0], -n_nwu[1], -n_nwu[2]])

    d_nwu = np.array(
        [
            cos_rake * cos_strike + sin_rake * cos_dip * sin_strike,
            -cos_rake * sin_strike + sin_rake * cos_dip * cos_strike,
            sin_rake * sin_dip,
        ]
    )
    d = np.array([d_nwu[0], -d_nwu[1], -d_nwu[2]])
    if n[2] > 0:  # 保证n朝上
        n = -n
        d = -d
    return n, d


def plane2tbp(strike, dip, rake) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    :param strike: unit: degree
    :param dip: unit: degree
    :param rake: unit: degree
    :return: [np.array(t), np.array(b), np.array(p)]
    """
    n, d = plane2nd(strike, dip, rake)
    t = 1 / np.sqrt(2) * (n + d)
    p = 1 / np.sqrt(2) * (n - d)
    b = np.cross(t, p)
    if t[2] < 0:
        t = -t
    if b[2] < 0:
        b = -b
    if p[2] < 0:
        p = -p
    return t, b, p
