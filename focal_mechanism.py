import warnings
from typing import Tuple, List

import numpy as np

from magnitude import moment_from_moment_tensor
from geo import cartesian_2_spherical
from others import cal_kronecker, fibonacci_sphere


def check_convert_fm(focal_mechanism) -> List:
    """

    :param focal_mechanism:
    :return: [M11, M12, M13, M22, M23, M33]
    """
    if len(focal_mechanism) == 3:
        mt = plane2mt(1, focal_mechanism[0],
                      focal_mechanism[1], focal_mechanism[2])
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


def convert_mt_axis(mt, convert_flag) -> List:
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


def mt2plane(mt):
    """

    :param mt: in NED axis, [M11, M12, M13, M22, M23, M33].
    :return: [[strike1, dip1, rake1], [strike2, dip2, rake2],
    n1, d1, n2, d2, t, b, p, eigenvalues]

    n is the normal vector of the plane, in NED axis.
    d is the rupture vector on the plane, in NED axis.

    n points in the negative direction of D, i.e. upwards

    When the dip angle is 0,
    """
    M = mt2full_mt_matrix(mt)
    [eigenvalues, eigenvectors] = np.linalg.eig(M)

    index = eigenvalues.argsort()
    eigenvectors = eigenvectors[:, index]
    eigenvalues = eigenvalues[index]
    p = eigenvectors[:, 0]
    b = eigenvectors[:, 1]
    t = eigenvectors[:, 2]
    n = 1 / np.sqrt(2) * (t + p)
    d = 1 / np.sqrt(2) * (t - p)

    threshold = 1e-5

    def ignore_small_angle_vector(vector):
        for i in range(3):
            if np.abs(vector[i]) < threshold:
                vector[i] = 0
            if np.abs(vector[i] - 1) < threshold:
                vector[i] = 1
            if np.abs(vector[i] + 1) < threshold:
                vector[i] = -1
        return vector

    n = ignore_small_angle_vector(n)
    d = ignore_small_angle_vector(d)

    # print(eigenvalues)
    # print(eigenvectors)
    # print(t, b, p)
    # print(n,d)

    def nd2plane(n_in, d_in):
        # 保证n朝上
        if n_in[2] > 0:
            n_in = -n_in
            d_in = -d_in
        delta = np.arccos(-n_in[2])
        if np.abs(delta - 0) <= threshold:
            delta = 0
        if np.abs(delta - np.pi / 2) <= threshold:
            delta = np.pi / 2

        if n_in[1] == 0:
            if delta == 0:
                warnings.warn("n is vertical. Strike is set as 0.")
                phi = 0
            elif delta == np.pi / 2:
                warnings.warn(
                    "n is horizontal. The part directed by n is set as the hanging wall of the fault."
                )
                if n_in[0] > 0:
                    phi = np.pi * 3 / 2
                else:
                    phi = np.pi / 2
            else:
                if n_in[0] > 0:
                    phi = np.pi * 3 / 2
                else:
                    phi = np.pi / 2
        else:
            if delta == np.pi / 2:
                warnings.warn(
                    "n is horizontal. The part directed by n is set as the hanging wall of the fault."
                )
                if n_in[0] == 0:
                    if n_in[1] == 1:
                        phi = 0
                    elif n_in[1] == -1:
                        phi = np.pi
                    else:
                        raise ValueError(
                            "n is horizontal,n[0] is 0, but n[1] is not 1 or -1."
                        )
                else:
                    phi = np.arctan(-n_in[0] / n_in[1])
            else:
                phi = np.arctan(-n_in[0] / n_in[1])
        if (n_in[0] <= 0) and (n_in[1] > 0):
            pass
        elif (n_in[0] <= 0) and (n_in[1] < 0):
            phi = phi + np.pi
        elif (n_in[0] > 0) and (n_in[1] < 0):
            phi = phi + np.pi
        elif (n_in[0] > 0) and (n_in[1] > 0):
            phi = phi + 2 * np.pi

        cos_lambda = d_in[0] * np.cos(phi) + d_in[1] * np.sin(phi)
        sin_lambda_cos_delta = d_in[0] * np.sin(phi) - d_in[1] * np.cos(phi)
        lambda_ = np.arccos(cos_lambda)
        if sin_lambda_cos_delta < 0:
            lambda_ = -lambda_
        # elif np.abs(sin_lambda_cos_delta) <= threshold:
        #     if d_in[0] != 0:
        #         lambda_ = np.arctan(d_in[1] / d_in[0])
        #     else:
        #         if np.abs(d_in[1] - 1) <= threshold:
        #             lambda_ = 0
        #         else:
        #             lambda_ = np.arccos(-d_in[1])

        pl = np.array([phi, delta, lambda_])
        pl = pl * 180 / np.pi
        return pl, n_in, d_in

    pl1, n1, d1 = list(nd2plane(n, d))
    pl2, n2, d2 = list(nd2plane(d, n))

    return [pl1, pl2, n1, d1, n2, d2, t, b, p, eigenvalues]


def plane2mt(M0, strike, dip, rake) -> np.ndarray:
    """

    :param M0: scalar moment, unit: Nm
    :param strike: strike angle, unit: degree
    :param dip: dip angle, unit: degree
    :param rake: rake angle, unit: degree
    :return: mt : numpy array
        in NEZ(NED) axis, np.array([M11, M12, M13, M22, M23, M33]).
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
        sin_dip * cos_lambda * cos_2strike + 1 /
        2 * sin_2dip * sin_lambda * sin_2strike
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


def mt2full_mt_matrix(mt) -> np.ndarray:
    """
    create full moment tensor matrix from 6 components.
    :param mt: in NED axis, [M11, M12, M13, M22, M23, M33].
    :return: full moment tensor matrix, in NED axis.
    np.array([[M11, M12, M13],
              [M12, M22, M23],
              [M13, M23, M33]])
    """
    mpq = np.zeros((3, 3))
    mpq[0, 0] = mt[0]
    mpq[0, 1] = mt[1]
    mpq[0, 2] = mt[2]
    mpq[1, 0] = mpq[0, 1]
    mpq[1, 1] = mt[3]
    mpq[1, 2] = mt[4]
    mpq[2, 0] = mpq[0, 2]
    mpq[2, 1] = mpq[1, 2]
    mpq[2, 2] = mt[5]
    return mpq


def cal_A_FP_mt(mt, az, takeoff):
    """
    calculate radiation pattern of P wave
    :param mt: [M11,M12,M13,M22,M23,M33]
    :param az: azimuth in degree
    :param takeoff: takeoff angle in degree
    :return: radiation pattern of P wave at (az, takeoff)
    """
    mt = check_convert_fm(mt)
    mt = np.array(mt) / moment_from_moment_tensor(mt)
    n = create_n_vector(az, takeoff)
    M = mt2full_mt_matrix(mt)
    pattern = np.dot(n, np.dot(M, n.T))[0][0]
    return pattern


def cal_A_FS_mt(mt, az, takeoff):
    """
    calculate radiation pattern of S wave
    Quantitative Seismology (4.29)，p77
    :param mt: [M11,M12,M13,M22,M23,M33]
    :param az: azimuth in degree
    :param takeoff: takeoff angle in degree
    :return: radiation pattern of S wave at (az, takeoff)
    """
    mt = check_convert_fm(mt)
    mt = np.array(mt) / moment_from_moment_tensor(mt)
    n = create_n_vector(az, takeoff)[0]
    M = mt2full_mt_matrix(mt)
    pattern_n = [0, 0, 0]
    for k in range(3):
        for p in range(3):
            for q in range(3):
                pattern_n[k] = (
                    pattern_n[k] - (n[k] * n[p] -
                                    cal_kronecker(k, p)) * n[q] * M[p, q]
                )
    pattern = np.linalg.norm(np.array(pattern_n))
    return pattern


def convert_axis_ned2source(mt, az, takeoff):
    """
    convert axis from NED to source sphere
    get the vector of axis in source sphere which the radiation pattern is normal
    :param mt: [M11,M12,M13,M22,M23,M33]
    :param az: azimuth in degree
    :param takeoff: takeoff angle in degree
    :return:  phi_,theta_ in deg
    """
    [_, _, _, _, _, _, t, b, p, _] = mt2plane(mt)
    x1 = np.array([np.sqrt(2) / 2 * (t + p)])
    x2 = np.array([b])
    x3 = np.array([np.sqrt(2) / 2 * (t - p)])
    A = np.array(np.concatenate([x1.T, x2.T, x3.T], axis=1))
    theta = takeoff * np.pi / 180
    phi = az * np.pi / 180
    r_ned = np.array(
        [np.sin(theta) * np.cos(phi), np.sin(theta)
         * np.sin(phi), np.cos(theta)]
    )
    r_source = np.dot(np.linalg.inv(A), r_ned)
    if (r_source[1] >= 0) and (r_source[0] >= 0):
        phi_ = np.arctan(r_source[1] / r_source[0]) * 180 / np.pi
    elif (r_source[1] >= 0) and (r_source[0] <= 0):
        phi_ = 180 + np.arctan(r_source[1] / r_source[0]) * 180 / np.pi
    elif (r_source[1] <= 0) and (r_source[0] <= 0):
        phi_ = 180 + np.arctan(r_source[1] / r_source[0]) * 180 / np.pi
    elif (r_source[1] <= 0) and (r_source[0] >= 0):
        phi_ = np.arctan(r_source[1] / r_source[0]) * 180 / np.pi
    else:
        raise "az_angle,takeoff error"
    theta_ = np.arccos(r_source[2]) * 180 / np.pi

    return phi_, theta_


def cal_A_FP_DP(strike, dip, rake, az, takeoff):
    """
    calculate radiation pattern of P wave (double couple)
    :param strike: strike in degree
    :param dip: dip in degree
    :param rake: rake in degree
    :param az: azimuth in degree
    :param takeoff: takeoff angle in degree
    :return: pattern
    """
    mt = plane2mt(1, strike, dip, rake)
    phi, theta = convert_axis_ned2source(mt, az, takeoff)
    pattern = np.sin(2 * theta * np.pi / 180) * np.cos(phi * np.pi / 180)
    return pattern


def cal_A_FSV_DP(strike, dip, rake, az, takeoff):
    """
    calculate radiation pattern of P wave (double couple)
    :param strike: strike in degree
    :param dip: dip in degree
    :param rake: rake in degree
    :param az: azimuth in degree
    :param takeoff: takeoff angle in degree
    :return: pattern
    """
    mt = plane2mt(1, strike, dip, rake)
    phi, theta = convert_axis_ned2source(mt, az, takeoff)
    pattern = np.cos(2 * theta * np.pi / 180) * np.cos(phi * np.pi / 180)
    return pattern


def cal_A_FSH_DP(strike, dip, rake, az, takeoff):
    """
    calculate radiation pattern of P wave (double couple)
    :param strike: strike in degree
    :param dip: dip in degree
    :param rake: rake in degree
    :param az: azimuth in degree
    :param takeoff: takeoff angle in degree
    :return: pattern
    """
    mt = plane2mt(1, strike, dip, rake)
    phi, theta = convert_axis_ned2source(mt, az, takeoff)
    pattern = -np.cos(theta * np.pi / 180) * np.sin(phi * np.pi / 180)
    return pattern


def cal_kagan(mt1, mt2):
    mt1 = np.array(mt1) / moment_from_moment_tensor(mt1)
    mt2 = np.array(mt2) / moment_from_moment_tensor(mt2)
    # mt1 和 mt2 相差不大的时候
    kagan = np.arccos(
        np.sum(
            mt1[0] * mt2[0]
            + 2 * mt1[1] * mt2[1]
            + 2 * mt1[2] * mt2[2]
            + mt1[3] * mt2[3]
            + 2 * mt1[4] * mt2[4]
            + mt1[5] * mt2[5]
        )
        / 2
    )
    return kagan


def create_vr_fault_n(az_vr, vr_fault, strike, dip, rake):
    """

    :param az_vr: azimuth of the rupture direction on the surface, unit: degree
    :param vr_fault: rupture velocity on surface, unit: m/s
    :param strike: strike of the fault, unit: degree
    :param dip: dip of the fault, unit: degree
    :param rake: rake of the fault, unit: degree
    :return: vr_fault_n, np.array([vr_fault_x, vr_fault_y, vr_fault_z]), NED axis, unit: m/s
    """
    az_vr = az_vr * np.pi / 180
    n_vr_surf_n = np.array([np.cos(az_vr), np.sin(az_vr), 0])
    n, _ = plane2nd(strike, dip, rake)
    n_vr_fault0 = n_vr_surf_n[0]
    n_vr_fault1 = n_vr_surf_n[1]
    if n[2] != 0:
        n_vr_fault2 = -(n_vr_surf_n[0] * n[0] + n_vr_surf_n[1] * n[1]) / n[2]
    else:
        n_vr_fault2 = 0
    n_vr_fault = np.array([n_vr_fault0, n_vr_fault1, n_vr_fault2])
    n_vr_fault = n_vr_fault / np.linalg.norm(n_vr_fault)
    vr_fault_n = vr_fault * n_vr_fault
    return vr_fault_n


def create_n_vector(az, takeoff):
    """
    create n vector at source sphere from azimuth and takeoff angle.
    axis is NEZ(NED).
    :param az: azimuth in degree.
    :param takeoff: takeoff angle in degree.
    :return:
    """
    theta_ = takeoff * np.pi / 180
    phi_ = az * np.pi / 180
    x1 = np.sin(theta_) * np.cos(phi_)
    x2 = np.sin(theta_) * np.sin(phi_)
    x3 = np.cos(theta_)
    n = np.array([[x1, x2, x3]])
    return n


def convert_vr_surf2vr_fault_n(az_vr, vr_surf, strike, dip, rake):
    """

    :param az_vr: azimuth of the rupture direction on the surface, unit: degree
    :param vr_surf: rupture velocity on surface, unit: m/s
    :param strike: strike of the fault, unit: degree
    :param dip: dip of the fault, unit: degree
    :param rake: rake of the fault, unit: degree
    :return: vr_fault_n, np.array([vr_fault_x, vr_fault_y, vr_fault_z]), NED axis, unit: m/s
    """
    n_vr_fault = create_vr_fault_n(az_vr, 1, strike, dip, rake)
    vr_fault_n = (
        vr_surf / (np.sqrt(n_vr_fault[0] ** 2 +
                   n_vr_fault[1] ** 2)) * n_vr_fault
    )
    return vr_fault_n


def convert_vr_fault2vr_surf_n(az_vr, vr_fault, strike, dip, rake):
    """
    :param az_vr: azimuth of the rupture direction on the surface, unit: degree
    :param vr_fault: rupture velocity on surface, unit: m/s
    :param strike: strike of the fault, unit: degree
    :param dip: dip of the fault, unit: degree
    :param rake: rake of the fault, unit: degree
    :return: vr_surf_n, np.array([vr_fault_x, vr_fault_y, vr_fault_z]), NED axis, unit: m/s
    """
    vr_fault_n = create_vr_fault_n(az_vr, vr_fault, strike, dip, rake)
    vr_surf_n = np.array([vr_fault_n[0], vr_fault_n[1], 0])
    return vr_surf_n


def judge_fm_type(rake):
    """
    the range rake is (-180, 180] deg
    focal mechanism type
    1 normal
    0 strike-slip
    -1 thrust
    :param rake: unit degree
    :return: fm_type
    """
    fm_type = None
    if -180 < rake <= -135:
        fm_type = 0
    elif -135 <= rake < -45:
        fm_type = 1
    elif -45 <= rake < 45:
        fm_type = 0
    elif 45 <= rake < 135:
        fm_type = -1
    elif 135 <= rake <= 180:
        fm_type = 0
    else:
        raise ValueError("rake vaule must in (-180,180] deg")
    return fm_type


def dc_partion(mt):
    s3, s2, s1 = mt2plane(mt)[-1]
    dc = 1 - 2 * abs(s2) / abs(s1)
    return dc


def epsilon_dc_clvd(mt):
    s = mt2plane(mt)[-1]
    s = np.abs(s)
    s = np.sort(s)
    return s[0] / s[-1]


def cal_local_min_max_radiation_pattern(
    mt, az, takeoff, delta_omega=5, delta_n=10, phase="P"
):
    """

    :param mt: ned
    :param az: azimuth in deg
    :param takeoff: takeoff in deg
    :param delta_omega: the radius degree of sub sphere in deg
    :param delta_n: the num of points in sub sphere
    :param phase:
    :return:
    """
    n = create_n_vector(az=az, takeoff=takeoff).flatten()
    omega = np.deg2rad(delta_omega)
    area_ratio = (1 - np.cos(omega)) / 2
    delta_r = 2 * np.sin(omega / 2)
    N = round(np.ceil(delta_n / area_ratio))
    # print("N", N)
    # print("delta_r", delta_r)
    sphere = fibonacci_sphere(1, N)
    FR_sq_local_min = np.inf
    FR_sq_local_max = -np.inf
    for i in range(N):
        ni = sphere[i, :].flatten()
        delta_ri = np.sqrt(np.sum((n - ni) ** 2))
        if delta_ri <= delta_r:
            _, azi, takeoffi = cartesian_2_spherical(*ni)
            azi, takeoffi = np.rad2deg(azi), np.rad2deg(takeoffi)
            if phase == "P":
                FRi = cal_A_FP_mt(mt=mt, az=azi, takeoff=takeoffi)
            elif phase == "S":
                FRi = cal_A_FS_mt(mt=mt, az=azi, takeoff=takeoffi)
            elif phase == "SH":
                strike, dip, rake = mt2plane(mt=mt)[0]
                FRi = cal_A_FSH_DP(strike, dip, rake, az=azi, takeoff=takeoffi)
            elif phase == "SV":
                strike, dip, rake = mt2plane(mt=mt)[0]
                FRi = cal_A_FSV_DP(strike, dip, rake, az=azi, takeoff=takeoffi)
            else:
                raise ValueError("phase type wrong")
            if FRi**2 < FR_sq_local_min:
                FR_sq_local_min = FRi**2
            if FRi**2 > FR_sq_local_max:
                FR_sq_local_max = FRi**2
    return np.sqrt(FR_sq_local_min), np.sqrt(FR_sq_local_max)


def cal_local_mean_radiation_pattern(
    mt, az, takeoff, delta_omega=5, delta_n=10, phase="P"
):
    """

    :param mt: ned
    :param az: azimuth in deg
    :param takeoff: takeoff in deg
    :param delta_omega: the radius degree of sub sphere in deg
    :param delta_n: the num of points in sub sphere
    :param phase:
    :return:
    """
    n = create_n_vector(az=az, takeoff=takeoff).flatten()
    omega = np.deg2rad(delta_omega)
    area_ratio = (1 - np.cos(omega)) / 2
    delta_r = 2 * np.sin(omega / 2)
    N = round(np.ceil(delta_n / area_ratio))
    # print("N", N)
    # print("delta_r", delta_r)
    sphere = fibonacci_sphere(1, N)
    FR_sq_local_mean = 0
    num = 0
    for i in range(N):
        ni = sphere[i, :].flatten()
        delta_ri = np.sqrt(np.sum((n - ni) ** 2))
        if delta_ri <= delta_r:
            _, azi, takeoffi = cartesian_2_spherical(*ni)
            azi, takeoffi = np.rad2deg(azi), np.rad2deg(takeoffi)
            if phase == "P":
                FRi = cal_A_FP_mt(mt=mt, az=azi, takeoff=takeoffi)
            elif phase == "S":
                FRi = cal_A_FS_mt(mt=mt, az=azi, takeoff=takeoffi)
            elif phase == "SH":
                strike, dip, rake = mt2plane(mt=mt)[0]
                FRi = cal_A_FSH_DP(strike, dip, rake, az=azi, takeoff=takeoffi)
            elif phase == "SV":
                strike, dip, rake = mt2plane(mt=mt)[0]
                FRi = cal_A_FSV_DP(strike, dip, rake, az=azi, takeoff=takeoffi)
            else:
                raise ValueError("phase type wrong")
            FR_sq_local_mean = FR_sq_local_mean + FRi**2
            num = num + 1
    FR_local_mean = np.sqrt(FR_sq_local_mean / num)
    # print("num", num)
    # print("FR_local_sq_mean", FR_sq_local_mean / num)
    return FR_local_mean


def cal_mean_radiation_pattern_sq(mt, phase="P", N=100):
    sphere = fibonacci_sphere(1, N)
    FR_sq_mean = 0
    for i in range(N):
        _, az, takeoff = cartesian_2_spherical(*sphere[i, :])
        az, takeoff = np.rad2deg(az), np.rad2deg(takeoff)
        if phase == "P":
            FRi = cal_A_FP_mt(mt, az, takeoff)
        elif phase == "S":
            FRi = cal_A_FS_mt(mt, az, takeoff)
        else:
            raise ValueError("phase wrong")
        FR_sq_mean = FR_sq_mean + FRi**2
    FR_sq_mean = FR_sq_mean / N
    return FR_sq_mean


def cal_rescaled_FP(cf, ray, sta_info):
    FP = cal_A_FP_mt(
        mt=cf.focal_mechanism,
        az=sta_info["az"],
        takeoff=ray.takeoff_angle(),
    )
    if cf.rescale_rp_type == 0:
        FP_rescale = cal_local_mean_radiation_pattern(
            mt=cf.focal_mechanism,
            az=sta_info["az"],
            takeoff=ray.takeoff_angle(),
            delta_omega=cf.rescale_rp_paras[0],
            delta_n=int(cf.rescale_rp_paras[1]),
            phase="P",
        )
    elif cf.rescale_rp_type == 1:
        FP_rescale, _ = cal_local_min_max_radiation_pattern(
            mt=cf.focal_mechanism,
            az=sta_info["az"],
            takeoff=ray.takeoff_angle(),
            delta_omega=cf.rescale_rp_paras[0],
            delta_n=int(cf.rescale_rp_paras[1]),
            phase="P",
        )
    elif cf.rescale_rp_type == 2:
        _, FP_rescale = cal_local_min_max_radiation_pattern(
            mt=cf.focal_mechanism,
            az=sta_info["az"],
            takeoff=ray.takeoff_angle(),
            delta_omega=cf.rescale_rp_paras[0],
            delta_n=int(cf.rescale_rp_paras[1]),
            phase="P",
        )
    elif cf.rescale_rp_type == 3:
        FP_min, FP_max = cal_local_min_max_radiation_pattern(
            mt=cf.focal_mechanism,
            az=sta_info["az"],
            takeoff=ray.takeoff_angle(),
            delta_omega=cf.rescale_rp_paras[0],
            delta_n=int(cf.rescale_rp_paras[1]),
            phase="P",
        )
        if cf.fm_rup_type[:2] == "SS":
            FP_rescale = FP_max
        else:
            FP_rescale = FP_min
    elif cf.rescale_rp_type == 4:
        FP_min, FP_max = cal_local_min_max_radiation_pattern(
            mt=cf.focal_mechanism,
            az=sta_info["az"],
            takeoff=ray.takeoff_angle(),
            delta_omega=cf.rescale_rp_paras[0],
            delta_n=int(cf.rescale_rp_paras[1]),
            phase="P",
        )
        if cf.fm_rup_type[:2] == "SS":
            FP_rescale = FP_max
        else:
            FP_rescale = FP
    else:
        raise ValueError(
            "rescale_rp_type is %s, rescale_rp_type wrong!" % str(
                cf.rescale_rp_type)
        )
    return FP, FP_rescale


if __name__ == "__main__":
    pass
