from typing import Tuple
import warnings

import numpy as np


def convert_mt_axis(mt, convert_flag):
    """Convert six moment components between NED and spherical RTP axes.

    Parameters
    ----------
    mt : array_like
        Six NED moment components [Mnn, Mne, Mnd, Mee, Med, Mdd] in N m unless a coordinate flag specifies otherwise.
    convert_flag : str
        Either ned2rtp or rtp2ned. RTP ordering is [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp] with outward radial, colatitude and longitude axes.

    Returns
    -------
    converted : list or array_like
        Six components in the target ordering, retaining input physical units.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions. An unrecognized flag currently returns the input unchanged; use only the documented flags.
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


def tensor2full_tensor_matrix(mt, flag="ned") -> np.ndarray:
    """Expand six symmetric tensor components into a full matrix.

    Parameters
    ----------
    mt : array_like
        Six NED moment components [Mnn, Mne, Mnd, Mee, Med, Mdd] in N m unless a coordinate flag specifies otherwise.
    flag : str, optional
        ned uses [Mnn, Mne, Mnd, Mee, Med, Mdd]; rtp uses [Mrr, Mtt, Mpp, Mrt, Mrp, Mtp]. Default: 'ned'.

    Returns
    -------
    matrix : numpy.ndarray
        Shape (3, 3), with rows and columns in the selected axis order.

    Raises
    ------
    ValueError
        flag is neither ned nor rtp.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions.
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
    """Compute scalar moment from the Frobenius norm of a NED tensor.

    Parameters
    ----------
    mt : array_like
        Six NED moment components [Mnn, Mne, Mnd, Mee, Med, Mdd] in N m unless a coordinate flag specifies otherwise.

    Returns
    -------
    moment : float or numpy.ndarray
        sqrt((Mnn^2 + Mee^2 + Mdd^2 + 2*Mne^2 + 2*Mnd^2 + 2*Med^2)/2),
        in the input moment units. Component-first batches return one value per tensor.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions.
    """
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
    """Convert a supported focal mechanism to six NED moment components.

    Parameters
    ----------
    focal_mechanism : array_like
        Either [strike, dip, rake] in degrees; [M0, strike, dip, rake]; six NED components [Mnn, Mne, Mnd, Mee, Med, Mdd]; or [M0, six components]. Three angles imply unit moment; seven entries normalize the six-component shape to M0. Moments are in N m.

    Returns
    -------
    mt : list of float
        [Mnn, Mne, Mnd, Mee, Med, Mdd] in N m.

    Raises
    ------
    ValueError
        The input length is not 3, 4, 6 or 7.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions. A six-component input retains its magnitude. A seven-entry input uses the first value as M0 and normalizes the remaining shape; the shape must have nonzero scalar moment.
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
    """Compute scalar moment from the Frobenius norm of a NED tensor.

    Parameters
    ----------
    mt : array_like
        Six NED moment components [Mnn, Mne, Mnd, Mee, Med, Mdd] in N m unless a coordinate flag specifies otherwise.

    Returns
    -------
    moment : float or numpy.ndarray
        sqrt((Mnn^2 + Mee^2 + Mdd^2 + 2*Mne^2 + 2*Mnd^2 + 2*Med^2)/2),
        in the input moment units. Component-first batches return one value per tensor.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions.
    """
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


def mt2plane(mt):
    """Extract two nodal planes and principal axes from a NED moment tensor.

    Parameters
    ----------
    mt : array_like
        Six NED moment components [Mnn, Mne, Mnd, Mee, Med, Mdd] in N m unless a coordinate flag specifies otherwise.

    Returns
    -------
    result : list
        [plane1, plane2, n1, d1, n2, d2, t, b, p, eigenvalues]. Each plane is
        [strike, dip, rake] in degrees; each vector has shape (3,) in NED.
        Eigenvalues retain the input moment units.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions. For non-double-couple tensors, the planes describe the extracted orientation; they do not reconstruct arbitrary isotropic or CLVD contributions. Repeated eigenvalues make orientation nonunique.
    """
    M = tensor2full_tensor_matrix(mt)
    [eigenvalues, eigenvectors] = np.linalg.eig(M)

    index = eigenvalues.argsort()
    eigenvectors = eigenvectors[:, index]
    eigenvalues = eigenvalues[index]
    p = eigenvectors[:, 0]
    b = eigenvectors[:, 1]
    t = eigenvectors[:, 2]
    n = 1 / np.sqrt(2) * (t + p)
    d = 1 / np.sqrt(2) * (t - p)

    threshold = 1e-15

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
                # warnings.warn(
                #     "n is horizontal. The part directed by n is set as the hanging wall of the fault."
                # )
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
        if cos_lambda > 1:
            lambda_ = np.pi / 2
        elif cos_lambda < -1:
            lambda_ = -np.pi / 2
        else:
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


def plane2mt(M0, strike, dip, rake):
    """Convert a double-couple mechanism to a NED moment tensor.

    Parameters
    ----------
    M0 : float
        Scalar seismic moment in N m.
    strike : float
        Strike angle in degrees clockwise from north.
    dip : float
        Dip angle in degrees from horizontal.
    rake : float
        Rake angle in degrees in the fault plane.

    Returns
    -------
    mt : numpy.ndarray
        Shape (6,), [Mnn, Mne, Mnd, Mee, Med, Mdd], in N m.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions.
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
    """Compute the fault normal and slip unit vectors in NED coordinates.

    Parameters
    ----------
    strike : float
        Strike angle in degrees clockwise from north.
    dip : float
        Dip angle in degrees from horizontal.
    rake : float
        Rake angle in degrees in the fault plane.

    Returns
    -------
    n, d : numpy.ndarray
        Two shape-(3,) vectors; n is oriented upward (nonpositive down component).

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions.
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
    """Compute tension, neutral and pressure axes for a double couple.

    Parameters
    ----------
    strike : float
        Strike angle in degrees clockwise from north.
    dip : float
        Dip angle in degrees from horizontal.
    rake : float
        Rake angle in degrees in the fault plane.

    Returns
    -------
    t, b, p : numpy.ndarray
        Three shape-(3,) NED unit vectors, each oriented into the lower hemisphere.

    Notes
    -----
    See the scientific conventions guide. NED means north, east, down; waveform output uses different axis conventions.
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
