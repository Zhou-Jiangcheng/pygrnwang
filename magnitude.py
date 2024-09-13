import numpy as np


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


def moment_from_stf(stf, srate):
    m0 = np.trapz(stf, dx=1 / srate)
    return m0


def moment_from_moment_mag(mw):
    m0 = 10 ** (3 / 2 * mw + 9.1)
    return m0


def moment_magnitude(m0):
    Mw = 2 / 3 * (np.log10(m0) - 9.1)
    return Mw


def energy_magnitude(energy):
    Me = 2 / 3 * np.log10(energy) - 2.9
    return Me


def energy_from_energy_mag(me):
    Er = 10 ** (3 / 2 * (me + 2.9))
    return Er


def moment_magnitude_from_fault_length(fault_length):
    """
    calculate moment magnitude from fault length
     (Wells and Coppersmith, 1994)
    :param fault_length: fault length, unit: m
    :return: moment magnitude
    """
    Mw = 1.16 * np.log10(fault_length / 1000) + 5.08
    return Mw


def fault_length_from_moment_magnitude(mw):
    """
    calculate fault length from moment magnitude
    (Wells and Coppersmith, 1994)
    :param mw: moment magnitude
    :return: fault length, unit: m
    """
    L = 10 ** ((mw - 5.08) / 1.16)
    L = L * 1000
    return L


def cal_m0_from_sub_stfs(sub_stfs, srate):
    """
    calculate scalar moment from sub-fault source time functions
    :param sub_stfs: sub-fault source time functions, np.array, (N,M,length_in_counts) unit: Nm/s
    :param srate: sampling rate, unit: Hz
    :return: source_m0, scalar moment, unit: Nm
    """
    m0 = 0
    for i in range(sub_stfs.shape[0]):
        for j in range(sub_stfs.shape[1]):
            m0 += moment_from_stf(sub_stfs[i, j, :], srate)
    return m0


if __name__ == "__main__":
    L_ = fault_length_from_moment_magnitude(8)
    print(L_)
