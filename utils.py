import numpy as np


def get_number_in_line(line):
    """
    获取这一行里的数字
    """
    numbers = []
    for item in line.split():
        if item != "":
            try:
                numbers.append(float(item))
            except ValueError:
                pass
    return numbers


def group(inp_list, num_in_each_group):
    group_list = []
    for i in range(len(inp_list) // num_in_each_group):
        group_list.append(inp_list[i * num_in_each_group : (i + 1) * num_in_each_group])
    rest = len(inp_list) % num_in_each_group
    if rest != 0:
        group_list.append(inp_list[-rest:])
    return group_list


def check_convert_fm(focal_mechanism):
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
        M0_temp = cal_m0_from_mt(temp)
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


def rotate_rtz_to_enz(az_in_deg, r, t, z):
    az = np.deg2rad(az_in_deg)
    e = r * np.sin(az) + t * np.cos(az)
    n = r * np.cos(az) - t * np.sin(az)
    return [e, n, z]


def find_nearest_dichotomy(value, value_list):
    # 二分法找到列表中最接近给定值的项
    start, end = 0, len(value_list) - 1
    while start <= end:
        mid_index = (start + end) // 2
        # print(value, value_list[mid_index], value_list[mid_index + 1])
        if value == value_list[start]:
            return [value_list[start], start]
        elif value == value_list[end]:
            return [value_list[end], end]
        elif value_list[mid_index] <= value <= value_list[mid_index + 1]:
            if (value - value_list[mid_index]) <= (value_list[mid_index + 1] - value):
                return [value_list[mid_index], mid_index]
            else:
                return [value_list[mid_index + 1], mid_index + 1]
        elif value > value_list[mid_index]:
            start = mid_index + 1
        elif value < value_list[mid_index]:
            end = mid_index - 1
    else:
        raise ValueError("value " + str(value) + " not in the value_list range")


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


def convert_earth_model_nd2inp(path_nd, path_output):
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            lines_new.append(temp)
    for i in range(len(lines_new)):
        # print(lines_new[i])
        lines_new[i] = "  ".join([str(i + 1)] + lines_new[i]) + "\n"
    with open(path_output, "w") as fw:
        fw.writelines(lines_new)
    return lines_new


if __name__ == "__main__":
    path_nd = "/home/zjc/python_works/pygrnwang/ak135.nd"
    path_output = "ak135.dat"
    convert_earth_model_nd2inp(path_nd, path_output)
