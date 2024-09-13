import os.path
import re
import numpy as np


def is_float(s):
    try:
        float(s)
        return True
    except ValueError:
        return False


def is_int(s):
    try:
        int(s)
        return True
    except ValueError:
        return False


def create_lt_mat(vector):
    """
    Create a left-triangular matrix from a vector.
    :param vector:
    :return: matrix A, shape=(len(vector),len(vector))
    """
    N = len(vector)
    A = np.zeros([N, N])
    for n in range(N):
        A[n] = np.concatenate([vector[: n + 1][::-1], np.zeros(N - (n + 1))])
    return A


def create_lt_mat_complex(vector):
    """
    Create a left-triangular matrix from a vector.
    :param vector:
    :return: matrix A, shape=(len(vector),len(vector))
    """
    N = len(vector)
    A = np.zeros([N, N], dtype=complex)
    for n in range(N):
        A[n] = np.concatenate([vector[: n + 1][::-1], np.zeros(N - (n + 1))])
    return A


def create_conv_mat(vector):
    """
    Create a full conv matrix from a vector.
    :param vector:
    :return: matrix A, shape=(2*len(vector)-1,len(vector))
    """
    N = len(vector)
    A = np.zeros([2 * N - 1, N])
    for n in range(N):
        A[n: N + n, n] = vector
    return A


def create_cyclic_conv_mat(vector):
    """
    Create a full cyclic conv matrix from a vector.
    The vector should be padded with zeros already.
    :param vector:
    :return: matrix A, shape=(len(vector),len(vector))
    """
    N = len(vector)
    A = np.zeros([N, N])
    for n in range(N):
        A[:, n] = np.roll(vector, n)
    return A


def find_nearest_dichotomy(value, value_list):
    """

    :param value:
    :param value_list:
    :return: [value_list[ind], ind]
    """
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
        raise ValueError("value " + str(value) +
                         " not in the value_list range")


def normal_distribution_2D(mu1, mu2, sigma1, sigma2, rho, x, y):
    """
    calculate the value of 2D normal distribution,
    distribution.shape= (len(x), len(y))

    :param mu1:
    :param mu2:
    :param sigma1:
    :param sigma2:
    :param rho:
    :param x: np.array, X-axis coordinate point
    :param y: np.array, Y-axis coordinate point
    :return:
    distribution: np.array, shape=(len(x), len(y))
    """
    X, Y = np.meshgrid(x, y)
    distribution = (
        1
        / (2 * np.pi * sigma1 * sigma2 * np.sqrt(1 - rho**2))
        * np.exp(
            -1
            / (2 * (1 - rho**2))
            * (
                (X - mu1) ** 2 / sigma1**2
                + (Y - mu2) ** 2 / sigma2**2
                - 2 * rho * (X - mu1) * (Y - mu2) / (sigma1 * sigma2)
            )
        )
    )
    return distribution


def fibonacci_sphere(r0=1, N=100):
    points = np.zeros((N, 3))
    phi = np.pi * (3 - np.sqrt(5))  # golden angle in radians
    y = 1 - np.arange(N) / (N - 1) * 2
    radius = np.sqrt(1 - y**2)
    theta = phi * np.arange(N)
    points[:, 0] = radius * np.cos(theta)
    points[:, 1] = y
    points[:, 2] = radius * np.sin(theta)
    points = points * r0
    return points


def get_file_path(path_dir, net, sta, loc, cha):
    # if loc == "0":
    #     loc = "00"
    # elif loc == "1":
    #     loc = "01"
    path = os.path.join(path_dir, "%s.%s.%s.%s" % (net, sta, loc, cha))
    if os.path.exists(path):
        return path
    else:
        if cha[-1] == "1":
            cha[-1] = "E"
        if cha[-1] == "2":
            cha[-1] = "N"
        path = os.path.join(path_dir, "%s.%s.%s.%s" % (net, sta, loc, cha))
        if os.path.exists(path):
            return path
        else:
            raise ValueError("%s\nthis path does not exist" % path)


def convert_str2float_list(strings: str):
    list_temp = strings.replace("[", "").replace(
        "]", "").replace(" ", "").split(",")
    list_final = []
    for i in range(len(list_temp)):
        if list_temp[i] != "":
            list_final.append(float(list_temp[i]))
    return list_final


def convert_str2int_list(strings: str):
    list_temp = strings.replace("[", "").replace(
        "]", "").replace(" ", "").split(",")
    list_final = []
    for i in range(len(list_temp)):
        if list_temp[i] != "":
            list_final.append(int(list_temp[i]))
    return list_final


def convert_multi_str2float_list(strings):
    strings = strings.replace("\n", ",")
    strings = strings[1:-1]
    list_temp = strings.replace(" ", "").split(",")
    list_final = []
    for i in range(len(list_temp)):
        if "[" in list_temp[i]:
            list_final.append([])
            list_temp[i] = list_temp[i][1:]
        if "]" in list_temp[i]:
            list_temp[i] = list_temp[i][:-1]
        if list_temp[i] != "":
            list_final[-1].append(float(list_temp[i]))
    return list_final


def convert_str2str_list(strings: str):
    list_temp = (
        strings.replace("[", "")
        .replace("]", "")
        .replace(" ", "")
        .replace("\n", "")
        .split(",")
    )
    list_final = []
    for i in range(len(list_temp)):
        if list_temp[i] != "":
            list_final.append(list_temp[i])
    return list_final


def convert_str2item_list(strings: str):
    list_temp = (
        strings.replace("[", "")
        .replace("]", "")
        .replace(" ", "")
        .replace("\n", "")
        .split(",")
    )
    list_final = []
    vdict = {"None": None, "True": True, "False": False}
    for i in range(len(list_temp)):
        if list_temp[i] != "":
            if is_int(list_temp[i]):
                list_final.append(int(list_temp[i]))
            elif is_float(list_temp[i]):
                list_final.append(float(list_temp[i]))
            elif list_temp[i] in vdict.keys():
                list_final.append(vdict[list_temp[i]])
            else:
                list_final.append(list_temp[i])
    return list_final


def get_number_in_line(line):
    numbers = re.findall(
        r"[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?", line)
    # print(numbers)
    numbers = [float(item) for item in numbers if item != ""]
    return numbers


def reshape_surface_3D(data: np.ndarray) -> np.ndarray:
    """
    to reshape slips for ax.plot_surface
    The new slips point is located on the frame of the plane
    :param data:
    :return:
    """
    data_new = np.zeros([data.shape[0] + 1, data.shape[1] + 1, 3])
    P00 = data[0, 0, :] + 1 / 2 * (data[0, 0, :] - data[1, 1, :])
    P0M = data[0, -1, :] + 1 / 2 * (data[0, -1, :] - data[1, -2, :])
    PN0 = data[-1, 0, :] + 1 / 2 * (data[-1, 0, :] - data[-2, 1, :])
    PNM = data[-1, -1, :] + 1 / 2 * (data[-1, -1, :] - data[-2, -2, :])
    for i in range(3):
        data_new[0, :, i] = np.linspace(P00[i], P0M[i], data_new.shape[1])
        data_new[-1, :, i] = np.linspace(PN0[i], PNM[i], data_new.shape[1])

        for m in range(data_new.shape[1]):
            data_new[:, m, i] = np.linspace(
                data_new[0, m, i], data_new[-1, m, i], data_new.shape[0]
            )

    return data_new


def cal_kronecker(i, j):
    if i == j:
        return 1
    else:
        return 0


def cal_max_dist_from_2d_points(A: np.ndarray, B: np.ndarray):
    """

    :param A: (m,2)
    :param B: (n,2)
    :return: max_distance
    """
    # Calculate the differences in each dimension (broadcasting)
    differences = A[:, np.newaxis, :] - B[np.newaxis, :, :]

    # Square the differences and sum across columns (to get squared distances)
    squared_distances = np.sum(differences**2, axis=2)

    # Take the square root to get Euclidean distances
    distances = np.sqrt(squared_distances)

    # Find the maximum distance
    max_distance = np.max(distances)
    return max_distance


def cal_weight_std(data, mean, weights):
    return np.sqrt(np.sum((data - mean) ** 2 * weights))
