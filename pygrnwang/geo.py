import numpy as np
import pandas as pd

d2m = 111194.92664455874
d2km = 111.19492664455874


def rotate_2d_points(points: np.ndarray, degree: float) -> np.ndarray:
    """
    Rotate a 2D point counterclockwise by a given number of degrees.

    :param points: np.array([[x1, y1],[x2,y2],...])
    :param degree: The angle in degrees by which the point is to be rotated (anti-clockwise).
    :return: rotated_points.
    """
    # Convert degrees to radians
    radians = np.radians(degree)

    # Rotation matrix for counterclockwise rotation
    rotation_matrix = np.array(
        [[np.cos(radians), -np.sin(radians)], [np.sin(radians), np.cos(radians)]]
    ).T

    # Rotated point
    rotated_points = np.dot(points, rotation_matrix)
    # print("rotated_points", rotated_points.shape)

    return rotated_points


def rotate_rtz_to_enz(az_in_deg, r, t, z):
    az = np.deg2rad(az_in_deg)
    e = r * np.sin(az) - t * np.cos(az)
    n = r * np.cos(az) + t * np.sin(az)
    return np.array([e, n, z])


def rotate_vector_x(v, gamma):
    """

    :param v:
    :param gamma: deg
    :return:
    """
    gamma = np.deg2rad(gamma)
    R = np.array(
        [
            [1, 0, 0],
            [0, np.cos(gamma), -np.sin(gamma)],
            [0, np.sin(gamma), np.cos(gamma)],
        ]
    )
    rotated_vec = np.dot(R, v)

    return rotated_vec


def rotate_vector_y(v, gamma):
    """

    :param v:
    :param gamma: deg
    :return:
    """
    gamma = np.deg2rad(gamma)
    R = np.array(
        [
            [np.cos(gamma), 0, np.sin(gamma)],
            [0, 1, 0],
            [-np.sin(gamma), 0, np.cos(gamma)],
        ]
    )

    rotated_vec = np.dot(R, v)

    return rotated_vec


def rotate_vector_z(v, gamma):
    """

    :param v:
    :param gamma: deg
    :return:
    """
    gamma = np.deg2rad(gamma)
    R = np.array(
        [
            [np.cos(gamma), -np.sin(gamma), 0],
            [np.sin(gamma), np.cos(gamma), 0],
            [0, 0, 1],
        ]
    )

    rotated_vec = np.dot(R, v)

    return rotated_vec


def spherical_2_cartesian(r, phi, theta):
    """

    :param r:
    :param phi: deg
    :param theta: deg
    :return:
    """
    phi, theta = np.deg2rad(phi), np.deg2rad(theta)
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    return np.array([x, y, z])


def cartesian_2_spherical(x, y, z):
    """

    :param x:
    :param y:
    :param z:
    :return: r, phi, theta ( , rad, rad)
    """
    r = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arccos(z / r)
    if x != 0:
        phi_x = np.arctan(y / x)
        if (x > 0) and (y >= 0):
            phi = phi_x
        elif (x < 0) and (y >= 0):
            phi = phi_x + np.pi
        elif (x < 0) and (y < 0):
            phi = phi_x + np.pi
        elif (x > 0) and (y < 0):
            phi = 2 * np.pi + phi_x
        else:
            raise ValueError
    else:
        if y > 0:
            phi = np.pi / 2
        elif y < 0:
            phi = 3 / 2 * np.pi
        else:
            phi = 0
    return r, phi, theta


def geo_2_r_earth(lat, lon, dep, r0=6371000):
    """

    :param lat: deg
    :param lon: deg
    :param dep: m
    :param r0: m
    :return: r_earth (unit m)
    """
    lat, lon = np.deg2rad(lat), np.deg2rad(lon)
    r_earth = np.array(
        [
            [
                (r0 - dep) * np.cos(lat) * np.cos(lon),
                (r0 - dep) * np.cos(lat) * np.sin(lon),
                (r0 - dep) * np.sin(lat),
            ]
        ]
    ).flatten()
    return r_earth


def r_earth_2_geo(r_earth, r0=6371000):
    """

    :param r_earth: m
    :param r0: m
    :return: lat, lon, depth (deg, deg, m)
    """
    r, lon, co_lat = cartesian_2_spherical(r_earth[0], r_earth[1], r_earth[2])
    depth = r0 - r
    lon = np.rad2deg(lon)
    if lon > 180:
        lon = lon - 360
    lat = np.rad2deg(np.pi / 2 - co_lat)
    return np.array([lat, lon, depth])


def convert_axis_delta_geo2ned(lat0, lon0, dep0, lat1, lon1, dep1):
    """

    :param lat0:
    :param lon0:
    :param dep0:
    :param lat1:
    :param lon1:
    :param dep1:
    :return: r (in ned axis), unit m
    """
    r_earth0 = geo_2_r_earth(lat0, lon0, dep0)
    r_earth1 = geo_2_r_earth(lat1, lon1, dep1)
    delta_r_earth = r_earth1 - r_earth0
    lat0 = np.deg2rad(lat0)
    lon0 = np.deg2rad(lon0)
    A = np.array(
        [
            [-np.sin(lat0) * np.cos(lon0), -np.sin(lat0) * np.sin(lon0), np.cos(lat0)],
            [-np.sin(lon0), np.cos(lon0), 0],
            [-np.cos(lat0) * np.cos(lon0), -np.cos(lat0) * np.sin(lon0), -np.sin(lat0)],
        ]
    )
    delta_r_ned = np.dot(A, delta_r_earth).flatten()
    return delta_r_ned


def convert_axis_delta_ned2geo(lat0, lon0, dep0, r_ned):
    """

    :param lat0: deg
    :param lon0: deg
    :param dep0: m
    :param r_ned: np.ndarray, m
    :return: lat, lon, depth (deg, deg, m)
    """
    r_ned = np.array(r_ned).flatten()
    r_earth0 = geo_2_r_earth(lat0, lon0, dep0)
    lat0, lon0 = np.deg2rad(lat0), np.deg2rad(lon0)
    A = np.array(
        [
            [-np.sin(lat0) * np.cos(lon0), -np.sin(lat0) * np.sin(lon0), np.cos(lat0)],
            [-np.sin(lon0), np.cos(lon0), 0],
            [-np.cos(lat0) * np.cos(lon0), -np.cos(lat0) * np.sin(lon0), -np.sin(lat0)],
        ]
    )
    delta_r_earth = np.dot(A.T, r_ned).flatten()
    r_earth1 = r_earth0 + delta_r_earth
    lat, lon, depth = r_earth_2_geo(r_earth1)
    return np.array([lat, lon, depth])


def convert_sub_faults_geo2ned(sub_faults, source_point, approximate=True):
    """
    :param sub_faults: lat(deg), lon(deg), dep(km)
    :param source_point: lat(deg), lon(deg), dep(km)
    :param approximate:
    :return: sub_faults_ned (unit m)
    """
    sub_faults_c = sub_faults.copy()
    source_point_c = source_point.copy()
    sub_faults_c[:, 2] = sub_faults_c[:, 2] * 1e3
    source_point_c[2] = source_point_c[2] * 1e3
    sub_faults_ned = np.zeros_like(sub_faults_c)
    if not approximate:
        for n in range(sub_faults_c.shape[0]):
            sub_faults_ned[n, :] = convert_axis_delta_geo2ned(
                *source_point_c.copy(), *sub_faults_c[n, :].copy()
            ).flatten()
    else:
        for n in range(sub_faults_c.shape[0]):
            x = (sub_faults_c[n, 0] - source_point_c[0]) * d2m
            y = (sub_faults_c[n, 1] - source_point_c[1]) * d2m * np.cos(np.deg2rad(sub_faults_c[n, 0]))
            z = sub_faults_c[n, 2] - source_point_c[2]
            sub_faults_ned[n, :] = np.array([x, y, z])
    return sub_faults_ned


def select_df_geo(df: pd.DataFrame, lat_range, lon_range, time_range) -> pd.DataFrame:
    """
    :param df: must have columns "lat","lon","time"
    :param lat_range: [lat_min, lat_max]
    :param lon_range: [lon_min, lon_max]
    :param time_range: [time_min, time_max] in format "2000-01-01T00:00:00Z"
    return df_filtered
    """
    time_range[0] = pd.to_datetime(time_range[0])
    time_range[1] = pd.to_datetime(time_range[1])
    df["time"] = pd.to_datetime(df["time"])
    df_filtered = df[
        (df["lon"] >= lon_range[0])
        & (df["lon"] <= 180)
        & (df["lon"] >= -180)
        & (df["lon"] <= lon_range[1])
        & (df["lat"] >= lat_range[0])
        & (df["lat"] <= lat_range[1])
        & (df["time"] >= time_range[0])
        & (df["time"] <= time_range[1])
    ]
    return df_filtered


def geographic_centroid(points):
    """
    Calculate the approximate geographic center.

    :param points: NumPy array of shape (N, 2), where each row is [latitude_deg, longitude_deg]
    :return: (centroid_lat_deg, centroid_lon_deg)
    """
    lat_rad = np.deg2rad(points[:, 0])
    lon_rad = np.deg2rad(points[:, 1])

    x = np.cos(lat_rad) * np.cos(lon_rad)
    y = np.cos(lat_rad) * np.sin(lon_rad)
    z = np.sin(lat_rad)

    x_mean = x.mean()
    y_mean = y.mean()
    z_mean = z.mean()

    centroid_lat_rad = np.arctan2(z_mean, np.sqrt(x_mean**2 + y_mean**2))
    centroid_lon_rad = np.arctan2(y_mean, x_mean)

    centroid_lat_deg = np.rad2deg(centroid_lat_rad)
    centroid_lon_deg = np.rad2deg(centroid_lon_rad)

    return centroid_lat_deg, centroid_lon_deg


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


def create_rotate_z_mat(gamma):
    """
    Generates a rotation matrix about the Z-axis.
    From y to x.
    Parameters:
        gamma : float
            Rotation angle in radians.

    Returns:
        R : numpy.ndarray
            A 3x3 rotation matrix.
    """
    R = np.array(
        [
            [np.cos(gamma), -np.sin(gamma), 0],
            [np.sin(gamma), np.cos(gamma), 0],
            [0, 0, 1],
        ]
    )
    return R


def rotate_symmetric_tensor_series(tensor, gamma):
    """
    Rotates a series of symmetric tensors without using an explicit loop.

    Parameters:
        tensor: numpy array of shape (n, 6)
            Each row is [xx, xy, xz, yy, yz, zz] representing a symmetric tensor.
        gamma: float
            Rotation angle around the z-axis (in radians) used to create the rotation matrix.

    Returns:
        rotated_tensor: numpy array of shape (n, 6)
            Rotated tensor components in the same order as the input.
    """
    # Create the 3x3 rotation matrix (assumed to be defined elsewhere).
    R = create_rotate_z_mat(gamma)
    n = tensor.shape[0]

    # Construct full symmetric matrices from the condensed tensor representation.
    A = np.empty((n, 3, 3), dtype=tensor.dtype)
    A[:, 0, 0] = tensor[:, 0]
    A[:, 0, 1] = tensor[:, 1]
    A[:, 0, 2] = tensor[:, 2]
    A[:, 1, 0] = tensor[:, 1]
    A[:, 1, 1] = tensor[:, 3]
    A[:, 1, 2] = tensor[:, 4]
    A[:, 2, 0] = tensor[:, 2]
    A[:, 2, 1] = tensor[:, 4]
    A[:, 2, 2] = tensor[:, 5]

    # Rotate each tensor using batch matrix multiplication:
    # Compute rotated_A = R.T @ A @ R for each tensor.
    rotated_A = np.einsum("ij,njk,kl->nil", R.T, A, R)

    # Extract the independent components from the rotated tensors.
    rotated_tensor = np.empty((n, 6), dtype=tensor.dtype)
    rotated_tensor[:, 0] = rotated_A[:, 0, 0]
    rotated_tensor[:, 1] = rotated_A[:, 0, 1]
    rotated_tensor[:, 2] = rotated_A[:, 0, 2]
    rotated_tensor[:, 3] = rotated_A[:, 1, 1]
    rotated_tensor[:, 4] = rotated_A[:, 1, 2]
    rotated_tensor[:, 5] = rotated_A[:, 2, 2]

    return rotated_tensor
