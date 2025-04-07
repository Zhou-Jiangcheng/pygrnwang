import os
import shutil

import numpy as np
import pandas as pd

d2m = 111194.92664455874


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

    # 对向量进行矩阵乘法
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

    # 对向量进行矩阵乘法
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


def cal_travel_time_taup(
    event_depth_km, dist_km, phase, receiver_depth_km=0, model_name="ak135"
):
    cmd = "taup time -ph %s -model %s -h %f --stadepth %f -km %f --time" % (
        phase,
        model_name,
        event_depth_km,
        receiver_depth_km,
        dist_km,
    )
    result = os.popen(cmd).read()
    t = float(result.strip())
    return t


def cal_first_p_s(event_depth_km, dist_km, receiver_depth_km=0, model_name="ak135"):
    if event_depth_km < receiver_depth_km:
        event_depth_km, receiver_depth_km = receiver_depth_km, event_depth_km
    if shutil.which("taup"):
        cmd = "taup time -ph %s -model %s -h %f --stadepth %f -km %f --time" % (
            "p,P,Pn,PKP,Pdiff",
            model_name,
            event_depth_km,
            receiver_depth_km,
            dist_km,
        )
        result = os.popen(cmd).read()
        first_P = min(list(map(float, result.split())))

        cmd = "taup time -ph %s -model %s -h %f --stadepth %f -km %f --time" % (
            "s,S,Sn,SKS",
            model_name,
            event_depth_km,
            receiver_depth_km,
            dist_km,
        )
        result = os.popen(cmd).read()
        first_S = min(list(map(float, result.split())))
        return first_P, first_S
    else:
        raise ValueError(
            "Taup must be installed and placed in the PATH environment variable!"
        )


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
    :param sub_faults:
    :param source_point:
    :param approximate:
    :return: sub_faults_ned (unit m)
    """
    origin_point = source_point.copy()
    sub_faults_ned = np.zeros_like(sub_faults)
    if not approximate:
        for n in range(sub_faults.shape[0]):
            sub_faults_ned[n, :] = convert_axis_delta_geo2ned(
                *origin_point, *sub_faults[n, :]
            ).flatten()
    else:
        for n in range(sub_faults.shape[0]):
            x = (sub_faults[n, 0] - origin_point[0]) * d2m
            y = (sub_faults[n, 1] - origin_point[1]) * d2m
            z = sub_faults[n, 2] - origin_point[2]
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


if __name__ == "__main__":
    pass
