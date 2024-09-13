import os.path
import numpy as np
import pandas as pd
from pygrnwang.create_qseis import convert2bin
from pygrnwang.read_qseis import read_time_series, synthesize


def read_qseis06_grn_raw(path_qseis_data, az, fm, num):
    ex_z_df = pd.read_csv(os.path.join(path_qseis_data, "ex.tz"), sep="\\s+")
    ex_r_df = pd.read_csv(os.path.join(path_qseis_data, "ex.tr"), sep="\\s+")
    ss_z_df = pd.read_csv(os.path.join(path_qseis_data, "ss.tz"), sep="\\s+")
    ss_r_df = pd.read_csv(os.path.join(path_qseis_data, "ss.tr"), sep="\\s+")
    ss_t_df = pd.read_csv(os.path.join(path_qseis_data, "ss.tt"), sep="\\s+")
    ds_z_df = pd.read_csv(os.path.join(path_qseis_data, "ds.tz"), sep="\\s+")
    ds_r_df = pd.read_csv(os.path.join(path_qseis_data, "ds.tr"), sep="\\s+")
    ds_t_df = pd.read_csv(os.path.join(path_qseis_data, "ds.tt"), sep="\\s+")
    cl_z_df = pd.read_csv(os.path.join(path_qseis_data, "cl.tz"), sep="\\s+")
    cl_r_df = pd.read_csv(os.path.join(path_qseis_data, "cl.tr"), sep="\\s+")
    time_series = np.concatenate(
        [
            ex_z_df.iloc[:, num].values,
            ex_r_df.iloc[:, num].values,
            ss_z_df.iloc[:, num].values,
            ss_r_df.iloc[:, num].values,
            ss_t_df.iloc[:, num].values,
            ds_z_df.iloc[:, num].values,
            ds_r_df.iloc[:, num].values,
            ds_t_df.iloc[:, num].values,
            cl_z_df.iloc[:, num].values,
            cl_r_df.iloc[:, num].values,
        ]
    ).T
    time_series = time_series.reshape(10, len(ex_z_df))
    # z,t,r
    seismograms = synthesize(az_in_deg=az, time_series=time_series, focal_mechanism=fm)
    # r = ss_r_df.iloc[:, num].values
    # t = ss_t_df.iloc[:, num].values
    # z = ss_z_df.iloc[:, num].values
    return seismograms


def read_qseis06_grn_bin(path_qseis_data, az, fm, num, sampling_num):
    time_series = read_time_series(
        path_greenfunc=path_qseis_data,
        start_count=num * sampling_num * 10,
        sampling_num=sampling_num,
    )
    # z,t,r
    seismograms = synthesize(az_in_deg=az, time_series=time_series, focal_mechanism=fm)
    return seismograms


def read_qseis06_output(path_qseis_data, num):
    df_r = pd.read_csv(os.path.join(path_qseis_data, "out.tr"), sep="\\s+")
    df_t = pd.read_csv(os.path.join(path_qseis_data, "out.tt"), sep="\\s+")
    df_z = -pd.read_csv(os.path.join(path_qseis_data, "out.tz"), sep="\\s+")
    # print(df_r.iloc[:, num].values.shape)
    return df_z.iloc[:, num].values, df_t.iloc[:, num].values, df_r.iloc[:, num].values


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    path_output0 = "/e/all_green_lib/qseis_lib_test/test/"
    r_0, t_0, z_0 = read_qseis06_grn_raw(path_output0, 60.0, [243.0, 86.0, 7.0], 1)
    path_output1 = "/e/all_green_lib/qseis_lib_test/test_alias/"
    r_1, t_1, z_1 = read_qseis06_grn_raw(path_output1, 60.0, [243.0, 86.0, 7.0], 1)

    plt.figure()
    plt.plot(r_0)
    plt.plot(r_1)
    plt.show()

    plt.figure()
    plt.plot(t_0)
    plt.plot(t_1)
    plt.show()

    plt.figure()
    plt.plot(z_0)
    plt.plot(z_1)
    plt.show()
