import os
import struct
import numpy as np
from matplotlib import pyplot as plt
from obspy.taup import TauPyModel
from obspy.geodetics import kilometers2degrees
from pygrnwang.utils import (
    cal_m0_from_mt,
    rotate_rtz_to_enz,
    find_nearest_dichotomy,
    plane2mt,
)
from pygrnwang.signal_process import taper, filter_butter, linear_interp, resample


class Spgrn(object):
    def __init__(
        self,
        path_green_lib,
        event_depth_in_km=None,
        az_in_deg=None,
        dist_in_km=None,
        focal_mechanism=None,
        before_tp=120,
        time_window=None,
        srate=None,
        freq_min=0,
        freq_max=0,
        filter_zerophase=False,
        taper_length=0,
        tension_compression_flag=False,
    ):
        """
        init the read class of spgrn
        :param path_green_lib: the path of green's function lib,
        the "GreenFunc" and "GreenSpec" dirs should in this dir,
        and the grn_d file and other files are stored in these 2 dirs
        by recv_depth in sub dirs 1.0, 2.0, 3.0 ... separately
        if the green lib is created by create_spgrn_bulk.py
        :param event_depth_in_km:
        :param az_in_deg: from n to e (clockwise)
        :param dist_in_km: distance in km
        :param focal_mechanism: list, [M0, strike, dip, rake]/
                                      [M11,M12,M13,M22,M23,M33]/
                                      [M0,M11,M12,M13,M22,M23,M33]
        :param before_tp: the first count in green lib before tP,
                          must be the minus of the value 'before_p'
                          set in create_spgrn.py/create_spgrn_bulk.py
                          in seconds
        :param time_window: the length of seismograms, in seconds.
        the seismograms will be cut from begin_time to begin_time+time_window
        :param srate: in Hz
        :param freq_min: Hz
        :param freq_max: Hz
        :param filter_zerophase: True/False
        :param taper_length: length of taper at both ends in counts,
                             5 percentage if set to None
        :param tension_compression_flag: tension or compression the waveform
        to make it match the arrival time of first P and S wave
        """
        if focal_mechanism is None:
            focal_mechanism = [1, 1, 0, 0, 1, 0, 1]
        self.path_green = os.path.join(path_green_lib, "GreenFunc")
        self.event_depth = event_depth_in_km
        self.az_in_deg = az_in_deg
        self.dist_in_km = dist_in_km
        self.focal_mechanism = focal_mechanism
        self.before_tP = before_tp  # s
        self.time_window = time_window  # s

        self.freq_min = freq_min
        self.freq_max = freq_max
        self.filter_zerophase = filter_zerophase
        self.taper_length = taper_length
        self.tension_compression_flag = tension_compression_flag

        self.dist_in_deg = kilometers2degrees(dist_in_km)
        self.green_depth_list = self.find_green_depth_list()
        self.green_depth = find_nearest_dichotomy(
            self.event_depth, self.green_depth_list
        )[0]
        self.green_info = self.read_green_info()
        # {
        #     "time_window": time_window,
        #     "sampling_interval": sampling_interval,
        #     "samples_num": int(samples_num),
        #     "dist_list": dist_list
        # }
        self.green_tpts_table = self.read_tpts_table()
        # {
        #     "p_onset": tp_onset,
        #     "p_takeoff": tp_takeoff,
        #     "p_slowness": tp_slowness,
        #     "s_onset": ts_onset,
        #     "s_takeoff": ts_takeoff,
        #     "s_slowness": ts_slowness
        #
        # }
        if srate is not None:
            self.srate = srate
        else:
            self.srate = 1 / self.green_info["sampling_interval"]
        self.waveforms = None

    @staticmethod
    def get_number_in_line(line):
        # 获取这一行greeninfo.dat里的数字
        numbers = []
        for item in line.split():
            if item != "":
                numbers.append(float(item))
        return numbers

    def find_green_depth_list(self):
        # 获得格林函数库中已计算的深度
        green_depth_list: list = os.listdir(self.path_green)
        for i in range(len(green_depth_list)):
            green_depth_list[i] = float(green_depth_list[i])
        green_depth_list.sort()
        return green_depth_list

    def read_green_info(self):
        # 读取格林函数库说明文件
        with open(
            os.path.join(
                self.path_green,
                "%.1f" % self.green_depth,
                "GreenInfo%.1f.dat" % self.green_depth,
            ),
            "r",
        ) as fr:
            lines = fr.readlines()
        [time_window, sampling_interval, samples_num] = self.get_number_in_line(
            lines[6]
        )
        number_of_distance = self.get_number_in_line(lines[9])[0]
        dist_list = []
        for i in range(int(np.ceil(number_of_distance / 5))):
            temp = self.get_number_in_line(lines[10 + i])
            for item in temp:
                dist_list.append(item)
        return {
            "time_window": time_window,
            "sampling_interval": sampling_interval,
            "samples_num": int(samples_num),
            "dist_list": dist_list,
        }

    def read_tpts_table(self):
        #   {green_function_distance,
        #      {onset [s], takeoff [deg], slowness [s/m]}}
        #   读取格林函数库的理论到时
        dist_green_num = find_nearest_dichotomy(
            self.dist_in_km, self.green_info["dist_list"]
        )[1]
        start_count = 4 + dist_green_num * 12

        def seek_time_table(start_count_in, phase):
            if phase == "P":
                fname = "tptable.dat"
            elif phase == "S":
                fname = "tstable.dat"
            else:
                raise ValueError("phase must be P or S")
            fr = open(
                os.path.join(self.path_green, "%.1f" % self.green_depth, fname), "rb"
            )
            fr.seek(start_count_in)
            t_onset = struct.unpack("f", fr.read(4))[0]
            t_takeoff = struct.unpack("f", fr.read(4))[0]
            t_slowness = struct.unpack("f", fr.read(4))[0]
            fr.close()
            return [t_onset, t_takeoff, t_slowness]

        [tp_onset, tp_takeoff, tp_slowness] = seek_time_table(start_count, "P")
        [ts_onset, ts_takeoff, ts_slowness] = seek_time_table(start_count, "S")

        return {
            "p_onset": tp_onset,
            "p_takeoff": tp_takeoff,
            "p_slowness": tp_slowness,
            "s_onset": ts_onset,
            "s_takeoff": ts_takeoff,
            "s_slowness": ts_slowness,
        }

    def read_time_series(self):
        dist_green_num = find_nearest_dichotomy(
            self.dist_in_km, self.green_info["dist_list"]
        )[1]
        length_each = 3 + (1 + self.green_info["samples_num"] + 1) * 10
        start_count = dist_green_num * length_each

        fr = open(
            os.path.join(
                self.path_green,
                "%.1f" % self.green_depth,
                "grn_d%.1f" % self.green_depth,
            ),
            "rb",
        )
        time_series = np.fromfile(
            file=fr, dtype=np.float32, count=length_each, sep="", offset=start_count * 4
        )
        time_series = time_series[3:].reshape(10, self.green_info["samples_num"] + 2)
        time_series = time_series[:, 1:-1]
        fr.close()
        return time_series

    def synthesize(self):
        # 根据震源机制合成地震图
        t = self.read_time_series()
        if len(self.focal_mechanism) == 3:
            mt = plane2mt(
                1,
                self.focal_mechanism[0],
                self.focal_mechanism[1],
                self.focal_mechanism[2],
            )
            [M11, M12, M13, M22, M23, M33] = list(mt)
        elif len(self.focal_mechanism) == 4:
            mt = plane2mt(
                self.focal_mechanism[0],
                self.focal_mechanism[1],
                self.focal_mechanism[2],
                self.focal_mechanism[3],
            )
            [M11, M12, M13, M22, M23, M33] = list(mt)
        elif len(self.focal_mechanism) == 6:
            [M11, M12, M13, M22, M23, M33] = self.focal_mechanism
        elif len(self.focal_mechanism) == 7:
            M0 = self.focal_mechanism[0]
            temp = np.array(self.focal_mechanism[1:])
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
        # print([M11, M12, M13, M22, M23, M33])
        # expl.,strike-slip,dip-slip,clvd
        exp = (M11 + M22 + M33) / 3
        clvd = M33 - exp
        ss12 = -M12
        ss11 = (M11 - M22) / 2
        ds31 = M13
        ds23 = -M23

        az = (180 - self.az_in_deg) * np.pi / 180
        sin_az, cos_az = np.sin(az), np.cos(az)
        sin_2az, cos_2az = np.sin(2 * az), np.cos(2 * az)
        m1 = [
            exp,
            clvd,
            (ss12 * sin_2az) + ss11 * cos_2az,
            (ds31 * cos_az) + ds23 * sin_az,
        ]
        m2 = [(ss12 * cos_2az - ss11 * sin_2az), (ds31 * sin_az - ds23 * cos_az)]
        z = t[0] * m1[0] + t[8] * m1[1] + t[2] * m1[2] + t[5] * m1[3]
        r = t[1] * m1[0] + t[9] * m1[1] + t[3] * m1[2] + t[6] * m1[3]
        t = t[4] * m2[0] + t[7] * m2[1]

        seismograms = [r, t, z]
        return seismograms

    def cal_traveltime(self, model_name="ak135f_no_mud"):
        """
        cal the traveltime of first P and first S
        :param model_name: default ak135f_no_mud
        :return: [first P, first S] in seconds
        """
        dist_in_deg = kilometers2degrees(self.dist_in_km)
        model = TauPyModel(model_name)
        P_group = model.get_travel_times(
            source_depth_in_km=self.event_depth,
            distance_in_degree=dist_in_deg,
            phase_list=["P"],
        )
        if len(P_group) != 0:
            first_P = P_group[0].time
        else:
            P_group = model.get_travel_times(
                source_depth_in_km=self.event_depth,
                distance_in_degree=dist_in_deg,
                phase_list=["Pdiff"],
            )
            if len(P_group) != 0:
                first_P = P_group[0].time
            else:
                P_group = model.get_travel_times(
                    source_depth_in_km=self.event_depth,
                    distance_in_degree=dist_in_deg,
                    phase_list=["PKP", "PKIKP", "PKiKP"],
                )
                temp = np.inf
                for i in range(len(P_group)):
                    if P_group[i].time < temp:
                        temp = P_group[i].time
                first_P = temp
        S_group = model.get_travel_times(
            source_depth_in_km=self.event_depth,
            distance_in_degree=dist_in_deg,
            phase_list=["S"],
        )
        if len(S_group) != 0:
            first_S = S_group[0].time
        else:
            S_group = model.get_travel_times(
                source_depth_in_km=self.event_depth,
                distance_in_degree=dist_in_deg,
                phase_list=["SKS", "SKIKS", "SKiKS"],
            )
            temp = np.inf
            for i in range(len(S_group)):
                if S_group[i].time < temp:
                    temp = S_group[i].time
            first_S = temp
        return [first_P, first_S]

    def cut(self, seismograms):
        begin_count = 0
        end_count = round((0 + self.time_window) * self.srate)
        for i in range(3):
            seismograms[i] = seismograms[i][begin_count:end_count]
        return seismograms

    def shift(self, seismograms):
        # 将计算好的一定震中距的理论格林函数转换为指定震中距的
        [first_P, first_S] = self.cal_traveltime()
        sampling_interval = self.green_info["sampling_interval"]
        first_P = round(first_P / sampling_interval)
        first_S = round(first_S / sampling_interval)
        first_P_green = round(self.green_tpts_table["p_onset"] / sampling_interval)
        first_S_green = round(self.green_tpts_table["s_onset"] / sampling_interval)
        for i in range(3):
            tP_in_count = round(self.before_tP * self.srate)
            before_P = seismograms[i][:tP_in_count]
            P_S = linear_interp(
                seismograms[i][
                    tP_in_count : tP_in_count + first_S_green - first_P_green
                ],
                first_S - first_P,
            )
            after_S = seismograms[i][tP_in_count + first_S_green - first_P_green :]
            seismograms[i] = np.concatenate([before_P, P_S, after_S])
        return seismograms

    def read_rtz(self):
        seismograms = self.synthesize()
        if self.filter_zerophase:
            for i in range(len(seismograms)):
                seismograms[i] = taper(seismograms[i])
                seismograms[i] = filter_butter(
                    seismograms[i],
                    round(1 / self.green_info["sampling_interval"]),
                    [self.freq_min, self.freq_max],
                    zerophase=True,
                )
        else:
            for i in range(len(seismograms)):
                seismograms[i] = filter_butter(
                    seismograms[i],
                    round(1 / self.green_info["sampling_interval"]),
                    [self.freq_min, self.freq_max],
                )
        for i in range(3):
            seismograms[i] = resample(
                seismograms[i],
                srate_old=round(1 / self.green_info["sampling_interval"]),
                srate_new=self.srate,
            )
        seismograms = self.cut(seismograms)
        if self.tension_compression_flag:
            seismograms = self.shift(seismograms)

        self.waveforms = seismograms
        return seismograms

    def read_enz(self):
        [r, t, z] = self.read_rtz()
        [e, n, z] = rotate_rtz_to_enz(self.az_in_deg, r, t, z)
        self.waveforms = [e, n, z]
        return [e, n, z]

    def read_rtz_raw(self):
        """
        read the raw seismograms, don't do cut, filter, resample etc.
        srate is the same as the green function
        time_window is the same as the green function
        :return: [r, t, z]
        """
        seismograms = self.synthesize()
        if self.tension_compression_flag:
            seismograms = self.shift(seismograms)
        self.waveforms = seismograms
        return seismograms

    def read_enz_raw(self):
        """
        read the raw seismograms, don't do cut, filter, resample etc.
        srate is the same as the green's function
        time_window is the same as the green function
        :return: [e, n, z]
        """
        [r, t, z] = self.read_rtz_raw()
        [e, n, z] = rotate_rtz_to_enz(self.az_in_deg, r, t, z)
        self.waveforms = [e, n, z]
        return [e, n, z]

    def plot_waveforms(self, begin_time=0, length_in_s=None):
        if length_in_s is None:
            length = len(self.waveforms[0][:])
        else:
            length = length_in_s * self.srate
        b = round(begin_time * self.srate)
        t = np.arange(0, round((length - b) / self.srate), 1 / self.srate)
        plt.figure()
        plt.subplot(3, 1, 1)
        plt.title(
            "recv_depth:%.2f az:%.2f dist:%.2f"
            % (self.green_depth, self.az_in_deg, self.dist_in_km)
        )
        plt.plot(t, self.waveforms[0][b:length], color="red", linewidth=0.2)
        plt.subplot(3, 1, 2)
        plt.plot(t, self.waveforms[1][b:length], color="green", linewidth=0.2)
        plt.subplot(3, 1, 3)
        plt.plot(t, self.waveforms[2][b:length], color="blue", linewidth=0.2)
        plt.xlabel("t [s]")
        plt.ylabel("u [m]")
        plt.show()


if __name__ == "__main__":
    # print(gps2dist_azimuth(0, 0, -10, 10))
    path = "/e/spgrn2020_lib"
    az_ = 30
    dep_ = 50
    dist_in_km_ = 8000
    fm_ = [1e19] + [30, 40, 50]
    srate_ = 2
    print(dist_in_km_)
    seismograms1 = Spgrn(
        path_green_lib=path,
        event_depth_in_km=50,
        az_in_deg=az_,
        dist_in_km=dist_in_km_,
        focal_mechanism=[1e19] + [30, 40, 50],
        time_window=4096,
    )
    seismograms1.read_enz()
    seismograms1.plot_waveforms()
