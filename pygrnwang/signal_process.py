import numpy as np
from scipy import signal


def taper(data, taper_length=None, max_percentage=0.05) -> np.ndarray:
    data = data.copy()
    if taper_length is None:
        taper_length = max(2, round(len(data) * max_percentage))
    taper_window = signal.windows.hann(2 * taper_length)
    data[:taper_length] = data[:taper_length] * taper_window[:taper_length]
    data[-taper_length:] = data[-taper_length:] * taper_window[-taper_length:]
    return data


def cal_sos(srate, freq_band, butter_order=4):
    fn = srate / 2
    if (freq_band[0] == 0) and (freq_band[1] != 0) and (freq_band[1] / fn < 1):
        sos = signal.butter(
            butter_order, freq_band[1] / fn, btype="lowpass", output="sos"
        )
    elif (freq_band[0] != 0) and ((freq_band[1] == 0) or (freq_band[1] / fn >= 1)):
        sos = signal.butter(
            butter_order, freq_band[0] / fn, btype="highpass", output="sos"
        )
    elif (freq_band[0] != 0) and (freq_band[1] != 0) and (freq_band[1] / fn < 1):
        sos = signal.butter(
            butter_order,
            [freq_band[0] / fn, freq_band[1] / fn],
            btype="bandpass",
            output="sos",
        )
    else:
        sos = None
    return sos


def filter_butter(data: np.ndarray, srate, freq_band, butter_order=4, zero_phase=False):
    data = data.copy()
    sos = cal_sos(srate, freq_band, butter_order)
    if sos is not None:
        if zero_phase:
            data = signal.sosfiltfilt(sos, data)
        else:
            data = signal.sosfilt(sos, data)
    return data


def resample(data, srate_old: float, srate_new: float, zero_phase=True):
    if zero_phase:
        if float(srate_old).is_integer() and float(srate_new).is_integer():
            srate_old = int(srate_old)
            srate_new = int(srate_new)
            gcd = np.gcd(srate_new, srate_old)
            p = srate_new // gcd
            q = srate_old // gcd
            data = signal.resample_poly(data, p, q)
        else:
            data = signal.resample(data, round(len(data) * srate_new / srate_old))
            data = filter_butter(
                data=data,
                srate=srate_new,
                freq_band=[0, srate_old / 2],
                zero_phase=zero_phase,
            )
    else:
        q = srate_old / srate_new
        if srate_new < srate_old:
            if q.is_integer():
                data = signal.decimate(
                    data, q=int(q), ftype="fir", zero_phase=zero_phase
                )
            else:
                data = filter_butter(
                    data=data,
                    srate=srate_old,
                    freq_band=[0, srate_new / 2],
                    zero_phase=zero_phase,
                )
                data = signal.resample(x=data, num=round(len(data) * q))
        elif srate_new > srate_old:
            data = signal.resample(data, round(len(data) * q))
            data = filter_butter(
                data=data,
                srate=srate_new,
                freq_band=[0, srate_old / 2],
                zero_phase=zero_phase,
            )
    return data


def linear_interp(data, N_new) -> np.ndarray:
    points_loc = np.arange(0, len(data))
    points_loc_new = np.linspace(0, len(data), N_new, endpoint=False)
    data_new = np.interp(points_loc_new, points_loc, data)
    return data_new
