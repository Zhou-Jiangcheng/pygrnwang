import numpy as np
from scipy import signal
from scipy.signal import find_peaks
from scipy.interpolate import interp1d


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


def filter_butter(data: np.ndarray, srate, freq_band, butter_order=4, zerophase=False):
    data = data.copy()
    sos = cal_sos(srate, freq_band, butter_order)
    if sos is not None:
        if zerophase:
            data = signal.sosfiltfilt(sos, data)
        else:
            data = signal.sosfilt(sos, data)
    else:
        pass
    return data


def cal_factors(N):
    factors = []
    i = 2
    M = N
    while i < N:
        if M % i == 0:
            factors.append(i)
            M = M / i
        else:
            i = i + 1
    return factors


def resample(data: np.ndarray, srate_old: float, srate_new: float, zero_phase=False):
    data = data.copy()
    if srate_new < srate_old:
        q = srate_old / srate_new
        if q.is_integer():
            q = int(q)
            if q > 10:
                factors = cal_factors(q)
                # print(factors)
                for factor in factors:
                    data = signal.decimate(data, q=factor, zero_phase=zero_phase)
            else:
                data = signal.decimate(data, q=q, zero_phase=zero_phase)
        else:
            data = filter_butter(
                data=data,
                srate=srate_old,
                freq_band=[0, srate_new / 2],
                zerophase=zero_phase,
            )
            data = signal.resample(x=data, num=round(len(data) * srate_new / srate_old))
    elif srate_new > srate_old:
        data = signal.resample(data.copy(), round(len(data) * srate_new / srate_old))
        data = filter_butter(
            data=data,
            srate=srate_new,
            freq_band=[0, srate_old / 2],
            zerophase=zero_phase,
        )
    else:
        pass
    return data


def linear_interp(data, N_new) -> np.ndarray:
    points_loc = np.arange(0, len(data))
    points_loc_new = np.linspace(0, len(data), N_new, endpoint=False)
    data_new = np.interp(points_loc_new, points_loc, data)
    return data_new


def smooth(data: np.ndarray, smooth_window: int):
    """
    it should be noted the length of data after smooth
    will be different with origin data
    :param data:
    :param smooth_window:
    :return:
    """
    if smooth_window is None:
        smooth_window = max(3, round(0.02 * len(data)))
    if smooth_window == 0:
        return data
    data = np.concatenate(
        [
            np.ones(smooth_window) * data[0],
            data,
            np.ones(smooth_window) * data[-1],
        ]
    )
    data_smooth = np.convolve(
        data, np.ones(smooth_window) / smooth_window, mode="same"
    )[smooth_window:-smooth_window]
    # data_smooth = np.convolve(data, np.ones(smooth_window) / smooth_window, mode="full")
    return data_smooth


def smooth_cyclic(data: np.ndarray, smooth_window: int):
    if smooth_window is None:
        smooth_window = max(3, round(0.02 * len(data)))
    if smooth_window == 0:
        return data
    data_smooth = np.convolve(data, np.ones(smooth_window) / smooth_window, mode="same")
    return data_smooth


def smooth_except_both_ends(v, n):
    v_temp = np.zeros_like(v)
    for i in range(n):
        v_temp += np.roll(v, i - n // 2, axis=0)
    v_temp /= n
    v_new = np.zeros_like(v)
    v_new[: n // 2] = v[: n // 2]
    v_new[n // 2 : -n // 2] = v_temp[n // 2 : -n // 2]
    v_new[-n // 2 :] = v[-n // 2 :]
    return v_new


def cal_sliding_average(
    time: np.ndarray,
    data: np.ndarray,
    window_length: int,
    overlap: int = 0,
) -> np.ndarray:
    """

    :param time: s
    :param data:
    :param window_length: s
    :param overlap: s
    :return:
    """
    data = data.copy()
    data_new = []
    n = 0
    inds = [0]
    for i in range(len(time)):
        if time[i] < (n + 1) * window_length:
            pass
        else:
            inds.append(i)
            data_new.append(np.mean(data[inds[n] : inds[n + 1]]))
            n = n + 1
    data_new = np.array(data_new)
    return data_new


def Gauss(f, a):
    # 高斯滤波器,a=2.5 fc=1.2, a=1.25 fc=0.6, a=0.625 fc=0.3
    a = np.ones_like(f) * a
    return np.e ** (-((2 * np.pi * f) ** 2) / (4 * a**2))


def add_noise(data: np.ndarray, snr: float):
    """
    :param data:
    :param snr: unit None
    :return: slips
    """
    data = data.copy()
    noise = np.random.normal(0, 1, (len(data),))
    P_signal = np.sum(data**2)
    P_noise = P_signal / snr
    noise = np.sqrt(P_noise) * noise / np.sqrt(np.sum(noise**2))
    data = data + noise
    return data


def add_noise_db(data: np.ndarray, snr_db: float):
    """
    :param data:
    :param snr_db: unit db
    :return: slips
    """
    data = data.copy()
    noise = np.random.normal(0, 1, (len(data),))
    P_signal = np.sum(data**2)
    P_noise = P_signal / (10 ** (snr_db / 10))
    noise = np.sqrt(P_noise) * noise / np.sqrt(np.sum(noise**2))
    data = data + noise
    return data


def intype2int(intype):
    if intype == "DISP":
        intype = 0
    elif intype == "VEL":
        intype = 1
    elif intype == "ACC":
        intype = 2
    else:
        raise ValueError("type wrong")
    return intype


def convert_by_type(input_type: str, output_type: str, data: np.ndarray, srate: float):
    data = data.copy()

    input_type = intype2int(input_type)
    output_type = intype2int(output_type)
    delta = input_type - output_type
    if delta > 0:
        for i in range(delta):
            data = np.cumsum(data) * 1 / srate
    elif delta < 0:
        for i in range(-delta):
            data = np.gradient(data) * srate
    else:
        pass
    return data


def find_not_zeros_at_end(w: np.ndarray):
    """
    :param w:
    :return: ind
    """
    w = w.copy()[::-1]
    end_point = np.nonzero(w)[0][0]
    end_point = len(w) - end_point
    return end_point


def merge_spectrum(f_band_1, spec1, f_band_2, spec2, delta_f, freq_gap):
    spec1 = linear_interp(spec1, round((f_band_1[-1] - f_band_1[0]) / delta_f))
    spec2 = linear_interp(spec2, round((f_band_2[-1] - f_band_2[0]) / delta_f))
    N_freq_gap = round(freq_gap / delta_f)
    weight = np.linspace(0, 1, 2 * N_freq_gap)
    spec_merge = np.zeros(len(spec1) + len(spec2) - 2 * N_freq_gap)
    spec_merge[: len(spec1) - 2 * N_freq_gap] = spec1[: len(spec1) - 2 * N_freq_gap]
    spec_merge[len(spec1) - 2 * N_freq_gap : len(spec1)] = (
        spec1[len(spec1) - 2 * N_freq_gap :] * (1 - weight)
        + spec2[: 2 * N_freq_gap] * weight
    )
    spec_merge[len(spec1) :] = spec2[2 * N_freq_gap :]
    f_merge = np.linspace(
        min(f_band_1[0], f_band_2[0]), max(f_band_1[-1], f_band_2[-1]), len(spec_merge)
    )
    return f_merge, spec_merge


def upper_envelope(data, distance=None):
    if distance is not None:
        peaks, _ = find_peaks(data, distance=distance)
    else:
        peaks, _ = find_peaks(data)
    # print(peaks)
    # Extract the peak values and their corresponding indices
    peak_values = data[peaks]
    peak_indices = np.arange(len(data))[peaks]

    if distance is not None:
        peak_values = np.concatenate(
            [data[: peaks[0] - 1][::distance].flatten(), peak_values]
        )
        peak_indices = np.concatenate(
            [np.arange(peaks[0] - 1)[::distance], peak_indices]
        )
    else:
        peak_values = np.concatenate([data[: peaks[0] - 1].flatten(), peak_values])
        peak_indices = np.concatenate([np.arange(peaks[0] - 1), peak_indices])

    # Create an interpolation function for the peaks
    interp_func = interp1d(
        peak_indices, peak_values, kind="linear", fill_value="extrapolate"
    )

    # Generate the envelope values for the entire dataset
    envelope = interp_func(np.arange(len(data)))
    return envelope, peak_indices, peak_values


def time_shift_b(a: np.ndarray, b: np.ndarray):
    """
    Shift array b to the position where its correlation with a is maximized.
    If shifted backwards, prepend zeros; if shifted forwards, truncate the front.
    :return shifted_b, best_shift, max_correlation
    """
    a = a.copy()
    b = b.copy()
    max_correlation = -np.inf
    best_shift = 0

    for shift in range(-len(b) + 1, len(a)):
        # print("shift len(a), len(b)", shift, len(a), len(b))
        if shift < 0:
            shifted_b = b[-shift:]  # Truncate the front of b
            shifted_b = np.concatenate([shifted_b, np.zeros(-shift)])
        else:
            shifted_b = np.pad(b, (shift, 0), "constant", constant_values=0)[
                : len(a)
            ]  # Pad and trim b to match the length of a
        # print("shift len(a), len(b)", shift, len(a), len(b))
        if np.max(a) == 0 or np.max(shifted_b) == 0:
            correlation = 0
        else:
            correlation = np.corrcoef(a / np.max(a), shifted_b / np.max(shifted_b))[
                0, 1
            ]
        if correlation > max_correlation:
            max_correlation = correlation
            best_shift = shift

    # Apply the best shift to b
    if best_shift < 0:
        shifted_b = b[-best_shift:]
    else:
        shifted_b = np.pad(b, (best_shift, 0), "constant", constant_values=0)[: len(a)]

    return shifted_b, best_shift, max_correlation
