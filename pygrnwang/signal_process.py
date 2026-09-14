import numpy as np
from scipy import signal


def taper(data, taper_length=None, max_percentage=0.05) -> np.ndarray:
    """Apply a Hann taper to both ends of a one-dimensional trace.

    Parameters
    ----------
    data : numpy.ndarray
        One-dimensional uniformly sampled signal; input is not modified.
    taper_length : int or None, optional
        Number of samples tapered at each end. None uses max(2, round(len(data)*max_percentage)); keep within the signal length. Default: None.
    max_percentage : float, optional
        Fraction of trace length used at each tapered end when taper_length is None. Default: 0.05.

    Returns
    -------
    tapered : numpy.ndarray
        Same shape and units as the input.

    Raises
    ------
    ValueError
        The taper does not fit the trace.

    Notes
    -----
    The function copies the input. For very short traces, supply a compatible explicit taper_length.
    """
    data = data.copy()
    if taper_length is None:
        taper_length = max(2, round(len(data) * max_percentage))
    taper_window = signal.windows.hann(2 * taper_length)
    data[:taper_length] = data[:taper_length] * taper_window[:taper_length]
    data[-taper_length:] = data[-taper_length:] * taper_window[-taper_length:]
    return data


def cal_sos(srate, freq_band, butter_order=4):
    """Design lowpass, highpass or bandpass Butterworth second-order sections.

    Parameters
    ----------
    srate : float
        Positive output sampling rate in Hz.
    freq_band : sequence of float or None
        Two cutoff frequencies [low, high] in Hz. None disables filtering in readers; a missing corner selects lowpass or highpass.
    butter_order : int, optional
        Butterworth filter order. Default: 4.

    Returns
    -------
    sos : numpy.ndarray or None
        Shape (number_of_sections, 6); None means that no filtering is needed.

    Raises
    ------
    ValueError
        Requested filter order or nonzero cutoff frequencies are invalid.

    Notes
    -----
    Supply a two-element freq_band. None or zero at a corner means no cutoff there. A high corner at or above Nyquist is ignored; a nonzero low corner then selects highpass.
    """
    fn = srate / 2
    low = 0 if freq_band[0] is None else freq_band[0]
    high = 0 if freq_band[1] is None else freq_band[1]
    if (low == 0) and (high != 0) and (high / fn < 1):
        sos = signal.butter(butter_order, high / fn, btype="lowpass", output="sos")
    elif (low != 0) and ((high == 0) or (high / fn >= 1)):
        sos = signal.butter(butter_order, low / fn, btype="highpass", output="sos")
    elif (low != 0) and (high != 0) and (high / fn < 1):
        sos = signal.butter(
            butter_order,
            [low / fn, high / fn],
            btype="bandpass",
            output="sos",
        )
    else:
        sos = None
    return sos


def filter_butter(data: np.ndarray, srate, freq_band, butter_order=4, zero_phase=False):
    """Filter an array along its last axis with a Butterworth filter.

    Parameters
    ----------
    data : numpy.ndarray
        Signal array; filtering acts along the final axis. The input is copied.
    srate : float
        Positive output sampling rate in Hz.
    freq_band : sequence of float or None
        Two cutoff frequencies [low, high] in Hz. None disables filtering in readers; a missing corner selects lowpass or highpass.
    butter_order : int, optional
        Butterworth filter order. Default: 4.
    zero_phase : bool, optional
        True applies forward/backward filtering; False uses causal filtering. Default: False.

    Returns
    -------
    filtered : numpy.ndarray
        Same shape and units as data; a copy is returned even when filtering is disabled.

    Raises
    ------
    ValueError
        Filter parameters are invalid or the trace is too short for zero-phase padding.

    Notes
    -----
    Unlike the reader wrappers, this function requires a two-element freq_band; use [None, None] to disable it. Forward/backward filtering requires enough samples for padding.
    """
    data = data.copy()
    sos = cal_sos(srate, freq_band, butter_order)
    if sos is not None:
        if zero_phase:
            data = signal.sosfiltfilt(sos, data)
        else:
            data = signal.sosfilt(sos, data)
    return data


def resample(data, srate_old: float, srate_new: float, zero_phase=True):
    """Resample a one-dimensional signal with rate-dependent antialias handling.

    Parameters
    ----------
    data : numpy.ndarray
        One-dimensional uniformly sampled signal; input is not modified.
    srate_old : float
        Original positive sampling rate in Hz.
    srate_new : float
        Desired positive sampling rate in Hz.
    zero_phase : bool, optional
        True applies forward/backward filtering; False uses causal filtering. Default: True.

    Returns
    -------
    resampled : numpy.ndarray
        One-dimensional output at srate_new; length is determined by the selected
        SciPy method (polyphase uses a ceiling, FFT uses the requested rounded length).

    Raises
    ------
    ValueError
        Rates or resulting filter parameters are invalid, or zero-phase padding cannot fit.

    Notes
    -----
    With zero_phase=True, integer sampling rates use polyphase resampling; other rates use FFT resampling followed by filtering. With zero_phase=False, integer downsampling uses causal FIR decimation. See the signal-processing guide.
    """
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
                data = signal.resample(x=data, num=round(len(data) / q))
        elif srate_new > srate_old:
            data = signal.resample(data, round(len(data) / q))
            data = filter_butter(
                data=data,
                srate=srate_new,
                freq_band=[0, srate_old / 2],
                zero_phase=zero_phase,
            )
    return data


def linear_interp(data, N_new) -> np.ndarray:
    """Resample a trace by linear interpolation while preserving both endpoints.

    Parameters
    ----------
    data : numpy.ndarray
        One-dimensional uniformly sampled signal; input is not modified.
    N_new : int
        Number of output samples; must be positive.

    Returns
    -------
    interpolated : numpy.ndarray
        Shape (N_new,), retaining the signal units.

    Raises
    ------
    ValueError
        The input is empty or N_new is invalid.

    Notes
    -----
    This function supplies no antialias lowpass filter; use resample for sampled waveform rate changes.
    """
    points_loc = np.arange(0, len(data))
    points_loc_new = np.linspace(0, len(data) - 1, N_new)
    data_new = np.interp(points_loc_new, points_loc, data)
    return data_new
