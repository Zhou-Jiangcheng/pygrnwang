import numpy as np
from typing import Union

from .read_spgrn2020 import seek_spgrn2020
from .read_spgrn2012 import seek_spgrn2012
from .read_qssp2020 import seek_qssp2020
from .read_qseis2025 import seek_qseis2025


def read_syn(
    method: str,
    path_green: str,
    event_depth_km: float,
    receiver_depth_km: float,
    az_deg: float,
    dist_km: float,
    focal_mechanism: Union[np.ndarray, list],
    srate: float,
    output_type: str = "disp",
    rotate: bool = True,
    before_p: Union[float, None] = None,
    pad_zeros: bool = False,
    shift: bool = False,
    only_seismograms: bool = True,
    model_name: str = "ak135fc",
    green_info: Union[dict, None] = None,
    interpolate_type: int = 0,
    freq_band=None,
    butter_order: int = 4,
    zero_phase: bool = False,
):
    """Read dynamic synthetic waveforms through an explicit backend selector.

    Parameters
    ----------
    method : str
        One of qseis2025, qssp2020, spgrn2020 or spgrn2012; QSEIS06 and EDCMP have separate readers.
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    event_depth_km : float
        Requested source depth in km, positive down.
    receiver_depth_km : float
        Requested receiver depth in km, positive down.
    az_deg : float
        Source-to-receiver azimuth in degrees clockwise from north.
    dist_km : float
        Epicentral distance in km; query within the stored distance grid.
    focal_mechanism : array_like
        Either [strike, dip, rake] in degrees; [M0, strike, dip, rake]; six NED components [Mnn, Mne, Mnd, Mee, Med, Mdd]; or [M0, six components]. Three angles imply unit moment; seven entries normalize the six-component shape to M0. Moments are in N m.
    srate : float
        Positive output sampling rate in Hz.
    output_type : str, optional
        Requested observable; supported values and units are listed in Notes. Default: 'disp'.
    rotate : bool, optional
        Rotate vector output to east, north, up when True; False retains radial, transverse, up. Tensor layouts are specified in Notes. Default: True.
    before_p : float or None, optional
        Seconds before the library P onset at the new first sample. None preserves the native window. Default: None.
    pad_zeros : bool, optional
        Shift to source-origin time using zero padding. Use separately from before_p. Default: False.
    shift : bool, optional
        Correct the time axis using P/S arrivals recomputed for the requested geometry and model. Default: False.
    only_seismograms : bool, optional
        Return just the waveform array when True; False returns the array and six metadata fields. Default: True.
    model_name : str, optional
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library. Default: 'ak135fc'.
    green_info : dict or None, optional
        Preloaded green_lib_info.json mapping; None loads it from path_green. Default: None.
    interpolate_type : int, optional
        0 selects nearest neighbor; 1 interpolates source depth, receiver depth and distance. Returned grid metadata remains nearest neighbor. Default: 0.
    freq_band : sequence of float or None, optional
        Two cutoff frequencies [low, high] in Hz. None disables filtering in readers; a missing corner selects lowpass or highpass. Default: None.
    butter_order : int, optional
        Butterworth filter order. Default: 4.
    zero_phase : bool, optional
        True applies forward/backward filtering; False uses causal filtering. Default: False.

    Returns
    -------
    seismograms : numpy.ndarray
        Shape (C, N): C=3 for vectors, 6 for tensors or 1 for scalar outputs,
        with physical units and component conventions of the selected backend.
    metadata : tuple, conditional
        With only_seismograms=False returns the seven-tuple (seismograms,
        tpts_table, first_p, first_s, grn_dep_source, grn_dep_receiver, grn_dist).
        Arrival times are seconds after origin; first_p/first_s are None unless
        shift=True. Grid coordinates are nearest stored nodes in km.

    Raises
    ------
    ValueError
        method is unsupported or a selected reader rejects its parameters.
    OSError
        The selected reader cannot access required library files.
    KeyError
        Metadata does not describe the selected backend.

    Notes
    -----
    Dispatches unchanged keyword arguments to the selected seek function. Vector rotate=True means east, north, up; native tensor layouts and time origins differ by backend. Consult seek_qseis2025, seek_qssp2020, seek_spgrn2020 or seek_spgrn2012 for supported observables, units and tpts_table availability. No automatic backend detection, grid conversion or physical normalization is performed.
    """
    kwargs = {
        "path_green": path_green,
        "event_depth_km": event_depth_km,
        "receiver_depth_km": receiver_depth_km,
        "az_deg": az_deg,
        "dist_km": dist_km,
        "focal_mechanism": focal_mechanism,
        "srate": srate,
        "output_type": output_type,
        "rotate": rotate,
        "before_p": before_p,
        "pad_zeros": pad_zeros,
        "shift": shift,
        "only_seismograms": only_seismograms,
        "model_name": model_name,
        "green_info": green_info,
        "interpolate_type": interpolate_type,
        "freq_band": freq_band,
        "butter_order": butter_order,
        "zero_phase": zero_phase,
    }

    if method == "qseis2025":
        return seek_qseis2025(**kwargs)
    elif method == "qssp2020":
        return seek_qssp2020(**kwargs)
    elif method == "spgrn2020":
        return seek_spgrn2020(**kwargs)
    elif method == "spgrn2012":
        return seek_spgrn2012(**kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")


if __name__ == "__main__":
    pass
