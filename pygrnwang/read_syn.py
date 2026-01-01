import numpy as np
from typing import Union

from .read_spgrn import seek_spgrn2020
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
):
    """
    Unified interface to read synthetic seismograms using different methods.

    :param path_green: Root directory of the data.
    :param event_depth_km: Event depth in km.
    :param receiver_depth_km: Receiver depth in km.
    :param az_deg: Azimuth in degrees.
    :param dist_km: Epicentral distance in km.
    :param focal_mechanism: [strike, dip, rake] or [M11, M12, M13, M22, M23, M33].
    :param srate: Sampling rate in Hz.
    :param output_type: disp | velo | strain | strain_rate |
            stress | stress_rate | rota | rota_rate.
    :param before_p: Time before P-wave.
    :param pad_zeros: Pad with zeros.
    :param shift: Shift seismograms based on tpts.
    :param rotate: Rotate rtz2ned.
    :param only_seismograms: Return only seismograms.
    :param model_name: Model name.
    :param green_info: Green's function library info.
    :param interpolate_type:
            0 for nearest neighbor,
            1 for trilinear interpolation (Source Depth, Receiver Depth, Distance).
    :return: (
            seismograms_resample,
            tpts_table,
            first_p,
            first_s,
            grn_dep_source,
            grn_dep_receiver,
            grn_dist,
        )
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
    }

    if method == "qseis2025":
        return seek_qseis2025(**kwargs)
    elif method == "qssp2020":
        return seek_qssp2020(**kwargs)
    elif method == "spgrn2020":
        return seek_spgrn2020(**kwargs)
    else:
        raise ValueError(f"Unknown method: {method}")


if __name__ == '__main__':
    pass
