import os
import subprocess


def read_qseis_by_fortran(
    path_green_lib,
    event_depth_in_km=None,
    az_in_deg=None,
    dist_in_km=None,
    focal_mechanism=None,
    srate=1.0,
    zero_phase=False,
    rotate=True,
    time_reduction_slowness=8,
    before_p=None,
    pad_zeros=False,
    shift=False,
    only_seismograms=True,
    model_name="ak135fc",
    nrmax=200,
):
    pass
