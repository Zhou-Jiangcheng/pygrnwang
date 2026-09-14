"""Build full-wavefield SPGRN2020 traces and plot native origin-time windows."""
from pathlib import Path

import numpy as np
from scipy.io import FortranFile

from common import (MECHANISM, MOMENT_NM, REGIONAL_SAMPLING_INTERVAL_S, REGIONAL_MAX_FREQUENCY_HZ, REGIONAL_SOURCE_DURATION_S, REGIONAL_STF, finish, parser_for, prepare,
                    require_library_settings, save_waveforms)
from pygrnwang.create_spgrn2020_bulk import (
    pre_process_spgrn2020, create_grnlib_spgrn2020_sequential)
from pygrnwang.read_spgrn2020 import seek_spgrn2020
from spectral_settings import verify_spherical_spectrum


def main():
    args = parser_for("spgrn2020").parse_args()
    output, library, model, report, started = prepare(args, "SPGRN2020")
    dt, window, before_p = REGIONAL_SAMPLING_INTERVAL_S, 1020.0, 40.0
    if not args.reuse:
        pre_process_spgrn2020(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], spec_time_window=4092.0,
            sampling_interval=dt, max_frequency=REGIONAL_MAX_FREQUENCY_HZ, max_slowness=0.0,
            anti_alias=0.01, gravity_fc=0.0, gravity_harmonic=0,
            cal_sph=1, cal_tor=1, source_radius=0.0, cal_gf=1,
            time_window=window, green_before_p=before_p, source_duration=REGIONAL_SOURCE_DURATION_S,
            dist_range=[300.0, 900.0], delta_dist_range=[300.0, 300.0],
            path_nd=model, earth_model_layer_num=None, physical_dispersion=0,
        )
        create_grnlib_spgrn2020_sequential(library)
    info = require_library_settings(library, max_slowness=0.0,
                                    sampling_interval=dt, time_window=window, max_frequency=REGIONAL_MAX_FREQUENCY_HZ, spec_time_window=4092.0,
                                    source_duration=REGIONAL_SOURCE_DURATION_S, green_before_p=before_p,
                                    dist_list=[300.0, 600.0, 900.0])
    distances = info["dist_list"]
    # The native solver rounds P - before_p to whole seconds. Read its actual
    # start-time records instead of assigning every trace a P-relative axis.
    starts = []
    native_path = Path(library) / "GreenFunc" / "10.00" / "0.00" / "grn_d10.00"
    with FortranFile(native_path, "r") as native:
        for distance in distances:
            starts.append(float(native.read_reals(np.float32)[0]))
            for _ in range(10):
                native.read_reals(np.float32)
    arrays = [MOMENT_NM * seek_spgrn2020(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="disp", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances]
    save_waveforms(output, report, "disp", arrays, distances, dt,
                   ["E", "N", "U"], "m", start_times=starts)
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  spec_time_window_s=4092.0, max_frequency_hz=REGIONAL_MAX_FREQUENCY_HZ, green_before_p_s=before_p, source_duration_s=REGIONAL_SOURCE_DURATION_S,
                  max_slowness_s_km=0.0, full_wavefield=True, trace_start_times_s=starts)
    report["physical_source_time_function"] = dict(REGIONAL_STF)
    report["spectral_settings"] = verify_spherical_spectrum(library, "spgrn2020")
    report["source_radius_km"] = 0.0
    finish(output, report, started)


if __name__ == "__main__":
    main()
