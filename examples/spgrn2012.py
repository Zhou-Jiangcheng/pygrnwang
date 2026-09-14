"""Deprecated SPGRN2012 tutorial with an explicitly matched physical source.

Retained for existing long-period workflows; use SPGRN2020 for new calculations.
Here a full-period native impulse velocity is forward-convolved with the same
64 s moment-rate pulse used by the other regional examples, then integrated.
"""
from pathlib import Path

import numpy as np

from common import (MECHANISM, MOMENT_NM, REGIONAL_MAX_FREQUENCY_HZ,
                    REGIONAL_SAMPLING_INTERVAL_S, REGIONAL_SOURCE_DURATION_S,
                    REGIONAL_STF, finish, parser_for, prepare,
                    require_library_settings, save_waveforms)
from pygrnwang.create_spgrn2012_bulk import (
    pre_process_spgrn2012, create_grnlib_spgrn2012_sequential)
from pygrnwang.pytaup import create_tpts_table
from pygrnwang.read_spgrn2012 import seek_spgrn2012
from spectral_settings import verify_spherical_spectrum
from spherical_source_time_function import (apply_spgrn2012_stf, inspect_spgrn2012_native,
                                           save_spgrn2012_stf, validate_spgrn2012_stf)


def main():
    args = parser_for("spgrn2012").parse_args()
    output, library, model, report, started = prepare(args, "SPGRN2012")
    dt, native_window, output_window = REGIONAL_SAMPLING_INTERVAL_S, 4092.0, 1020.0
    t0, v0 = -40.0, 10.0
    output_samples = int(round(output_window / dt)) + 1
    if not args.reuse:
        pre_process_spgrn2012(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], spec_time_window=native_window,
            sampling_interval=dt, max_frequency=REGIONAL_MAX_FREQUENCY_HZ, max_slowness=0.0,
            anti_alias=0.01, gravity_fc=0.0, gravity_harmonic=0,
            cal_sph=1, cal_tor=1, source_radius=0.0, cal_gf=1,
            time_window=native_window, t0=t0, v0=v0, source_duration=0.0,
            dist_range=[300.0, 900.0], delta_dist_range=[300.0, 300.0],
            path_nd=model, earth_model_layer_num=None, physical_dispersion=0,
        )
        create_grnlib_spgrn2012_sequential(library)
    info = require_library_settings(
        library, event_depth_list=[10.0], receiver_depth_list=[0.0],
        spec_time_window=native_window, time_window=native_window, sampling_interval=dt,
        samples_num=1024, source_duration=0.0, max_frequency=REGIONAL_MAX_FREQUENCY_HZ,
        max_slowness=0.0, anti_alias=0.01, physical_dispersion=0,
        gravity_fc=0.0, gravity_harmonic=0, cal_sph=1, cal_tor=1,
        source_radius=0.0, t0=t0, v0=v0, dist_list=[300.0, 600.0, 900.0],
    )
    native = inspect_spgrn2012_native(library)
    if args.reuse:
        validate_spgrn2012_stf(library)
    distances = info["dist_list"]
    if not args.reuse:
        # The serial builder does not create these tables; the reader needs them.
        create_tpts_table(str(Path(library) / "GreenFunc"), 10.0, 0.0,
                          distances, info["path_nd_without_Q"], False)
    raw_velocity = np.array([MOMENT_NM * seek_spgrn2012(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="velo", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances])
    matched_velocity, source = apply_spgrn2012_stf(raw_velocity, native)
    # Preserve the whole FFT period during source convolution. Integrate once,
    # then export the first 256 samples on each trace's actual origin-time axis.
    arrays = np.cumsum(matched_velocity, axis=2)[:, :, :output_samples] * dt
    save_waveforms(output, report, "disp", arrays, distances, dt,
                   ["E", "N", "U"], "m", start_times=native["trace_start_times_s"],
                   expected_samples=output_samples)
    source = save_spgrn2012_stf(output, native, raw_velocity, matched_velocity, source)
    validate_spgrn2012_stf(library)
    report.update(sampling_interval_s=dt, time_window_s=native_window,
                  native_samples=1024, output_window_s=output_window,
                  output_time_ranges_s=[[start, start + output_window]
                                        for start in native["trace_start_times_s"]],
                  distances_km=distances, spec_time_window_s=native_window,
                  max_frequency_hz=REGIONAL_MAX_FREQUENCY_HZ,
                  last_retained_frequency_hz=native["last_retained_frequency_hz"],
                  nyquist_bin_zero=True, t0_s=t0, v0_km_s=v0,
                  trace_start_times_s=native["trace_start_times_s"],
                  native_source_duration_s=0.0, source_duration_s=REGIONAL_SOURCE_DURATION_S,
                  effective_source_duration_s=REGIONAL_SOURCE_DURATION_S,
                  physical_source_time_function=dict(REGIONAL_STF),
                  full_wavefield=True, max_slowness_s_km=0.0,
                  source_radius_km=0.0,
                  spectral_settings=verify_spherical_spectrum(library, 'SPGRN2012'),
                  source_time_function=source)
    finish(output, report, started)


if __name__ == "__main__":
    main()
