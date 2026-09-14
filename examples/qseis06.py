"""Deprecated QSEIS06 tutorial; use qseis2025.py for new calculations.

Retained to build, read and plot existing QSEIS06 workflows. Migration requires
rebuilding the library and validating the replacement backend settings.
"""
from pathlib import Path

import numpy as np

from common import (MECHANISM, MOMENT_NM, finish, parser_for, prepare,
                    require_library_settings, save_waveforms)
from pygrnwang.create_qseis06_bulk import (
    pre_process_qseis06, create_grnlib_qseis06_sequential)
from pygrnwang.read_qseis06 import seek_qseis06
from source_time_function import prepare_qseis_stf, validate_qseis_stf


def main():
    parser = parser_for("qseis06")
    parser.add_argument("--regional", action="store_true",
                        help="Use 300/600/900 km, a 64 s wavelet and Earth flattening")
    parser.set_defaults(output_dir=None)
    args = parser.parse_args()
    if args.output_dir is None:
        directory = "qseis06-regional" if args.regional else "qseis06"
        args.output_dir = Path(__file__).resolve().parent / "output" / directory
    output, library, model, report, started = prepare(args, "QSEIS06")
    dt, window = (4.0, 4092.0) if args.regional else (0.5, 127.5)
    distances = [300.0, 600.0, 900.0] if args.regional else [30.0, 60.0, 90.0]
    wavelet_duration = 16 if args.regional else 4
    wavelet_type = 0 if args.regional else 2
    source_time_function = None
    output_end = 1020.0 if args.regional else window
    output_samples = int(round(output_end / dt)) + 1
    native_samples = int(round(window / dt)) + 1
    if not args.reuse:
        pre_process_qseis06(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], dist_range=[distances[0], distances[-1]],
            delta_dist=distances[0],
            N_each_group=3, time_window=window, sampling_interval=dt,
            wavelet_type=wavelet_type, wavelet_duration=wavelet_duration, time_reduction_velo=0,
            flat_earth_transform=args.regional, path_nd=model, earth_model_layer_num=24,
        )
        if args.regional:
            source_time_function = prepare_qseis_stf(
                library, duration_s=64.0, samples=1024)
        create_grnlib_qseis06_sequential(library, remove_pd=False)
    require_library_settings(
        library, event_depth_list=[10.0], receiver_depth_list=[0.0],
        grn_dist_range=[distances[0], distances[-1]], grn_delta_dist=distances[0],
        sampling_interval=dt, time_window=window, sampling_num=native_samples,
        wavelet_type=wavelet_type, wavelet_duration=wavelet_duration, time_reduction_velo=0,
        flat_earth_transform=args.regional, earth_model_layer_num=24,
        slowness_window=None, wavenumber_sampling_rate=12, anti_alias=0.01,
        free_surface=True,
    )
    if args.regional and args.reuse:
        source_time_function = validate_qseis_stf(library)
    arrays = [MOMENT_NM * seek_qseis06(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="velo" if args.regional else "disp", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances]
    if args.regional:
        # Custom type-0 Green functions are rates; integrate before cropping.
        arrays = [np.cumsum(values, axis=1) * dt for values in arrays]
    arrays = [values[:, :output_samples] for values in arrays]
    save_waveforms(output, report, "disp", arrays, distances, dt, ["E", "N", "U"], "m",
                   expected_samples=output_samples,
                   time_limits=(0.0, output_end) if args.regional else None)
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  output_time_range_s=[0.0, output_end], native_samples=native_samples,
                  earth_model_numeric_rows=24, wavelet_type=wavelet_type,
                  wavelet_duration_samples=wavelet_duration,
                  wavelet_duration_s=wavelet_duration * dt,
                  flat_earth_transform=args.regional, regional=args.regional)
    if args.regional:
        report["source_time_function"] = source_time_function
    finish(output, report, started)


if __name__ == "__main__":
    main()
