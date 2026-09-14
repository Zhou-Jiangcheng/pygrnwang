"""Build QSEIS2025 introductory or regional traces; optionally include tensors."""
from pathlib import Path

import numpy as np

from common import (MECHANISM, MOMENT_NM, REGIONAL_SAMPLING_INTERVAL_S, REGIONAL_STF, finish, parser_for, prepare,
                    require_library_settings, save_waveforms)
from pygrnwang.create_qseis2025_bulk import (
    pre_process_qseis2025, create_grnlib_qseis2025_sequential)
from pygrnwang.read_qseis2025 import get_outfile_name_list, seek_qseis2025
from source_time_function import prepare_qseis_stf, validate_qseis_stf
from spectral_settings import qseis_spectral_settings


def main():
    parser = parser_for("qseis2025")
    parser.add_argument("--regional", action="store_true",
                        help="Use 300/600/900 km, a 64 s wavelet and Earth flattening")
    parser.add_argument("--point-source", action="store_true",
                        help="Regional control run: disable the default Gaussian spatial source smoothing")
    parser.set_defaults(output_dir=None)
    parser.add_argument("--observables", choices=("disp", "all"), default="disp",
                        help="all also computes strain and stress")
    args = parser.parse_args()
    if args.point_source and not args.regional:
        parser.error("--point-source requires --regional")
    source_radius_ratio = 0.0 if args.point_source else 0.05
    if args.output_dir is None:
        directory = "qseis2025-regional" if args.regional else "qseis2025"
        if args.point_source:
            directory += "-point-source"
        args.output_dir = Path(__file__).resolve().parent / "output" / directory
    output, library, model, report, started = prepare(args, "QSEIS2025")
    dt, window = (REGIONAL_SAMPLING_INTERVAL_S, 4092.0) if args.regional else (0.5, 127.5)
    output_end = 1020.0 if args.regional else 100.0
    output_samples = int(round(output_end / dt)) + 1  # Include the final sample.
    native_samples = int(round(window / dt)) + 1
    distances = [300.0, 600.0, 900.0] if args.regional else [30.0, 60.0, 90.0]
    wavelet_duration = 16 if args.regional else 4
    wavelet_type = 0 if args.regional else 2
    source_time_function = None
    if not args.reuse:
        pre_process_qseis2025(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], dist_range=[distances[0], distances[-1]],
            delta_dist=distances[0],
            N_each_group=3, time_window=window, sampling_interval=dt, source_radius_ratio=source_radius_ratio,
            output_observables=([1, 0, 1, 1, 0] if args.observables == "all"
                                else [1, 0, 0, 0, 0]),
            wavelet_type=wavelet_type, wavelet_duration=wavelet_duration, time_reduction_velo=0,
            flat_earth_transform=args.regional, path_nd=model, earth_model_layer_num=24,
        )
        if args.regional:
            source_time_function = prepare_qseis_stf(
                library, duration_s=64.0, samples=1024)
        create_grnlib_qseis2025_sequential(library, remove_pd=False)
    require_library_settings(
        library, event_depth_list=[10.0], receiver_depth_list=[0.0],
        grn_dist_range=[distances[0], distances[-1]], grn_delta_dist=distances[0],
        sampling_interval=dt, time_window=window, sampling_num=native_samples,
        wavelet_type=wavelet_type, wavelet_duration=wavelet_duration, time_reduction_velo=0,
        flat_earth_transform=args.regional, earth_model_layer_num=24,
        slowness_window=None, wavenumber_sampling_rate=12, anti_alias=0.01,
        free_surface=0,
    )
    native_input = Path(library) / "10.00" / "0.00" / "0_0" / "grn.inp"
    lines = native_input.read_text(encoding="utf-8-sig").splitlines()
    headings = [i for i, line in enumerate(lines) if "WAVENUMBER INTEGRATION PARAMETERS" in line]
    if len(headings) != 1:
        raise ValueError("Expected one native wavenumber section")
    data = [line.split("#", 1)[0].strip() for line in lines[headings[0] + 1:]]
    data = [line for line in data if line]
    controls = [float(value) for value in data[1].split()]
    if controls != [1e-6, source_radius_ratio]:
        raise ValueError("Existing native spatial-source settings differ; use a fresh --output-dir")
    if args.regional and args.reuse:
        source_time_function = validate_qseis_stf(library)
    observables = ("disp", "strain", "stress") if args.observables == "all" else ("disp",)
    if args.reuse:
        # Backend metadata does not store observable flags. Check the requested
        # binary outputs before reading a displacement-only library as tensors.
        native_dir = Path(library) / "10.00" / "0.00" / "0_0"
        required = []
        for observable in observables:
            psv, sh = get_outfile_name_list(observable)
            required.extend(native_dir / ("grn_%s.bin" % name) for name in psv + sh)
        if any(not path.is_file() for path in required):
            raise ValueError("Requested outputs are absent. Recalculate in a fresh "
                             "--output-dir without --reuse using --observables all.")
    for observable in observables:
        read_type = ({"disp": "velo", "strain": "strain_rate", "stress": "stress_rate"}
                     [observable] if args.regional else observable)
        arrays = [MOMENT_NM * seek_qseis2025(
            path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
            az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
            srate=1 / dt, output_type=read_type, rotate=True,
            before_p=None, shift=False, pad_zeros=False,
        ) for distance in distances]
        if args.regional:
            # Integrate every custom-STF rate with the same origin-time rule.
            arrays = [np.cumsum(values, axis=1) * dt for values in arrays]
        arrays = [values[:, :output_samples] for values in arrays]
        labels = ["E", "N", "U"] if observable == "disp" else ["EE", "EN", "EU", "NN", "NU", "UU"]
        unit = {"disp": "m", "strain": "1", "stress": "Pa"}[observable]
        save_waveforms(output, report, observable, arrays, distances, dt, labels, unit,
                       expected_samples=output_samples, time_limits=(0.0, output_end))
    report.update(sampling_interval_s=dt, max_frequency_hz=0.5 / dt, time_window_s=window,
                  output_time_range_s=[0.0, output_end], distances_km=distances,
                  native_samples=native_samples, earth_model_numeric_rows=24,
                  wavelet_type=wavelet_type, wavelet_duration_samples=wavelet_duration,
                  wavelet_duration_s=wavelet_duration * dt,
                  flat_earth_transform=args.regional, regional=args.regional, source_radius_ratio=source_radius_ratio, point_source=args.point_source)
    if args.regional:
        report["source_time_function"] = source_time_function
        report["physical_source_time_function"] = dict(REGIONAL_STF)
        report["spectral_settings"] = qseis_spectral_settings(library)
    finish(output, report, started)


if __name__ == "__main__":
    main()
