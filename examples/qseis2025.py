"""Build, read and plot three QSEIS2025 traces; optionally include tensors."""
from common import (MECHANISM, MOMENT_NM, finish, parser_for, prepare,
                    save_waveforms)
from pygrnwang.create_qseis2025_bulk import (
    pre_process_qseis2025, create_grnlib_qseis2025_sequential)
from pygrnwang.read_qseis2025 import seek_qseis2025


def main():
    parser = parser_for("qseis2025")
    parser.add_argument("--observables", choices=("disp", "all"), default="disp",
                        help="all also computes strain and stress")
    args = parser.parse_args()
    output, library, model, report, started = prepare(args, "QSEIS2025")
    dt, window = 0.5, 127.5  # 256 samples; QSEIS durations below are sample counts.
    distances = [30.0, 60.0, 90.0]
    if not args.reuse:
        pre_process_qseis2025(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], dist_range=[30.0, 90.0], delta_dist=30.0,
            N_each_group=3, time_window=window, sampling_interval=dt,
            output_observables=([1, 0, 1, 1, 0] if args.observables == "all"
                                else [1, 0, 0, 0, 0]),
            wavelet_type=2, wavelet_duration=4, time_reduction_velo=0,
            flat_earth_transform=False, path_nd=model, earth_model_layer_num=24,
        )
        create_grnlib_qseis2025_sequential(library, remove_pd=False)
    observables = ("disp", "strain", "stress") if args.observables == "all" else ("disp",)
    for observable in observables:
        arrays = [MOMENT_NM * seek_qseis2025(
            path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
            az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
            srate=1 / dt, output_type=observable, rotate=True,
            before_p=None, shift=False, pad_zeros=False,
        ) for distance in distances]
        labels = ["E", "N", "U"] if observable == "disp" else ["EE", "EN", "EU", "NN", "NU", "UU"]
        unit = {"disp": "m", "strain": "1", "stress": "Pa"}[observable]
        save_waveforms(output, report, observable, arrays, distances, dt, labels, unit)
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  earth_model_numeric_rows=24, wavelet_type=2, wavelet_duration_samples=4)
    finish(output, report, started)


if __name__ == "__main__":
    main()
