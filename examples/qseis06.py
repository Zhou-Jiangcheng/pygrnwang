"""Deprecated QSEIS06 tutorial; use qseis2025.py for new calculations.

Retained to build, read and plot existing QSEIS06 workflows. Migration requires
rebuilding the library and validating the replacement backend settings.
"""
from common import (MECHANISM, MOMENT_NM, finish, parser_for, prepare,
                    save_waveforms)
from pygrnwang.create_qseis06_bulk import (
    pre_process_qseis06, create_grnlib_qseis06_sequential)
from pygrnwang.read_qseis06 import seek_qseis06


def main():
    args = parser_for("qseis06").parse_args()
    output, library, model, report, started = prepare(args, "QSEIS06")
    dt, window = 0.5, 127.5
    distances = [30.0, 60.0, 90.0]
    if not args.reuse:
        pre_process_qseis06(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], dist_range=[30.0, 90.0], delta_dist=30.0,
            N_each_group=3, time_window=window, sampling_interval=dt,
            wavelet_type=2, wavelet_duration=4, time_reduction_velo=0,
            flat_earth_transform=False, path_nd=model, earth_model_layer_num=24,
        )
        create_grnlib_qseis06_sequential(library, remove_pd=False)
    arrays = [MOMENT_NM * seek_qseis06(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="disp", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances]
    save_waveforms(output, report, "disp", arrays, distances, dt, ["E", "N", "U"], "m")
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  earth_model_numeric_rows=24, wavelet_type=2, wavelet_duration_samples=4)
    finish(output, report, started)


if __name__ == "__main__":
    main()
