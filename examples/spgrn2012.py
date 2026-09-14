"""Deprecated SPGRN2012 tutorial; use spgrn2020.py for new calculations.

Retained for existing long-period libraries with explicit travel-time tables.
Migration requires rebuilding the library and validating sampling and time origin.
"""
import json
from pathlib import Path

from common import (MECHANISM, MOMENT_NM, finish, parser_for, prepare,
                    save_waveforms)
from pygrnwang.create_spgrn2012_bulk import (
    pre_process_spgrn2012, create_grnlib_spgrn2012_sequential)
from pygrnwang.pytaup import create_tpts_table
from pygrnwang.read_spgrn2012 import seek_spgrn2012


def main():
    args = parser_for("spgrn2012").parse_args()
    output, library, model, report, started = prepare(args, "SPGRN2012")
    dt, window = 4.0, 1020.0
    t0, v0 = -40.0, 10.0
    if not args.reuse:
        pre_process_spgrn2012(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], spec_time_window=4092.0,
            sampling_interval=dt, max_frequency=0.0625, max_slowness=0.3,
            anti_alias=0.01, gravity_fc=0.0, gravity_harmonic=0,
            cal_sph=1, cal_tor=1, source_radius=0.0, cal_gf=1,
            time_window=window, t0=t0, v0=v0, source_duration=64.0,
            dist_range=[300.0, 900.0], delta_dist_range=[300.0, 300.0],
            path_nd=model, earth_model_layer_num=None, physical_dispersion=0,
        )
        create_grnlib_spgrn2012_sequential(library)
    info = json.loads((Path(library) / "green_lib_info.json").read_text(encoding="utf-8"))
    distances = info["dist_list"]  # SPGRN chooses and records its distance grid.
    if not args.reuse:
        # The serial builder does not create these tables; the reader needs them.
        create_tpts_table(str(Path(library) / "GreenFunc"), 10.0, 0.0,
                          distances, info["path_nd_without_Q"], False)
    arrays = [MOMENT_NM * seek_spgrn2012(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="disp", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances]
    save_waveforms(output, report, "disp", arrays, distances, dt,
                   ["E", "N", "U"], "m",
                   start_times=[t0 + distance / v0 for distance in distances])
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  spec_time_window_s=4092.0, max_frequency_hz=0.0625, t0_s=t0, v0_km_s=v0, source_duration_s=64.0)
    finish(output, report, started)


if __name__ == "__main__":
    main()
