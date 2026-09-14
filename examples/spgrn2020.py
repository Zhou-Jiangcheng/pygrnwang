"""A small long-period SPGRN2020 library with P-relative trace windows."""
import json
from pathlib import Path

from common import (MECHANISM, MOMENT_NM, finish, parser_for, prepare,
                    save_waveforms)
from pygrnwang.create_spgrn2020_bulk import (
    pre_process_spgrn2020, create_grnlib_spgrn2020_sequential)
from pygrnwang.read_spgrn2020 import seek_spgrn2020


def main():
    args = parser_for("spgrn2020").parse_args()
    output, library, model, report, started = prepare(args, "SPGRN2020")
    dt, window, before_p = 4.0, 1020.0, 40.0
    if not args.reuse:
        pre_process_spgrn2020(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], spec_time_window=4092.0,
            sampling_interval=dt, max_frequency=0.0625, max_slowness=0.3,
            anti_alias=0.01, gravity_fc=0.0, gravity_harmonic=0,
            cal_sph=1, cal_tor=1, source_radius=0.0, cal_gf=1,
            time_window=window, green_before_p=before_p, source_duration=64.0,
            dist_range=[300.0, 900.0], delta_dist_range=[300.0, 300.0],
            path_nd=model, earth_model_layer_num=None, physical_dispersion=0,
        )
        create_grnlib_spgrn2020_sequential(library)
    info = json.loads((Path(library) / "green_lib_info.json").read_text(encoding="utf-8"))
    distances = info["dist_list"]
    arrays = [MOMENT_NM * seek_spgrn2020(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="disp", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances]
    save_waveforms(output, report, "disp", arrays, distances, dt,
                   ["E", "N", "U"], "m", time_label="Time relative to library P arrival (s)",
                   start_times=[-before_p] * len(distances))
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  spec_time_window_s=4092.0, max_frequency_hz=0.0625, green_before_p_s=before_p, source_duration_s=64.0)
    finish(output, report, started)


if __name__ == "__main__":
    main()
