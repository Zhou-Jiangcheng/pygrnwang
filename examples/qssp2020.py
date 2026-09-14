"""Build QSSP2020 spectra first, then synthesize a small displacement library."""
from common import (MECHANISM, MOMENT_NM, REGIONAL_SAMPLING_INTERVAL_S, REGIONAL_MAX_FREQUENCY_HZ, REGIONAL_SOURCE_DURATION_S, REGIONAL_STF, finish, parser_for, prepare,
                    require_library_settings, save_waveforms)
from pygrnwang.create_qssp2020_bulk import (
    pre_process_qssp2020, create_grnlib_qssp2020_sequential)
from pygrnwang.read_qssp2020 import seek_qssp2020
from spectral_settings import verify_spherical_spectrum


def main():
    args = parser_for("qssp2020").parse_args()
    output, library, model, report, started = prepare(args, "QSSP2020")
    dt, window = REGIONAL_SAMPLING_INTERVAL_S, 1020.0
    distances = [300.0, 600.0, 900.0]
    if not args.reuse:
        pre_process_qssp2020(
            processes_num=1, path_green=library, event_depth_list=[10.0],
            receiver_depth_list=[0.0], spec_time_window=4092.0,
            sampling_interval=dt, max_frequency=REGIONAL_MAX_FREQUENCY_HZ, max_slowness=0.3,
            anti_alias=0.01, turning_point_filter=0, turning_point_d1=0.0,
            turning_point_d2=6371.0, free_surface_filter=1,
            gravity_fc=0.0, gravity_harmonic=0, cal_sph=1, cal_tor=1,
            min_harmonic=2000, max_harmonic=8000, source_radius=0.0,
            source_duration=REGIONAL_SOURCE_DURATION_S, output_observables=[1] + [0] * 10,
            time_window=window, time_reduction=0.0,
            dist_range=[300.0, 900.0], delta_dist=300.0,
            path_nd=model, earth_model_layer_num=None, physical_dispersion=0,
        )
        # Required for a fresh library: time-domain synthesis consumes GreenSpec.
        create_grnlib_qssp2020_sequential(library, cal_spec=True, remove_pd=False)
    require_library_settings(library, min_harmonic=2000, max_harmonic=8000,
                             sampling_interval=dt, time_window=window, max_frequency=REGIONAL_MAX_FREQUENCY_HZ, spec_time_window=4092.0,
                             source_duration=REGIONAL_SOURCE_DURATION_S, max_slowness=0.3,
                             grn_dist_range=[300.0, 900.0], grn_delta_dist=300.0)
    arrays = [MOMENT_NM * seek_qssp2020(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        srate=1 / dt, output_type="disp", rotate=True,
        before_p=None, shift=False, pad_zeros=False,
    ) for distance in distances]
    save_waveforms(output, report, "disp", arrays, distances, dt, ["E", "N", "U"], "m")
    report.update(sampling_interval_s=dt, time_window_s=window, distances_km=distances,
                  spec_time_window_s=4092.0, max_frequency_hz=REGIONAL_MAX_FREQUENCY_HZ, harmonic_bounds=[2000, 8000], source_duration_s=REGIONAL_SOURCE_DURATION_S)
    report["physical_source_time_function"] = dict(REGIONAL_STF)
    report["spectral_settings"] = verify_spherical_spectrum(library, "qssp2020")
    report["source_radius_km"] = 0.0
    finish(output, report, started)


if __name__ == "__main__":
    main()
