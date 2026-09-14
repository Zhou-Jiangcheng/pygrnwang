# Advanced implementation reference

The supported API is an explicit selection, not every name lacking an underscore. The following helpers remain available for inspecting formats, customizing native input files and reading old workflows. Their signatures are recorded here, but they do not carry the same compatibility guarantee as the supported reference.

Start with the backend tutorials and supported builders/readers. Low-level calls require the exact directory structure and binary layout produced by the corresponding program version. A helper from one backend must not be substituted into another backend solely because its arguments look similar.

Private names, backup modules, input-template string constants and command wrappers are excluded. The command-line guide describes installed executable wrappers. Historical utility aliases for create_rotate_z_mat and rotate_symmetric_tensor_series continue to point to the functions documented under geo.

## pygrnwang.create_edcmp

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_edcmp.py).

```python
create_inp_edcmp2(path_green: str, event_depth: float, obs_depth: float, dist_range: tuple[float, float], delta_dist: float, mt_ind: int, output_observables: tuple[int, int, int, int], layered: bool=True, lam: float=30516224000.0, mu: float=33701888000.0)
```

```python
call_edcmp2(event_depth, obs_depth, mt_ind, path_green, check_finished=False)
```

```python
convert_edcmp2(path_sub_dir, output_type_ind, remove=False)
```

## pygrnwang.create_edgrn

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_edgrn.py).

```python
create_inp_edgrn2(path_green, obs_depth, grn_dist_range, grn_delta_dist, grn_source_depth_range, grn_delta_source_depth, wavenumber_sampling_rate=12, path_nd=None, earth_model_layer_num=None)
```

```python
call_edgrn2(obs_depth, path_green, check_finished=False)
```

## pygrnwang.create_qseis06_bulk

**Deprecated backend: QSEIS06.** For new calculations use QSEIS2025, rebuild the library and revalidate settings and results.

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_qseis06_bulk.py).

```python
create_order_ind(order, diff_accu_order)
```

## pygrnwang.create_qseis06

**Deprecated backend: QSEIS06.** For new calculations use QSEIS2025, rebuild the library and revalidate settings and results.

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_qseis06.py).

```python
create_dir_qseis06(path_green, event_depth, receiver_depth, dist_range, delta_dist, N_each_group=100, order=0)
```

```python
create_inp_qseis06(path_sub_dir, event_depth, receiver_depth, dist_range, delta_dist, N_dist, N_dist_group, N_each_group, time_window, sampling_interval, slowness_int_algorithm=0, slowness_window=None, time_reduction_velo=0, wavenumber_sampling_rate=12, anti_alias=0.01, free_surface=True, wavelet_duration=4, wavelet_type=1, flat_earth_transform=True, path_nd=None, earth_model_layer_num=None, order=0)
```

```python
create_inp_qseis06_points(path_green, event_depth, receiver_depth, n_group, points, time_window, sampling_interval, slowness_int_algorithm=0, slowness_window=None, time_reduction_velo=0, wavenumber_sampling_rate=2, anti_alias=0.01, free_surface=True, wavelet_duration=4, wavelet_type=1, flat_earth_transform=True, path_nd=None, earth_model_layer_num=None, order=0)
```

```python
call_qseis06(event_depth, receiver_depth, n_group, order, path_green, check_finished=False)
```

```python
convert_pd2bin_qseis06(path_greenfunc, remove=False)
```

## pygrnwang.create_qseis2025

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_qseis2025.py).

```python
create_dir_qseis2025(path_green, event_depth, receiver_depth, dist_range, delta_dist, N_each_group=500)
```

```python
create_inp_qseis2025(path_green, event_depth, receiver_depth, dist_range, delta_dist, N_dist, N_dist_group, N_each_group, time_window, sampling_interval, output_observables, slowness_int_algorithm=0, eps_estimate_wavenumber=1e-06, source_radius_ratio=0.05, slowness_window=None, time_reduction_velo=0, wavenumber_sampling_rate=12, anti_alias=0.01, free_surface=0, wavelet_duration=4, wavelet_type=1, flat_earth_transform=True, path_nd=None, earth_model_layer_num=None)
```

```python
call_qseis2025(event_depth, receiver_depth, n_group, path_green, check_finished=False)
```

```python
convert_pd2bin_qseis2025(path_greenfunc, remove=False)
```

## pygrnwang.create_qssp2020_bulk

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_qssp2020_bulk.py).

```python
pre_process_spec(processes_num, path_green, event_depth_list, receiver_depth_list, spec_time_window, sampling_interval, max_frequency, max_slowness, anti_alias, turning_point_filter, turning_point_d1, turning_point_d2, free_surface_filter, gravity_fc, gravity_harmonic, cal_sph, cal_tor, min_harmonic, max_harmonic, source_radius, source_duration, time_window, time_reduction, dist_range, delta_dist, path_nd=None, earth_model_layer_num=None, physical_dispersion=0)
```

```python
pre_process_func(processes_num, path_green, event_depth_list, receiver_depth_list, spec_time_window, sampling_interval, max_frequency, max_slowness, anti_alias, turning_point_filter, turning_point_d1, turning_point_d2, free_surface_filter, gravity_fc, gravity_harmonic, cal_sph, cal_tor, min_harmonic, max_harmonic, source_radius, source_duration, output_observables: list, time_window, time_reduction, dist_range, delta_dist, path_nd=None, earth_model_layer_num=None, physical_dispersion=0)
```

```python
remove_dat_files(path_green)
```

## pygrnwang.create_qssp2020

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_qssp2020.py).

```python
create_dir_qssp2020(event_depth, receiver_depth, path_green)
```

```python
create_points(dist_range, delta_dist)
```

```python
create_locs(points, time_reduction)
```

```python
create_inp_qssp2020(mt_com, path_green, event_depth, receiver_depth, spec_time_window, sampling_interval, max_frequency, max_slowness, anti_alias, turning_point_filter, turning_point_d1, turning_point_d2, free_surface_filter, gravity_fc, gravity_harmonic, cal_sph, cal_tor, min_harmonic, max_harmonic, source_radius, cal_gf, source_duration, output_observables: list, time_window, time_reduction, dist_range, delta_dist, path_nd=None, earth_model_layer_num=None, physical_dispersion=0)
```

```python
call_qssp2020(event_depth, receiver_depth, mt_com, path_green, check_finished=False)
```

```python
convert_pd2bin_qssp2020(path_green, event_depth, receiver_depth, output_type_ind)
```

## pygrnwang.create_spgrn2012

**Deprecated backend: SPGRN2012.** For new calculations use SPGRN2020, rebuild the library and revalidate settings, time alignment and results.

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_spgrn2012.py).

```python
create_inp_spgrn2012(path_green, event_depth, receiver_depth, spec_time_window, sampling_interval, max_frequency, max_slowness, anti_alias, gravity_fc, gravity_harmonic, cal_sph, cal_tor, source_radius, cal_gf, time_window, t0, v0, source_duration, dist_range, delta_dist_range, path_nd=None, earth_model_layer_num=None, physical_dispersion=0)
```

```python
call_spgrn2012(event_depth, receiver_depth, path_green, check_finished=False)
```

## pygrnwang.create_spgrn2020_bulk

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_spgrn2020_bulk.py).

```python
update_green_info_lib_json(path_green, event_depth, receiver_depth)
```

## pygrnwang.create_spgrn2020

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/create_spgrn2020.py).

```python
create_dir_spgrn(event_depth, receiver_depth, path_green)
```

```python
create_inp_spgrn2020(path_green, event_depth, receiver_depth, spec_time_window, sampling_interval, max_frequency, max_slowness, anti_alias, gravity_fc, gravity_harmonic, cal_sph, cal_tor, source_radius, cal_gf, time_window, green_before_p, source_duration, dist_range, delta_dist_range, path_nd=None, earth_model_layer_num=None, physical_dispersion=0)
```

```python
call_spgrn2020(event_depth, receiver_depth, path_green, check_finished=False)
```

## pygrnwang.geo

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/geo.py).

```python
rotate_vector_x(v, gamma)
```

```python
rotate_vector_y(v, gamma)
```

```python
rotate_vector_z(v, gamma)
```

```python
spherical_2_cartesian(r, phi, theta)
```

```python
cartesian_2_spherical(x, y, z)
```

```python
convert_sub_faults_geo2ned(sub_faults, source_point, approximate=True)
```

```python
select_df_geo(df: pd.DataFrame, lat_range, lon_range, time_range)
```

```python
geographic_centroid(points)
```

```python
cal_max_dist_from_2d_points(A: np.ndarray, B: np.ndarray)
```

## pygrnwang.plot_seismograms

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/plot_seismograms.py).

```python
plot_seismograms(seismograms, srate, ylabel='u (m)', cut_length=None)
```

## pygrnwang.read_by_qssp

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_by_qssp.py).

```python
hash_read_pars(event_lat, event_lon, event_depth, receiver_lat, receiver_lon, receiver_depth)
```

```python
create_inp_qssp2020_read(read_name, path_green, event_lat, event_lon, event_depth, receiver_lat, receiver_lon, receiver_depth, focal_mechanism, green_info)
```

```python
call_qssp2020_read(path_green, path_inp, check_finished=False)
```

## pygrnwang.read_edcmp

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_edcmp.py).

```python
read_edcmp_raw(path_green, output_type, grn_event_depth, grn_obs_depth, mt_ind)
```

```python
interpolate_values(xmin, xmax, nx, ymin, ymax, ny, v_array, obs_array)
```

## pygrnwang.read_qseis06_diff

**Deprecated backend: QSEIS06.** Use QSEIS2025 direct tensor outputs for new calculations, after rebuilding and validating the library. This backend status does not deprecate generic tensor utilities.

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_qseis06_diff.py).

```python
diff_central_1order(v_array, diff_accu_order)
```

## pygrnwang.read_qseis06

**Deprecated backend: QSEIS06.** Use QSEIS2025 readers with a newly built and validated QSEIS2025 library for new calculations.

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_qseis06.py).

```python
read_time_series_qseis06_bin(path_greenfunc, start_count, sampling_num)
```

```python
read_time_series_qseis06_ascii(path_greenfunc, start_count)
```

```python
synthesize_qseis06(time_series_list, m1, m2)
```

## pygrnwang.read_qseis2025

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_qseis2025.py).

```python
get_outfile_name_list(output_type)
```

```python
read_time_series_qseis2025_ascii(path_greenfunc, start_count, output_type='disp')
```

```python
read_time_series_qseis2025_bin(path_greenfunc, start_count, output_type, sampling_num)
```

```python
synthesize_rzv(time_series, m1)
```

```python
synthesize_t(time_series, m2)
```

```python
get_sorted_grid_params(target, grid_list)
```

## pygrnwang.read_qssp2020

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_qssp2020.py).

```python
read_time_series_qssp2020(path_bin, ind, sampling_num)
```

```python
seek_raw_qssp2020(path_green, event_depth, receiver_depth, dist, output_type='disp', green_info=None)
```

```python
get_sorted_grid_params(target, grid_list)
```

## pygrnwang.read_spgrn2020

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/read_spgrn2020.py).

```python
synthesize_spgrn(az_in_deg, time_series, focal_mechanism)
```

```python
read_spgrn_data_by_index(path_grn_data, dist_index, green_info)
```

```python
read_spgrn_data_two_indices(path_grn_data, dist_idx_low, dist_idx_high, green_info)
```

```python
get_sorted_grid_params(target, grid_list)
```

## pygrnwang.utils

[Read the implementation](https://github.com/Zhou-Jiangcheng/pygrnwang/blob/main/pygrnwang/utils.py).

```python
read_source_array(source_inds, path_input, shift2corner=False, source_shapes=None)
```

```python
group(inp_list, num_in_each_group)
```

```python
shift_green2real_tpts(seismograms, tpts_table, green_before_p, srate, event_depth_km, dist_in_km, receiver_depth_km=0, model_name='ak135')
```

```python
group_planes(strike_array)
```

```python
reshape_sub_faults(sub_faults, num_strike, num_dip)
```

```python
call_exe(path_inp, path_finished, name)
```
