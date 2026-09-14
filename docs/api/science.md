# Models, mechanisms, geometry and signals

These utilities make axis, material, sampling and source conventions explicit. Model-building depths are generally kilometres, while the geographic Cartesian conversion functions take depths and offsets in **metres**. A six-component moment tensor uses NED; a rotated waveform vector uses ENU.

Only the objects listed here are included in the documented support surface. Additional helpers are catalogued in the advanced reference.

## focal_mechanism

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.check_convert_fm
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.plane2mt
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.convert_mt_axis
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.tensor2full_tensor_matrix
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.moment_from_moment_tensor
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.cal_m0_from_mt
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.mt2plane
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.plane2nd
```

```{eval-rst}
.. autofunction:: pygrnwang.focal_mechanism.plane2tbp
```

## geo

```{eval-rst}
.. autofunction:: pygrnwang.geo.rotate_2d_points
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.rotate_rtz_to_enz
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.create_rotate_z_mat
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.rotate_symmetric_tensor_series
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.geo_2_r_earth
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.r_earth_2_geo
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.convert_axis_delta_geo2ned
```

```{eval-rst}
.. autofunction:: pygrnwang.geo.convert_axis_delta_ned2geo
```

## signal_process

```{eval-rst}
.. autofunction:: pygrnwang.signal_process.taper
```

```{eval-rst}
.. autofunction:: pygrnwang.signal_process.cal_sos
```

```{eval-rst}
.. autofunction:: pygrnwang.signal_process.filter_butter
```

```{eval-rst}
.. autofunction:: pygrnwang.signal_process.resample
```

```{eval-rst}
.. autofunction:: pygrnwang.signal_process.linear_interp
```

## utils

```{eval-rst}
.. autofunction:: pygrnwang.utils.read_nd
```

```{eval-rst}
.. autofunction:: pygrnwang.utils.read_material_nd
```

```{eval-rst}
.. autofunction:: pygrnwang.utils.read_layerd_material
```

```{eval-rst}
.. autofunction:: pygrnwang.utils.convert_earth_model_nd2inp
```

```{eval-rst}
.. autofunction:: pygrnwang.utils.convert_earth_model_nd2nd_without_Q
```

```{eval-rst}
.. autofunction:: pygrnwang.utils.create_stf
```

```{eval-rst}
.. autofunction:: pygrnwang.utils.cal_grid
```

## create_nd_by_crust1_ak135

```{eval-rst}
.. autofunction:: pygrnwang.create_nd_by_crust1_ak135.create_nd_by_crust1_ak135
```

## crust1

```{eval-rst}
.. autoclass:: pygrnwang.crust1.CrustModel
```

```{eval-rst}
.. automethod:: pygrnwang.crust1.CrustModel.get_point
```
