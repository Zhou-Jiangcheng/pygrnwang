# QSEIS2025: displacement, strain and stress

QSEIS2025 calculates dynamic Green's functions in a layered half-space and
provides direct strain/stress and rotation kernels alongside vectors. Use
this tutorial for a complete first calculation and for the introductory
tensor workflow.

## Build and inspect displacement

From the repository root:

```console
python examples/qseis2025.py
```

The script prepares the model and input files, computes the library
sequentially, reads three distances and saves an ENU displacement figure.
E/N/U in example labels is the same east/north/up convention called ENZ
by the API.

In the default near-distance mode, source depth is 10 km, receiver depth
0 km, and distances are 30, 60 and 90 km. The 0.5 s sample interval and
127.5 s window produce 256 native samples.
After synthesis, every saved near-distance waveform and figure is cropped to
0–100 s inclusive (201 samples); the library keeps its full native window.
The model uses the first 24 numeric rows and disables the flat-Earth
transformation for this small half-space example. Moment is `10^15 N m`.

```{literalinclude} ../../examples/qseis2025.py
:language: python
:caption: Complete QSEIS2025 tutorial
```

```{figure} ../_static/examples/qseis2025.png
:alt: Three-component QSEIS2025 displacement at three epicentral distances.

Example displacement in metres with a fixed 0–100 s horizontal axis. Time is
measured from source origin because this example uses zero time reduction.
```

Results are under `examples/output/qseis2025/`. `disp.npz` contains the
three distance traces, component labels, seconds and units; its waveform
array has shape `(3, 3, 201)` and includes both 0 s and 100 s.
`summary.json` records checks and run details. The solver library still
contains 256 samples per trace.
Use `python examples/qseis2025.py --reuse` to read the completed library
again without launching the solver.

## Add tensor outputs

Use a separate output directory so the selected quantities match the library:

```console
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025_tensors
```

This enables displacement, strain and stress; it writes `strain.npz/png`
and `stress.npz/png` in addition to displacement. The tensor arrays have
shape `(3, 6, 201)` and components EE, EN, EU, NN, NU, UU. All three
observables are cropped to 0–100 s inclusive before saving arrays and plots.
Strain is dimensionless and stress is in Pa after the script's explicit
moment scaling.

```{figure} ../_static/examples/qseis2025-strain.png
:alt: Six ENU strain tensor components from the QSEIS2025 tensor tutorial.

Dimensionless tensor strain at 30, 60 and 90 km, cropped to 0–100 s since
source origin; shear entries are tensor strains rather than engineering shear.
```

```{figure} ../_static/examples/qseis2025-stress.png
:alt: Six ENU stress tensor components from the QSEIS2025 tensor tutorial.

Stress in Pa for the same moment and geometry, cropped to 0–100 s since source
origin. The surface receiver's traction components can be zero under the
free-surface boundary condition.
```

The five `output_observables` positions are:

| Zero-based position | Family | Reader requests |
| --- | --- | --- |
| 0 | Vector motion | `disp`, `velo`, `acce` |
| 1 | Fractional volume change | `volume` |
| 2 | Strain | `strain`, `strain_rate` |
| 3 | Stress | `stress`, `stress_rate` |
| 4 | Rotation | `rota`, `rota_rate` |

Thus `[1, 0, 1, 1, 0]` enables the demonstrated tensor calculation.
The physical stored quantity depends on `wavelet_type`. Here type 2 is
a tapered Heaviside, so the stored non-rate kernels include displacement,
strain and stress. `wavelet_duration=4` means **four samples**, or 2 s,
not four seconds.

## Regional waveforms at 300, 600 and 900 km

The regional option preserves the short default tutorial and selects
a separate calculation at 300, 600 and 900 km:

```console
python examples/qseis2025.py --regional
```

The source remains 10 km deep and receivers remain at the surface.
`flat_earth_transform=True` is enabled with the same 24 numeric model
rows. The interval is 4 s and the native window is 4092 s, giving
1024 library samples. The saved waveforms and plots cover 0–1020 s
inclusive, or 256 samples; their time zero is source origin.
Results default to `examples/output/qseis2025-regional/`, where
`disp.npz` has shape `(3, 3, 256)`. The regional QSEIS, SPGRN and QSSP
examples share the 4 s interval and 0.125 Hz Nyquist limit. QSEIS derives
its frequency range from sampling rather than a separate `max_frequency`
argument. The 1024-point FFT computes through
`511/4096 = 0.124755859375 Hz`; the Nyquist bin is zero.
All regional examples use strike/dip/rake 30°/45°/90°, moment `10^15 N m`
and azimuth 30°.

The source uses `wavelet_type=0, wavelet_duration=16` and 1024
custom moment-rate samples over 0–64 s. Its physical target is the same
normalized squared half-sinusoid used by SPGRN2020 and QSSP2020.
The example compensates for QSEIS's real-frequency implementation by
multiplying the input samples by `exp(2*pi*fi*t)`, where `fi<0` is
the solver's numerical-damping frequency. It does not renormalize the
written samples: their area is approximately 0.9647094, whereas the
effective pulse after damping correction has unit area and centroid 32 s.
See the [STF verification](../guides/backend-comparison.md#matching-the-effective-source-time-function).

Custom type 0 does not invoke automatic rate conversion in the ordinary
reader. The example explicitly reads velocity, then integrates once
using `cumsum * dt` to obtain displacement. QSEIS duration is still
measured in native time samples: sixteen 4 s intervals define the
64 s support, while the independent custom array has 1024 nodes.
The near-distance default remains type 2 with four samples at 0.5 s.

```{figure} ../_static/examples/qseis2025-regional.png
:alt: QSEIS2025 regional displacement at 300, 600 and 900 km.

Regional displacement in metres over 0–1020 s since source origin.
The calculation retains its full 4092 s native library window.
```

To include regional strain and stress, use:

```console
python examples/qseis2025.py --regional --observables all --output-dir examples/output/qseis2025-regional-tensors
```

This additionally saves `strain.npz/png` and `stress.npz/png`, with
tensor shape `(3, 6, 256)` and component order EE, EN, EU, NN, NU, UU.
Their units and moment normalization are the same as in the near-distance
tensor example. For the regional type-0 library, the script reads
`strain_rate` and `stress_rate` and integrates each exactly once using
`cumsum * dt` before saving the non-rate tensors. Requesting `strain`
or `stress` directly from a generic type-0 reader does not perform
that conversion automatically.

```{figure} ../_static/examples/qseis2025-regional-strain.png
:alt: Six QSEIS2025 regional ENU strain components at three distances.

Dimensionless regional strain over 0–1020 s. Shear entries are tensor strains.
```

```{figure} ../_static/examples/qseis2025-regional-stress.png
:alt: Six QSEIS2025 regional ENU stress components at three distances.

Regional stress in Pa over 0–1020 s for the same source and model.
```

The effective 64 s source pulse and frequency band are shared with the
spherical examples. The default QSEIS2025 regional calculation retains
`source_radius_ratio=0.05`, matching QSEIS06's native constant. This applies
frequency- and distance-dependent spatial smoothing; the spherical
examples use point sources. The 24-row half-space model and its
flat-Earth transformation also differ from the complete spherical model.
See [regional comparison limits](../guides/backend-comparison.md#qseis-at-the-same-regional-distances).
Reuse requires `--regional` and a library built with the current custom
source and the same selected observables and parameters. Rebuild libraries
created by the earlier type-2 regional example.

## Point-source control

Run a separate regional calculation to disable the Gaussian spatial
smoothing:

```console
python examples/qseis2025.py --regional --point-source
```

`--point-source` requires `--regional` and sets `source_radius_ratio=0`.
The default destination becomes
`examples/output/qseis2025-regional-point-source/`; the standard regional
library retains its 0.05 ratio. The time sampling, effective 64 s STF,
mechanism, moment, model and Earth flattening remain the same. Displacement
is saved as `disp.npz` with shape `(3, 3, 256)` and `disp.png`, together
with `summary.json` and the source-function records.

For a positive ratio, the native solver uses

```text
radius(f, r) = source_radius_ratio * min(sqrt(r**2 + (zs-zr)**2), Vp_source/(f+df))
multiplier(k) = exp(-(k*radius)**2/2)
```

Here `r` is epicentral distance, `zs-zr` is the source–receiver depth
separation in the backend coordinates, `Vp_source` is the source-layer
P-wave speed, `f` is frequency, `df` is the FFT frequency increment and
`k` is wavenumber. Use consistent length units. The backend applies any
selected Earth flattening before evaluating these quantities. The radius
changes with both receiver distance and frequency; it is not a fixed-radius
physical source disk. Setting the ratio to zero removes the Gaussian
multiplier and also changes the automatically estimated wavenumber cutoff.
This can substantially increase computation time.

```{figure} ../_static/examples/qseis2025-point-source.png
:alt: QSEIS2025 displacement with spatial Gaussian smoothing disabled.

Point-source control at 300, 600 and 900 km, with the same effective 64 s
STF, mechanism and 0–1020 s output window as the standard regional example.
```

The verified control changed only the spatial-source ratio among the
numerical parameters; it requested displacement only. Relative L2 differences
against the SPGRN2020 point-source calculation at 0.125 Hz were:

| Distance | QSEIS2025 ratio 0.05 | QSEIS2025 ratio 0 |
| --- | ---: | ---: |
| 300 km | 12.741% | 11.493% |
| 600 km | 19.574% | 17.076% |
| 900 km | 22.798% | 19.718% |

The comparison uses native origin-time coordinates, linearly interpolated
to a common 1 s grid from each SPGRN2020 start time (3, 40 or 78 s) through
500 s. Each value is the joint ENU norm of the difference divided by the
SPGRN2020 norm. No fitted amplitude scale or fitted time shift is applied.

```{figure} ../_static/examples/source-radius-comparison.png
:alt: Default and point-source QSEIS2025 waveforms compared with SPGRN2020.

Changing the spatial-source ratio reduces part of the discrepancy. A
substantial residual remains, so the control does not establish spatial
smoothing as the sole cause or demonstrate convergence of all settings.
```

The [recorded point-source run](../_static/examples/qseis2025-point-source.json)
took 647.470 s on the recorded Windows environment, versus approximately
247 s for the earlier default-ratio run. Its arrays were finite and the
native output contained no warnings or errors; runtimes depend on hardware
and concurrent work. See [validation](../validation.md) for run records.

Use `--regional --point-source --reuse` to reread a matching point-source
library. The example checks the native input's radius ratio and wavenumber
truncation tolerance as well as its other settings, and rejects a library
with incompatible spatial-source controls. QSEIS06 fixes the ratio at 0.05
in its native solver; its existing Python API has no equivalent flag.

## Parameters that control the calculation

`N_each_group` splits receiver distances into jobs; it is 3 in this example.
It controls job size and file layout, not the source model.
`wavenumber_sampling_rate`, `anti_alias`, `slowness_int_algorithm` and an
optional `slowness_window` control integration. QSEIS2025 additionally
exposes `eps_estimate_wavenumber` and `source_radius_ratio`.
Establish convergence when changing band, distance or depth.

`time_reduction_velo` is km/s; zero means no reduction. The raw
`free_surface` switch is the backend partial-solution switch, with 0
retaining free-surface effects. It is not the same Python boolean interface
as QSEIS06's `free_surface`.

The reader supports nearest or trilinear interpolation, ENZ rotation,
optional arrival adjustments and filtering/resampling. See
[reading](../guides/reading.md). An output family's source files must exist
before its derived rate/non-rate version can be requested.

## Outputs and limitations

Each depth pair contains distance-group directories with `grn.inp` and
basis files prefixed `ex`, `ss`, `ds` and `cl`. Components `e*` are strain
and `s*` are stress, despite swapped labels in some template comments.
The example retains ASCII and creates reader-compatible binary files.
Keep `green_lib_info.json` and travel-time tables with the data.

At `rotate=False`, six-component tensor output is a north-reference ENZ
basis, not the vector's RTZ convention. Use `rotate=True` and the documented
[component order](../conventions.md#six-component-tensors) for comparison.

The tutorial verifies executable output and reader behavior. A 30 km distance
grid and its integration defaults do not validate interpolation, high
frequencies, zero-distance behavior or interface stresses for your study.
