# QSEIS06 (deprecated): vector waveforms

```{admonition} Deprecated
:class: warning

The QSEIS06 backend is deprecated in pygrnwang. Use
[QSEIS2025](qseis2025.md) for new calculations. This tutorial and its
interfaces remain available for existing workflows. Migration requires
building a QSEIS2025 library and validating the numerical settings, source
time function and output conventions.
```

QSEIS06 provides the established layered half-space vector workflow. The
main reader returns displacement, velocity or acceleration. It remains
useful for reproducing QSEIS06 calculations and comparing an established
library with the newer direct-observable workflow.

## Complete calculation

```console
python examples/qseis06.py
```

The default near-distance mode builds one 10 km source depth and a
surface receiver at 30, 60 and 90 km, using a 0.5 s interval and 127.5 s
window. Like the default QSEIS2025 example,
it uses 24 numeric model rows and no flat-Earth transformation. The source
has strike/dip/rake 30°/45°/90° and moment `10^15 N m`.

```{literalinclude} ../../examples/qseis06.py
:language: python
:caption: Complete QSEIS06 tutorial
```

```{figure} ../_static/examples/qseis06.png
:alt: QSEIS06 east, north and up displacement for three distances.

Displacement in metres, plotted against seconds since source origin.
```

Expect `disp.npz` with waveform shape `(3, 3, 256)`, `disp.png`,
`summary.json` and a library beneath `examples/output/qseis06/`.
`--reuse` rereads that library.

## Regional waveforms at 300, 600 and 900 km

Run the same script with the regional option:

```console
python examples/qseis06.py --regional
```

This selects 300, 600 and 900 km while retaining a 10 km source and
surface receiver. It enables `flat_earth_transform=True` and keeps the
same 24 numeric model rows. A 4 s interval and 4092 s native window give
1024 library samples. After synthesis, the script saves 0–1020 s
inclusive, or 256 samples, under `examples/output/qseis06-regional/`.
The resulting `disp.npz` has shape `(3, 3, 256)` and a source-origin
time axis.

The regional source uses `wavelet_type=0, wavelet_duration=16` with
1024 custom moment-rate samples spanning 0–64 s. The physical target is
a normalized squared half-sinusoid. Before writing the input, the
example multiplies its samples by `exp(2*pi*fi*t)` to compensate for
QSEIS's numerical damping convention; `fi` is negative. The written
samples are not renormalized: their area is approximately 0.9647094,
while the effective physical pulse has unit area and a 32 s centroid.
The [STF verification](../guides/backend-comparison.md#matching-the-effective-source-time-function)
checks the actual pulse and its spectrum against the target.

A custom type-0 pulse does not trigger the ordinary reader's automatic
conversion between velocity and displacement. This example explicitly
reads `output_type="velo"`, integrates once using `cumsum * dt`, and
then saves displacement. The default near-distance example retains
type 2 with four 0.5 s samples; its behavior is unchanged.

```{figure} ../_static/examples/qseis06-regional.png
:alt: QSEIS06 regional displacement at 300, 600 and 900 km.

Regional displacement in metres, saved from 0 to 1020 s since source
origin. The native 4092 s library window is retained.
```

Matching the effective source pulse does not make this truncated
half-space model identical to the complete spherical model. The
flat-Earth transformation does not restore the omitted deep structure,
and the QSEIS numerical frequency range extends above the spherical
examples' 0.0625 Hz cutoff. See
[regional comparison limits](../guides/backend-comparison.md#qseis-at-the-same-regional-distances).
Use `--regional --reuse` only for a completed library with the current
custom source and matching regional parameters. A library from the
earlier type-2 regional example must be rebuilt.

## Source, sampling and boundary parameters

`wavelet_type=2` selects the tapered Heaviside. In the default
near-distance mode, `wavelet_duration=4` is four samples (2 s).
Type 1 instead stores a
velocity-like kernel; the reader integrates/differentiates as required.
`output_type` is `disp`, `velo` or `acce`.

`time_reduction_velo` is km/s and zero disables reduction.
`free_surface=True` retains free-surface effects; the writer translates
this boolean into the backend switch. `flat_earth_transform` controls
the geometric transformation; preserve it when reproducing a library.

`N_each_group` is the number of receiver distances per backend job.
`wavenumber_sampling_rate`, `slowness_int_algorithm`, `slowness_window`
and `anti_alias` affect integration and must be checked for convergence.
The ordinary reader chooses nearest source depth, receiver depth and
distance; it does not provide `interpolate_type`.

## Separate finite-difference tensor workflow

`pre_process_qseis06_strain_rate` creates extra distance/depth samples for
the derivative readers. Its stencil uses `diff_accu_order`, radial
increment ratio `k_dr` and depth increment `dz` in km.
An ordinary vector library does not contain those additional samples.

`seek_qseis06_strain_rate_diff` and
`seek_qseis06_stress_rate_diff` are separate advanced interfaces.
The latter uses elastic moduli for tensor conversion. They require
stencil convergence, layer/interface treatment and rotation checks
appropriate to the intended calculation. Keep derivative ASCII outputs
with `convert_pd2bin=False, remove_pd=False` because the current derivative
binary detection/call path is inconsistent; see
[known limitations](../guides/troubleshooting.md#known-implementation-limitations).
The [QSEIS2025 tensor
example](qseis2025.md#add-tensor-outputs) is the tested introductory route
to direct strain and stress.

## Stored data and interpretation

A source/receiver directory contains distance groups and basis response
files. Normal vector groups end in `_0`; derivative libraries can have
additional offset groups. The tutorial keeps ASCII while preparing the
binary representation for the reader. Retain metadata and matching
travel-time tables.

`rotate=True` returns E/N/up vectors. `only_seismograms=False` includes
nearest grid metadata; new `first_p`/`first_s` are still `None` unless
`shift=True`. Numerical agreement with another backend requires consistent
model, source wavelet, amplitude and time reference.
