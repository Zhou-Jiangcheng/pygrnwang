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

The [fresh comparison dated 2026-09-15](../guides/backend-comparison.md)
uses 0.25 s sampling, a 2 Hz Nyquist limit, a 1.25 s effective source pulse
and receivers 1 km deep at 300, 600 and 900 km. It is separate from the
near-distance and 64 s regional tutorials below; these commands do not
produce the new comparison figures.

That comparison retains stock QSEIS06's frequency-dependent Gaussian
smoothing and adds an **isolated `rd2r=0` point-source control build**.
The control is not the released backend or an option in its Python API.
For the tested settings, stock QSEIS06 and ratio-0.05 QSEIS2025 have
bit-identical saved displacement arrays, as do the two point-source
controls. See the comparison for the limits of this result.

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
time axis. The regional QSEIS, SPGRN and QSSP examples share the 4 s
interval and 0.125 Hz Nyquist limit. QSEIS obtains its frequency range
from the sampling interval rather than a separate `max_frequency`
argument. With 1024 FFT samples, the final computed positive-frequency
bin is `511/4096 = 0.124755859375 Hz`; the Nyquist bin is zero.

The regional source uses `wavelet_type=0, wavelet_duration=16` with
1024 custom moment-rate samples spanning 0–64 s. The physical target is
a normalized squared half-sinusoid. Before writing the input, the
example multiplies its samples by `exp(2*pi*fi*t)` to compensate for
QSEIS's numerical damping convention; `fi` is negative. The written
samples are not renormalized: their area is approximately 0.9647094,
while the effective physical pulse has unit area and a 32 s centroid.
The [STF verification](../guides/backend-comparison-64s.md#matching-the-effective-source-time-function)
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

The regional examples share the effective 64 s moment-rate pulse,
strike/dip/rake 30°/45°/90°, moment `10^15 N m` and azimuth 30°.
Matching those settings and the frequency band does not make the
half-space and spherical calculations identical.

QSEIS06 also applies frequency- and distance-dependent Gaussian spatial
smoothing. The native solver fixes its radius ratio at 0.05; the current
Python API cannot turn it off. The radius at each receiver and frequency
is 0.05 times the smaller of the source–receiver separation and
`Vp_source/(f + df)`, using the backend's model coordinates and source-layer
P-wave speed. Its wavenumber multiplier is `exp(-(k*radius)**2/2)`.
This differs from the spherical examples' point sources even when their
effective time functions agree. The default QSEIS2025 regional example
retains the same 0.05 ratio for comparison with QSEIS06.

In the archived 64 s comparison, the
[QSEIS2025 point-source control](qseis2025.md#point-source-control)
turns off this spatial smoothing. It reduces the discrepancy against
SPGRN2020, but the measured relative differences remain 11.493%, 17.076%
and 19.718% at 300, 600 and 900 km. Spatial smoothing therefore explains
only part of the discrepancy. See the
[regional comparison limits](../guides/backend-comparison-64s.md#qseis-at-the-same-regional-distances)
for the remaining differences in geometry, model and numerical treatment.

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
