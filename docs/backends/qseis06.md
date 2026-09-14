# QSEIS06: vector waveforms

QSEIS06 provides the established layered half-space vector workflow. The
main reader returns displacement, velocity or acceleration. It remains
useful for reproducing QSEIS06 calculations and comparing an established
library with the newer direct-observable workflow.

## Complete calculation

```console
python examples/qseis06.py
```

This builds one 10 km source depth and a surface receiver at 30, 60 and
90 km, using a 0.5 s interval and 127.5 s window. Like the QSEIS2025 example,
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

## Source, sampling and boundary parameters

`wavelet_type=2` selects the tapered Heaviside, and
`wavelet_duration=4` is four samples (2 s here). Type 1 instead stores a
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
