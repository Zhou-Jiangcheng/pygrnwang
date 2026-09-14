# SPGRN2012 (deprecated): reduced-time spherical waveforms

```{admonition} Deprecated
:class: warning

The SPGRN2012 backend is deprecated in pygrnwang. Use
[SPGRN2020](spgrn2020.md) for new calculations. This tutorial and its
interfaces remain available for existing workflows. Migration requires
building a SPGRN2020 library and validating its sampling and time origin;
SPGRN2020 uses windows referenced to P instead of SPGRN2012's reduction
rule.
```

SPGRN2012 computes waveforms for a radially layered spherical Earth. Its
time windows use a reduction offset and velocity, which makes it useful
for reproducing libraries based on `t0 + distance / v0`.

## Complete calculation

```console
python examples/spgrn2012.py
```

The example uses the full Earth model, a 10 km source, surface receivers
at approximately 300, 600 and 900 km, and one serial backend job. The
actual grid is read from backend metadata after calculation.

```{literalinclude} ../../examples/spgrn2012.py
:language: python
:caption: Complete SPGRN2012 tutorial, including travel-time tables
```

```{figure} ../_static/examples/spgrn2012.png
:alt: SPGRN2012 displacement at three distances with reduced start times.

Displacement in metres. Each trace uses its own reduced start time to
show time since source origin.
```

The output defaults to `examples/output/spgrn2012/`, with `disp.npz`,
`disp.png`, `summary.json` and `library/`. The waveform array has three
distances, three components and 256 samples.

## Spectra, output window and source

The example's spectral window is 4092 s and its output window is 1020 s,
both sampled at 4 s. These are distinct settings: a longer spectral window
helps separate repeated signals in the inverse transform while keeping
the useful saved waveform short.

`max_frequency=0.0625` Hz selects a long-period calculation;
`source_duration=64` s is a squared half-sinusoid duration.
`max_slowness` is s/km. `cal_sph` and `cal_tor` select P-SV and SH
contributions. `gravity_fc`/`gravity_harmonic` control the range of
self-gravitation; both are zero in this example. `physical_dispersion=0`
is also an explicit tutorial choice.

`cal_gf=1` requests spectral calculation. Reuse spectra with `0` only
when the model and spectral parameters are unchanged.
`delta_dist_range` gives smallest/largest distance increments in km;
equal values request a uniform grid. Use the actual `dist_list` in
`green_lib_info.json` for reading.

Fortran rounds `t0 + distance / v0` to the nearest integer second for
the native start time. Here `t0=-40` s and `v0=10` km/s, giving starts
of -10, 20 and 50 s. SPGRN2020 uses a different rule based on P onset.
Both example figures nevertheless express their time axes relative to
source origin, allowing comparisons over their common time interval.

## Travel-time tables are part of the workflow

The sequential builder updates library metadata but does not create the
Python travel-time tables consumed by `seek_spgrn2012`.
The script explicitly runs `create_tpts_table` after construction, using
the generated `noQ.nd` and backend's actual distance grid.

Keep the resulting `tp_table.bin`/`ts_table.bin` beneath the depth pair in
`GreenFunc/`. The reader loads them even when `shift=False`.
A successful Fortran run alone is therefore insufficient for a complete
Python-readable SPGRN2012 library.

## Reading and limits

The raw basis library represents velocity. `output_type="disp"` integrates
it; `"acce"` differentiates it. The example scales a unit mechanism by
`10^15 N m` and returns E/N/up displacement. Nearest and trilinear
waveform interpolation are available.

The [controlled spherical-backend comparison](../guides/backend-comparison.md)
found approximately 4% displacement differences from SPGRN2020's
full-wavefield result over the shared interval ending at 500 s.
SPGRN2012's older wavelet implementation evaluates the source spectrum
at real frequency, while SPGRN2020 and QSSP2020 include the imaginary
frequency used for numerical damping. With this example's parameters,
the older implementation produces an effective source-pulse area about
3.67% larger after damping correction. This is consistent with much of
the amplitude difference; it does not explain every residual difference.

The 64 s source and coarse distance grid do not establish convergence for
high-frequency regional phases or other source/receiver geometries.
Compare unshifted grid-point traces first, then validate interpolation
or arrival-based waveform adjustments for your application.
