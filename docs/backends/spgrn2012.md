# SPGRN2012: reduced-time spherical waveforms

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

The native start time is `t0 + distance / v0`. Here `t0=-40` s and
`v0=10` km/s, so the three nominal start times are -10, 20 and 50 s.
This is not the P-relative convention of SPGRN2020.

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

The 64 s source and coarse distance grid demonstrate a long-period workflow.
They do not establish convergence for high-frequency regional phases or
near-field static response. Compare unshifted grid-point traces first,
then validate interpolation or arrival-based waveform adjustments for
your application.
