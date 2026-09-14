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
python examples/spgrn2012.py --output-dir examples/output/spgrn2012-matched-band
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

The command above saves results beneath
`examples/output/spgrn2012-matched-band/`. `disp.npz` has shape
`(3 distances, 3 ENU components, 256 samples)` in metres; `disp.png` and
`summary.json` describe the exported result. `library/` contains the
native spectral and velocity files. Three additional archives make the
source calculation inspectable:

- `velocity-impulse.npz`: the complete native impulse response, read and
  rotated to ENU and scaled by the stated moment, in m/s;
- `velocity-matched.npz`: the complete velocity after forward source
  convolution, in m/s;
- `source_time_function.npz` and `.json`: the physical pulse, its analytic
  complex-frequency transform, validated native settings and archive hashes.

Both velocity archives have shape `(3, 3, 1024)` and retain the original
per-distance time axes. The script's default output directory remains
`examples/output/spgrn2012/`; the explicit directory above keeps this
calculation separate from earlier tutorial runs.

## Spectra, output window and source

Both the spectral window and the **native velocity output window** are
4092 s, sampled at 4 s. The native FFT has 1024 samples and period 4096 s.
The script retains all 1024 samples while applying the physical source,
then integrates once and exports the first 256 samples spanning 1020 s.
Convolving only the cropped displacement cannot reproduce this operation.

`max_frequency=0.125` Hz equals the sampling Nyquist frequency.
The actual spectrum header must contain `nfcut=512`, with frequency spacing
`1/4096` Hz. The solver explicitly zeroes the Nyquist bin, so its highest
computed frequency is `511/4096 = 0.124755859375` Hz.
`max_slowness=0` requests the full wavefield. `cal_sph=1` and `cal_tor=1`
include P-SV and SH contributions, and `source_radius=0` selects a point
source. `gravity_fc`/`gravity_harmonic` and `physical_dispersion` are zero;
the native additional Butterworth filter is disabled.

The physical source shared with the other regional examples is the
unit-area moment-rate pulse

```{math}
r(t)=\frac{2}{64}\sin^2\left(\frac{\pi t}{64}\right),\qquad 0\leq t\leq64\ \mathrm{s},
```

with zero rate outside that interval and centroid 32 s. This tutorial
sets the **native** `source_duration=0`, whose spectrum is unity, and
applies the target 64 s source in the example script. The library metadata
therefore correctly records zero native duration; the summary separately
records `native_source_duration_s=0` and `effective_source_duration_s=64`.
No solver kernel or public reader behavior is changed.

The native imaginary frequency is
`fi = log(0.01) / (2*pi*4096)`. The helper restores the numerical damping
of the complete impulse velocity, transforms it, multiplies by the exact
transform of `r(t)` at `f + i*fi`, transforms back and removes damping.
This is forward convolution: there is no division by an existing source
spectrum, fitted amplitude or time shift. The damped DC coefficient is
about 0.9647432; it must not be renormalized because the **physical** rate
already integrates to one. Independent quadrature checks the analytic
transform across the full 0–0.125 Hz band.

```{literalinclude} ../../examples/spherical_source_time_function.py
:language: python
:pyobject: apply_spgrn2012_stf
:caption: Full-period forward source convolution before displacement integration
```

`cal_gf=1` requests spectral calculation. `delta_dist_range` gives the
smallest/largest distance increments in km; equal values request a uniform
grid. The script reads the actual `dist_list` from `green_lib_info.json`.

Fortran rounds `t0 + distance / v0` to the nearest integer second for
the native start time. Here `t0=-40` s and `v0=10` km/s, giving starts
of -10, 20 and 50 s. The saved displacement axes therefore cover
-10–1010, 20–1040 and 50–1070 s relative to source origin. SPGRN2020
uses a different start rule based on P onset; compare their traces only
over the shared origin-time interval.

Use `--reuse` with the same explicit output directory only after a
successful run. The script checks native input records, all spectral
headers, complete velocity blocks, the physical pulse and archive hashes.
Older libraries with a native 64 s source, a 0.0625 Hz cutoff, a slowness
limit or a cropped native output are rejected and require a fresh build.

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

The raw basis library represents velocity. In the general reader,
`output_type="disp"` integrates it and `"acce"` differentiates it.
This example explicitly reads `"velo"`, applies the physical source above,
and uses `cumsum * dt` exactly once to obtain displacement. It scales a
unit mechanism by `10^15 N m` and returns E/N/up components. Nearest and
trilinear waveform interpolation remain available in the reader.

Earlier tutorial comparisons reported approximately 4% differences from
SPGRN2020. Those results used a different frequency cutoff and SPGRN2012's
native real-frequency wavelet, whose effective area became about 3.67%
larger after damping correction. They are historical results, not a
characterization of the source-matched example above. See the
[controlled backend comparison](../guides/backend-comparison.md) for the
current calculation, measured differences and remaining model/sampling
limitations. Matching the physical source and requested frequency band
does not establish exact equivalence between the solvers.

The 64 s source and coarse distance grid do not establish convergence for
high-frequency regional phases or other source/receiver geometries.
Compare unshifted grid-point traces first, then validate interpolation
or arrival-based waveform adjustments for your application.
