# SPGRN2020: windows referenced to P

SPGRN2020 builds spherical-Earth waveform libraries whose stored time
windows begin a chosen interval before direct P. Its native arrival
tables include onset, takeoff angle and slowness.

## Complete calculation

```console
python examples/spgrn2020.py
```

This builds a full-Earth long-period library for a 10 km source, surface
receiver and approximately 300, 600 and 900 km. It obtains the actual
distance list from the completed backend metadata.

```{literalinclude} ../../examples/spgrn2020.py
:language: python
:caption: Complete SPGRN2020 tutorial
```

```{figure} ../_static/examples/spgrn2020.png
:alt: SPGRN2020 displacement plotted relative to the library P arrival.

Displacement in metres. Zero on this figure is the library P arrival,
not source origin.
```

Outputs are under `examples/output/spgrn2020/`. `disp.npz` contains
three distances, E/N/up channels, 256 samples per trace and the plotted
time coordinates; the summary records the run.

## Parameter choices

`spec_time_window=4092` s exceeds `time_window=1020` s. The 4 s sample
interval, 0.0625 Hz cutoff and 64 s source duration define this small
long-period calculation. Source duration is in seconds. The example
selects spheroidal and toroidal modes, disables the configured
self-gravitation range and uses `cal_gf=1` for new spectra.

`green_before_p=40` means the stored trace begins 40 s before direct P.
The wrapper writes its negative as the Fortran start-time offset.
With no reader adjustment, the example therefore plots
`-40 + np.arange(n_samples) * 4` seconds relative to P.

`dist_range` and `delta_dist_range` are in km. The backend can choose
distance-dependent spacing; use the generated `dist_list` rather than
assuming every request was stored exactly. Recompute spectra if the
model, spectral sampling, frequency/slowness cutoffs or mode settings
change.

## Outputs and reading

`GreenSpec/` holds spectra and `GreenFunc/` holds the basis waveforms,
`GreenInfo*.dat` and native `tptable.dat`/`tstable.dat`. The metadata
update after calculation records `dist_list` and `samples_num`.

`seek_spgrn2020` returns displacement, velocity or acceleration, with
nearest or trilinear interpolation. It synthesizes from ten elementary
velocity traces and returns E/N/up by default. The travel-time dictionary
also exposes P/S takeoff angles in degrees and native slowness in **s/m**.
Multiply slowness by 1000 for s/km. It is neither the input cutoff's
s/km convention nor TauP ray parameter in s/radian.

The module also contains a cache and a precompute/fast-reader path for
repeated queries. Use the main reader as the reference while validating
a repeated-query optimization, including cache lifetime and timing options.

## Limits and comparison

A P-relative plot makes arrivals easy to compare, but source-origin
timing requires adding each trace's P onset. SPGRN2012 and SPGRN2020
can produce different-looking raw arrays simply because their start
times differ. Match source, model, units and time origins before comparing.

The tutorial's coarse grid and long source are not high-frequency or
near-field convergence checks. Native Fortran binary metadata must be
read with the matching package reader; do not treat it as headerless
travel-time arrays.
