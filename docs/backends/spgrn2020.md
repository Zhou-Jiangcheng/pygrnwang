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
:alt: SPGRN2020 displacement plotted in seconds since source origin.

Displacement in metres versus seconds since source origin. Each trace's
start time comes from its native binary record header.
```

Outputs are under `examples/output/spgrn2020/`. `disp.npz` contains
three distances, E/N/up channels, 256 samples per trace and the plotted
time coordinates; the summary records the run.

## Parameter choices

`spec_time_window=4092` s exceeds `time_window=1020` s. The 4 s sample
interval, 0.125 Hz cutoff and 64 s source duration define this small
long-period calculation. Source duration is in seconds. The example
selects spheroidal and toroidal modes, disables the configured
self-gravitation range and uses `cal_gf=1` for new spectra.

`max_slowness=0` selects SPGRN2020's existing full-wavefield branch,
which chooses a model-dependent slowness limit and a larger low-frequency
harmonic baseline. It does not restrict the calculation to zero slowness.
The previous positive cutoff of 0.3 s/km underestimated the required
low-frequency content at 300 km in this example. See the
[controlled comparison](../guides/backend-comparison.md) for the evidence
and the scope of the revised setting.

`green_before_p=40` requests a window beginning approximately 40 s
before direct P. The wrapper writes its negative as the Fortran start-time
offset; Fortran rounds the resulting start time to the nearest integer
second. The example reads that value from each native record header and
plots `t_start + np.arange(n_samples) * 4` seconds since source origin.
The three starts are 3, 40 and 78 s. Simply adding fractional P-table
onsets to a -40 s plotting axis would not reproduce those stored starts.

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

The storage window remains referenced to P even though the example now
plots source-origin time. Compare backends using the saved time coordinates,
matching source, model, units and processing. Raw sample indices do not
identify the same physical time across these libraries.

The revised example was checked against QSSP2020 at these three distances
and this frequency band. Agreement is a numerical cross-check, not an
absolute reference solution or validation for other depths, distances or
frequencies. Native Fortran binary metadata must be read with the matching
package reader; do not treat it as headerless travel-time arrays.
