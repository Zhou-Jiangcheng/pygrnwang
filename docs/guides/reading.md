# Read, interpolate and process results

A Green's library stores basis responses on a depth/distance grid. A reader
selects or interpolates those responses, combines them for a focal mechanism,
rotates components and applies requested time/signal processing.

## Choose an observable

Enable a quantity while building the library before requesting it from a
reader. QSEIS06/SPGRN vector readers can derive displacement, velocity and
acceleration from the stored vector kernels. QSEIS2025 can also derive
rate/non-rate versions of an enabled tensor/rotation family. QSSP writes
separate observable files, so an unselected output is not synthesized by
the reader.

For the exact QSEIS2025 flag order and the direct strain/stress example see
[QSEIS2025](../backends/qseis2025.md). QSSP has a different eleven-flag order
and a documented [conversion limitation](troubleshooting.md#known-implementation-limitations).

The dynamic readers' default `only_seismograms=True` returns an array.
Use `False` to obtain the seven-element tuple described in
[scientific conventions](../conventions.md#time-origin-reduction-and-arrivals).
Do not unpack the default return as if it included metadata.
EDCMP returns static component vectors and has no time samples.

## Select a grid point

Nearest-neighbor selection is the default (`interpolate_type=0`). Main
QSEIS2025, SPGRN2012, SPGRN2020 and QSSP2020 readers accept
`interpolate_type=1` for linear interpolation across source depth, receiver
depth and distance. A one-element depth dimension simply contributes that
one node. The QSEIS06 main reader does not expose this argument.

Use sorted depth grids and query within the computed distance/depth extent.
Some interpolation branches clamp endpoints, but this is not uniform
validation across all readers and nearest-distance lookup can attempt to read
beyond the last file. Do not rely on out-of-range extrapolation.

The interpolation blends basis waveforms before source synthesis. It does
not automatically align each neighboring trace by phase. Wide grid spacing
can smear arrivals or create artifacts, especially across interfaces or
phase-branch changes. Compare an interpolated waveform with an explicitly
calculated intermediate point to select a useful grid spacing.

Metadata still identifies nearest nodes when interpolation is enabled;
it is not a report of every contributing node or interpolation weight.
Travel-time metadata likewise comes from the selected library table rather
than a weighted average of all neighboring times.

For many static queries, `seek_edcmp2_bulk` accepts matching arrays of query
parameters and returns `(n_queries, n_components)` using nearest grid nodes.
It needs the combined binary grids from `convert_pd2bin_edcmp2_all` and
does not expose the dynamic readers' `interpolate_type` option. Supply one
mechanism row per query and follow its API array-shape contract.

## Rotate components

Use `rotate=True` to obtain ENZ vectors or tensors in the documented ENZ
six-component order. `rotate=False` is suitable when a radial/transverse
vector is desired, but the unrotated tensor convention differs by backend.
See the [component tables](../conventions.md#six-component-tensors).

The lower-level `rotate_rtz_to_enz` takes an azimuth in degrees.
`rotate_symmetric_tensor_series` takes a rotation angle in radians and
expects `(n_samples, 6)`, whereas reader output is `(6, n_samples)`.
It implements `R.T @ tensor @ R`. Transpose explicitly and verify the axis
convention if using it outside a reader. A tensor cannot be rotated by
treating its six entries as independent vector channels.

## Time axis and optional arrival adjustment

For an unshifted trace, construct the time axis from its native start time
and returned `srate`. `np.arange(n_samples) / srate` alone is time since
the array start, not necessarily time since source origin.

With `before_p=20`, the library P onset is approximately 20 s from array
start, subject to integer-sample rounding. For a plot relative to library P,
use `np.arange(n_samples) / srate - 20`. With `shift=False`,
`first_p`/`first_s` remain `None`; when requested, use the onsets in
`tpts_table` instead.

`shift=True` uses separately recomputed travel times and stretches portions
of the waveform. Treat it as an approximation and verify the time axis
against the actual returned trace. It is not a simple constant phase shift.
Do not use this option when either required arrival is non-finite.

## Filtering and resampling

Reader filter arguments are:

| Argument | Meaning |
| --- | --- |
| `freq_band=None` or `[None, None]` | No filter |
| `[low, high]` | Bandpass corners in Hz |
| `[low, None]` | Highpass |
| `[None, high]` | Lowpass |
| `butter_order` | Butterworth order, default 4 |
| `zero_phase=False` | One-pass causal filtering |
| `zero_phase=True` | Forward/backward filtering |

Choose positive corners with `low < high` and the upper corner below the
native Nyquist frequency. The lower-level helper has special behavior at
or above Nyquist, so invalid corner choices should be corrected by the
caller rather than used as an implicit filter-selection mechanism.
Forward/backward filtering changes the effective response and needs a
trace long enough for padding.

Readers filter at the library sampling rate, then perform optional time
adjustments and resample to `srate`. Integer old/new rates use polyphase
resampling; other rates use the signal-processing helper. QSEIS/SPGRN
integration or differentiation for the output type occurs after resampling.
Do not manually differentiate velocity a second time when requesting
`output_type="acce"`.

The low-level `resample` and `taper` helpers operate on one trace at a time;
`filter_butter` can act on arrays along their last axis. Taper a working
copy if filtering a sharply truncated trace. Retain enough time before and
after the useful window to inspect edge effects.

Resampling cannot recover frequencies absent from the computed Green's
functions. When reducing sample rate, keep the useful band below the new
Nyquist frequency and validate amplitudes in that band.

## Export with meaning

The examples save arrays and component labels in `.npz`, plots in `.png`
and run metadata in `summary.json`. For a downstream format such as SAC,
MiniSEED or an application-specific tensor table, also export units, source
moment, coordinates, channel order, sample rate, absolute/relative start
time, backend version and processing parameters. A generic channel label
such as "Z" alone does not encode the physical quantity or its normalization.
