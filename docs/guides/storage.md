# Library layout and data formats

Keep each backend and numerical configuration in a separate library root.
`green_lib_info.json` describes the grid, sample interval, time window and
other parameters required by the corresponding reader. It is part of the
data, not a disposable cache.

## Backend directories

| Backend | Important paths beneath its library root |
| --- | --- |
| QSEIS06/2025 | `<source:.2f>/<receiver:.2f>/<group>_<offset>/grn.inp` and basis output files |
| SPGRN2012/2020 | `GreenSpec/<source:.2f>/<receiver:.2f>/` and `GreenFunc/<source:.2f>/<receiver:.2f>/` |
| QSSP2020 | `GreenSpec/<source:.2f>/<receiver:.2f>/spec.inp`; `GreenFunc/<source:.2f>/<receiver:.2f>/<tensor-source>/` |
| EDGRN | `edgrn2/<receiver:.2f>/`, including `edgrn.ss`, `edgrn.ds` and `edgrn.cl` |
| EDCMP | `edcmp2/<source:.2f>/<receiver:.2f>/<basis-index>/`, with five basis indices |
| All preprocessed libraries | `green_lib_info.json`, task-group `.pkl` files and the model/`noQ.nd` |

Depth directory names have two decimal places. Choose depths that remain
distinct after this formatting; otherwise separate requested depths can
refer to the same directory. QSEIS distances are grouped to limit one
backend job's receiver count. Its ordinary group suffix is `_0`; derivative
libraries may have additional offset suffixes.

The actual sampled distance grid can extend beyond an input upper bound
when the interval does not divide the range exactly. Query the computed
grid recorded in metadata rather than inferring a different grid with an
unmatched `linspace`. SPGRN may use variable distance spacing and stores
its actual `dist_list` after the backend finishes.

## ASCII and binary

QSEIS ASCII files use source-type prefixes such as `ex`, `ss`, `ds`, `cl`
and observable extensions. Conversion writes `grn_<component>.bin`.
These binary arrays use `float32`, grouped by distance, basis source and
time as expected by the matching reader.

QSSP writes one file per selected observable component and six spherical
moment-source bases, `mrr, mtt, mpp, mrt, mrp, mtp`. Conversion removes
the time column and writes receiver-major `float32` data. Native SPGRN
files contain Fortran records and metadata; they are not plain arrays of
three equally sized displacement channels.

EDCMP `hs.disp`, `hs.strain`, `hs.stress` and `hs.tilt` are static output
tables. Conversion writes `float32` component arrays and, where the bulk
converter is used, combined grids for efficient bulk lookup. The main
reader can read the per-basis data.

Use the matching package converter and reader. Raw `numpy.fromfile`
without the exact layout can silently produce plausible but incorrectly
ordered data. Binary files are not self-describing, and native Fortran
record structure/byte order can depend on how a backend was built. Keep
the executable/package version and metadata with archived data.

## Conversion and cleanup

QSEIS creation functions generally convert after computing. QSSP's
sequential/parallel paths default to `convert_pd2bin=True, remove_pd=True`,
so ASCII files can disappear after conversion. During a first run, retaining
ASCII helps inspect units, headers and failures. Disable removal with
`remove_pd=False` where that argument exists.

EDCMP sequential creation does not perform the bulk conversion automatically;
call `convert_pd2bin_edcmp2_all` explicitly if it is required. Its tutorial
shows the complete sequence. QSSP has an
[observable-selection conversion issue](troubleshooting.md#known-implementation-limitations)
for strain/rotation, so retain original output when exploring those families.

Before removing ASCII, confirm that the binary reader reproduces the expected
array shape, finite values and representative traces. Do not move only the
`.bin` files while leaving the metadata or required travel-time tables behind.

## Travel-time files and trusted libraries

Python-created `tp_table.bin`/`ts_table.bin` store one `float32` time per
distance. SPGRN2020 `tptable.dat`/`tstable.dat` instead contain native
Fortran onset/takeoff/slowness records; use their designated reader.

The task-group `.pkl` files are internal serialized Python task lists.
Only run preprocessing/creation against libraries you created or trust,
because loading pickle can execute Python code. For sharing scientific
results, include neutral arrays, model text and JSON metadata alongside
the full library when practical.

## Tutorial artifacts

Each script uses `examples/output/<backend>` by default, with a library
subdirectory and saved observable arrays, PNG figures and `summary.json`.
Use `--output-dir` to select a separate location and `--reuse` to read a
completed tutorial run. Reuse writes `summary-reuse.json` and preserves
the original construction report. Figures show the saved arrays; the summary
reports shape, finite-value checks, elapsed time and output size.
