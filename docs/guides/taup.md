# Travel times with TauP

Travel-time queries use a Java subprocess when both `java` and `javac` are
on `PATH` and a bundled `TauP.jar` can be found. Otherwise the public
first-arrival functions use ObsPy. **JPype is not required.** Importing
`pygrnwang.pytaup` detects availability but does not start Java or compile
the bridge.

The package bundles TauP 2.6.1; commands or requirements from the latest
[TauP documentation](https://www.seis.sc.edu/TauP/) need not match this
bundled version. The package wrapper is the interface described here.

## Setup and selection

Install a JDK to use Java mode: a runtime containing `java` alone is
insufficient because the bridge must first be compiled with `javac`.
Use `java -version` and `javac -version` in the same environment where
Python will run. After changing `PATH`, start a new Python process because
backend selection occurs at module import.

Resource lookup checks the package's `exec/TauP.jar` first, then
`sys.exec_prefix/Scripts/TauP.jar` on Windows or
`sys.exec_prefix/bin/TauP.jar` elsewhere. Standard wheels contain the package
copy and install the environment copy.

The bridge is compiled lazily into a temporary directory once per Python
process. Each individual query starts a Java subprocess. Java compilation,
model-loading or query errors are reported as exceptions; automatic fallback
means selection when Java/JAR is unavailable, not retrying every Java error
with ObsPy.

## Single-location queries

`cal_first_p` and `cal_first_s` return seconds; `cal_first_p_s` returns
`(first_p, first_s)`. Their source/receiver depths and distance are in km.
The default model is `ak135`. For a custom model pass a `.nd` filename.
When the receiver is deeper than the source, these first-arrival helpers
swap the two depths to use reciprocity.

"First P" and "first S" refer to the earliest arrival in the package's
explicit phase lists:

| Family | Queried names |
| --- | --- |
| P | `p, P, pP, Pg, Pn, Pdiff, PKP` |
| S | `s, S, sS, pS, Sg, Sn, Sdiff, SKS` |

These are not a promise that every conceivable P-like or S-like branch is
included. No arrival yields `NaN`. Check finiteness before aligning a
waveform or converting a time to a sample index.

`taup_time_java` is the Java-specific interface for a user-provided
`phases_list`. It returns a dictionary with parallel lists
`phase`, `puristphase`, `time` and `rayparameter`. Empty lists mean no
arrival. Ray parameter is **seconds per radian**, preserving TauP's
underlying convention. Convert to seconds per degree by multiplying by
`pi / 180`; it is not already a s/km slowness.

## Custom models

`taup_create_npz_file(nd_file)` prepares a model for the selected backend.
In Java mode it returns the `.nd` path unchanged. In ObsPy mode it builds
and returns an adjacent `.npz`. Use its returned path when scripting both
backends. Java's direct file reader is not an ObsPy `.npz` reader.

The four-column travel-time model uses depth, P/S velocity and density.
The preprocessing helpers remove attenuation from the six-column wave
model. Prepare the custom model in the main process before creating workers
to avoid simultaneous model builds. Keep a writable model directory for
ObsPy and rebuild the model if the `.nd` contents change.

## Distance tables

`create_tpts_table` writes `tp_table.bin` and `ts_table.bin` in
`<path_green>/<source_depth:.2f>/<receiver_depth:.2f>/`. Each is a
headerless `float32` array in exactly the order of `dist_km_list`; its unit
is seconds. The function writes files and returns no array.

Java mode sends the complete distance list to a single Java process, reusing
one model during that query. It is preferable to repeated individual
queries on a distance grid. In ObsPy mode, `max_workers` controls the
process pool; lists shorter than 50 distances run serially.

`check_finished=True` skips calculation only when both output files exist.
It does not validate their length, distances or model identity. Reuse only
when all inputs match and verify the stored grid. The supplied tutorials
perform the relevant table preparation as part of a complete calculation.

SPGRN2020 also writes native Fortran `tptable.dat` and `tstable.dat`
including takeoff angles in degrees and slowness in s/m. Multiply this
native slowness by 1000 for s/km. Those files have their own reader and
are not interchangeable with the Python `tp_table.bin`/`ts_table.bin` arrays.
