# Troubleshooting and known limitations

Start by running the smallest tutorial for the relevant backend in a fresh
output directory. Retain the generated input, model, metadata and executable
output when reporting a failure.

## Installation or executable failures

A compatible wheel avoids a local Fortran build. A source installation
requires gfortran; compiled executables also need their platform runtime
libraries. A successful `import pygrnwang` alone does not demonstrate that
a backend executable can run.

On Linux/macOS, the current bulk calculation helper looks for
`<environment>/bin/<solver>.bin`, while standard wheels keep native
executables under `pygrnwang/exec/` and provide entry points without
the `.bin` suffix. For current bulk calculation tutorials, use the
documented source/editable installation, which places the expected files
in the environment. The direct wrapper command can find a package binary,
but this does not fix the bulk helper's separate lookup path.

Check the active Python environment, package version and backend command.
On Windows use an activated environment or `conda run -n ENV python ...`;
launching an environment's `python.exe` directly may omit DLL directories.
Use short, local, writable paths if a legacy Fortran path-length or input
filename issue is suspected. Inspect `grn.inp` or `spec.inp` for absolute
paths to the intended files.

Some backend calls change working directory. Resolve model/output paths
before the first call and avoid assuming the caller's working directory
remains unchanged.

## Java or model failures

Both `java` and `javac` are required for Java mode. Missing tools or JAR
select ObsPy when `pytaup` is imported. A selected Java backend that fails
to compile or load a model raises an error; it is not silently retried
with ObsPy. Start a new Python process after changing `PATH`.

For custom models, check numeric column count, velocity discontinuities,
boundary names, model extent and model-directory write access. ObsPy may
need to create an adjacent `.npz`. Match the travel-time model to the
wave-propagation model and regenerate tables after a change.

## Missing files, reshape errors or zero results

| Symptom | Checks |
| --- | --- |
| Missing `green_lib_info.json` | Run preprocessing for this backend in the correct root |
| Missing requested observable | Confirm its build-time output flag was enabled |
| Missing QSSP spectra | First run requires `cal_spec=True` |
| Missing SPGRN2012 `tp_table.bin` | Sequential/MPI paths need explicit Python table preparation |
| Unexpected array size/reshape error | Check native sampling, actual distance grid, matching metadata and completed conversion |
| No signal in selected time window | Check native start time, phase arrival, source mechanism, output band and component |
| Non-finite or implausible result | Inspect executable logs, model material, source normalization and numerical convergence |

A completed process or `.finished` marker is not a substitute for finite
values and scientifically plausible output. A component may be physically
zero for a particular mechanism/azimuth; inspect all components and a
second mechanism before concluding a solver failed.

## Arrival and waveform alignment

`first_p`/`first_s` equal `None` when a reader did not request a new
travel-time calculation (`shift=False`). This is expected. A TauP first
arrival of `NaN` means the phase query did not supply a usable arrival.

A trace beginning at index zero does not necessarily begin at source
origin: QSEIS, SPGRN and QSSP have different time-reduction conventions.
Check [time origin](../conventions.md#time-origin-reduction-and-arrivals).
Avoid combining `before_p` and `pad_zeros`. Validate `shift=True` against
a direct computation before using its piecewise resampling for scientific
measurements.

Short or sharply truncated traces can show filter-padding errors or edge
ringing. Keep the frequency band below Nyquist, retain time margins, and
inspect the unfiltered trace and source-time function.

## Known implementation limitations

These are observations from the current source audit, not changes to the
numerical code made by the documentation project.

**QSSP selected strain/rotation conversion.** The Fortran output flags
(zero-based) 3/4 mean strain/strain rate and 7/8 mean rotation/rotation rate.
The Python `create_qssp2020.output_type_list` instead places rotation at
3/4 and strain at 7/8. The automatic bulk converter uses that Python list
to decide which files to open. Selecting only one of these families can
therefore make conversion seek a file that was not generated. Displacement,
velocity, acceleration and stress/stress-rate entries agree between both
lists. The introductory QSSP tutorial exercises displacement. Retain native
ASCII and inspect enabled output files when investigating the affected
families; automatic conversion of arbitrary selections is not validated.

**QSEIS06 derivative libraries.** The finite-difference reader requires
the extra spatial samples and metadata produced by
`pre_process_qseis06_strain_rate`. An ordinary `pre_process_qseis06`
library cannot replace it. Retain ASCII when executing derivative jobs:
the binary detector looks for `grn_tz.npy` and its binary-reader call
omits the required `sampling_num` argument. Use
`convert_pd2bin=False, remove_pd=False` for that path. The derivative
tensor rotation uses `-az_deg` while the QSEIS2025 direct reader uses
`+az_deg`; compare signs against an independently checked calculation
rather than assuming equivalence. Prefer the tested direct QSEIS2025
strain/stress tutorial for an introductory tensor workflow; derivative
convergence and sign checks require their own calculation.

**QSEIS2025 template comments.** Some legacy input-template descriptions
swap the strain/stress labels on `e*`/`s*` filenames. The current Fortran
calculation and Python reader use `e*` for strain and `s*` for stress.
Follow the observable table in the tutorial and the actual output routine.

**SPGRN2012 sequential/MPI travel-time tables.** These creation paths update
library metadata but do not create the Python P/S tables that the reader
loads on every query. The complete tutorial explicitly creates those
tables using the actual computed distance grid.

**EDGRN source-depth sampling.** The backend requires at least two source
depth samples; the static tutorial builds 10 and 11 km while demonstrating
a query at 10 km. A single-depth request does not make that Fortran
constraint disappear.

**MPI postprocessing and task groups.** SPGRN MPI calls update shared
metadata from each rank without a dedicated postprocessing barrier.
QSSP MPI routines can access rank entries before excluding unused ranks
in uneven groups. Multi-node MPI requires a separate cluster validation,
including divisible task groups and metadata/postprocessing coordination.
The local tutorial checks do not establish multi-node correctness.

**Direct QSSP convenience calls.** The separate `read_by_qssp` pathway
includes source coordinates in its cache key but writes a source fixed
at latitude/longitude 0/0; a nonzero-coordinate request must not be
interpreted as a validated relocation of that source. Its completion
hash also omits the focal mechanism, so changing mechanism can reuse
old output when completion reuse is enabled. Use the precomputed-library
tutorial as the documented complete workflow.

**Boundary queries and reuse.** Reader endpoint clamping/validation is not
uniform, and completion flags do not hash models or parameters. Query
inside the computed grid and use a new directory after changing physical
or numerical settings.

## Report a reproducible issue

Include the package/commit version, operating system, Python version,
installation method, relevant compiler/JDK/MPI versions, backend name,
smallest driver, model, `green_lib_info.json` and executable error output.
State the expected physical quantity, component order, source amplitude
and time reference. Attach the smallest failing data subset or describe
how to generate it; a screenshot alone cannot establish a numerical
problem.
