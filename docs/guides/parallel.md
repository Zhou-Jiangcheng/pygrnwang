# Parallel execution and resuming

Run the sequential tutorial first, then enlarge the calculation. Parallel
workers distribute independent backend jobs; they do not change the physical
resolution or compute a single backend solve across multiple Python cores.

## Single-node multiprocessing

The preprocessors store `processes_num` and write task groups. Corresponding
`create_grnlib_*_parallel` functions use a process pool to execute those
jobs. Select a worker count that fits both CPU and memory: each Fortran
process has its own working arrays and can be memory-intensive.

### Memory

QSEIS06, QSEIS2025 and EDCMP2 use static Fortran arrays, so every process
needs a fixed amount of memory whatever the grid size: about 0.7, 1.3 and
0.95 GiB. Before starting, their `create_grnlib_*_parallel` functions issue a
`RuntimeWarning` when min(`processes_num`, number of jobs) processes of
`memory_per_job_gb` (default 0.75, 1.4 and 1.0 GiB) exceed the currently
available memory; on Linux this includes a Slurm or container cgroup limit.
The warning gives the largest worker count that fits. The run goes on: jobs
that run out of memory fail and are computed again after the others (see
below). To avoid those retries, reduce `processes_num` in the preprocessor.
`memory_per_job_gb=None` skips the warning. EDGRN2 needs about 10 MiB per
process. QSSP2020, SPGRN2020 and SPGRN2012 allocate their arrays from the
input, so their warning is off by default; pass the peak memory of one job
as `memory_per_job_gb` to enable it.

### Failed jobs, retries and Ctrl+C

A failed job does not stop the others, for example a job the operating
system killed or refused to start for lack of memory. After all jobs ran,
the builders compute the jobs that did not complete again, with the same
worker count, up to `max_retries=2` times. They then check the library and
raise `RuntimeError` listing the jobs that still fail, with their logs, and
the incomplete files. Ctrl+C stops the run at once: no new job starts and
the running backend processes are killed.

Check an existing library with the matching `check_grnlib_*` function:

```python
from pygrnwang.create_qseis2025_bulk import check_grnlib_qseis2025

problems = check_grnlib_qseis2025(path_green)  # [] when complete
```

Each check requires every file the readers use, for every prepared job, in
binary form with the expected size or as complete ASCII output, and the
travel-time tables where the backend has them. `check_values=True` also
reads the binary files and reports NaN or infinite values.

Every backend reads and writes its paths through 160-character variables.
The preprocessors raise `ValueError` when a job path would be longer, and
QSEIS preprocessors also when `N_each_group` exceeds 101 distances.

Keep the calculation driver under `if __name__ == "__main__":`.
This is required for Windows spawning and makes scripts portable. Use
absolute model and output paths because low-level calls can change working
directories. Do not share one job directory between simultaneous independent
runs.

A single source/receiver pair offers little parallelism for SPGRN. QSEIS
also distributes distance groups; QSSP distributes source/receiver pairs
for spectra and six tensor-source jobs for time-domain outputs. More
workers than available tasks provide no speedup.

## Multi-node MPI

MPI is optional. Install `mpi4py` and an MPI runtime compatible with the
cluster, and launch a dedicated driver with the site's launcher, for example
`mpiexec -n 4 python run_mpi.py`. A scheduler may require `srun` or other
site-specific allocation flags.

Use the `create_grnlib_*_parallel_multi_nodes` API, not the multiprocessing
API inside every MPI rank. Preprocess **once** before launching ranks, on
a filesystem visible with the same absolute paths on all nodes. Every node
needs the package, backend binaries and their runtime libraries.

The MPI routines associate ranks with the precomputed task groups. Match
the number of ranks to the size of the first group, which may be smaller
than requested `processes_num` when there are few jobs. Do not independently
rerun preprocessing on every rank.

QSSP exposes separate spectral and time-domain MPI functions. Complete the
spectral launch successfully before launching the time-domain driver.
Run format conversion once after all ranks finish. SPGRN2012 additionally
needs Python travel-time tables after its MPI run, and EDCMP's MPI path
needs explicit conversion if binary output is desired.

Every rank runs one backend process at a time. Before any job starts, one
rank per node warns when the ranks of its node may need more than the
available memory (`memory_per_job_gb`); the run goes on. A failed job does not stop the
others; after all ranks finish, the jobs that did not complete are shared
among the ranks and computed again, up to `max_retries` times. Rank 0 then
records SPGRN metadata, writes SPGRN2012 travel-time tables, checks the
library and raises `RuntimeError` if it is incomplete. QSSP MPI functions
accept more ranks than the group width; the extra ranks help with retries.
Multi-node MPI is not validated by the local tutorial runs or ordinary
documentation builds.

## Resuming a calculation

`check_finished=True` reuses backend directories that have a `.finished`
marker and complete output files; every other job is computed again. This
is an existence-based restart mechanism, not a parameter checksum. A marker
records a prior successful executable run: it is written only when the
executable exits with code 0 and without a Fortran `STOP` error message;
failed runs leave their log in `.failed` instead. A job that runs again
first deletes its earlier output, so a failure cannot leave old files that
pass the library check. Inspect the output and the numerical files before
treating a marker as scientific validation.

After a builder raised `RuntimeError`, fix the cause (for example reduce
`processes_num` if jobs ran out of memory) and rerun with
`check_finished=True` to compute only the unfinished jobs; an MPI build can
also be completed with the local parallel builder on one node.

Reuse only with identical model, grids, time window, sampling, source
function, enabled observables and solver settings. After changing any of
these, use a new output directory or deliberately rebuild the affected
stages. Do not preserve stale metadata with newly generated Green's files.

QSSP `cal_spec=False` is valid only when compatible spectral files already
exist. The first calculation must use `cal_spec=True`. SPGRN `cal_gf=1`
computes/updates spectra, while `0` requests compatible existing spectra.
Changes to spectral model/sampling settings require rebuilding spectra.

Travel-time tables have a separate `check_finished_tpts_table` setting
where exposed, and `create_tpts_table(check_finished=True)` only checks
whether its two files exist. New model or grid means new tables.

The tutorial scripts' `--reuse` mode reuses a completed tutorial library
and regenerates its saved results through the documented reader path,
writing `summary-reuse.json` while preserving `summary.json`.
Use the backend completion flags in your own driver for interrupted large
calculations; keep the original generated inputs and inspect incomplete jobs.

## Record a useful run

Save the host/platform, Python/package/dependency versions, worker count,
backend elapsed time, resulting bytes, input model and exact generated
files. The tutorials write `summary.json` with reproducibility and array
checks. For a production grid, also retain backend logs and review completed
file counts before deleting intermediate ASCII results.
