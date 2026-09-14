# Parallel execution and resuming

Run the sequential tutorial first, then enlarge the calculation. Parallel
workers distribute independent backend jobs; they do not change the physical
resolution or compute a single backend solve across multiple Python cores.

## Single-node multiprocessing

The preprocessors store `processes_num` and write task groups. Corresponding
`create_grnlib_*_parallel` functions use a process pool to execute those
jobs. Select a worker count that fits both CPU and memory: each Fortran
process has its own working arrays and can be memory-intensive.

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

The current SPGRN MPI routines update common metadata from each rank;
their end-of-run metadata handling should be reviewed for the target
cluster before a production multi-node run. QSSP MPI routines also have
limitations with uneven final groups; see
[known implementation limitations](troubleshooting.md#known-implementation-limitations).
Multi-node MPI is not validated by the local tutorial runs or ordinary
documentation builds.

## Resuming a calculation

`check_finished=True` reuses backend directories that have a `.finished`
marker and expected supporting files. This is an existence-based restart
mechanism, not a parameter checksum. A marker records a prior executable
run; inspect its output and the numerical files before treating it as
scientific validation.

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
