# Executable tutorials

**QSEIS06 and SPGRN2012 are deprecated in pygrnwang.** Their examples remain
available for existing workflows. For new calculations use `qseis2025.py`
and `spgrn2020.py`, respectively. Rebuild libraries and validate the new
backend's settings and time origin when migrating.

Install `pygrnwang` into the Python environment first, including its native
executables. Run these commands from the repository root:

```sh
python examples/qseis2025.py
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025-all
python examples/qseis06.py
python examples/qseis2025.py --regional --observables all
python examples/qseis06.py --regional
python examples/spgrn2012.py
python examples/spgrn2020.py
python examples/qssp2020.py
python examples/edgrn_edcmp.py
python examples/compare_backends.py
```

On Windows with Conda, activate the environment first or prefix each command
with `conda run -n YOUR_ENV`. Do not invoke an environment's Python executable
without activation. Each script has a `__main__` guard and uses the serial
builder. No MPI launcher is required.

`--output-dir PATH` changes the output root; relative paths are resolved against
the directory from which the script is started. The default is
`examples/output/<backend>` next to the scripts. Output paths are made absolute
before calling library functions that change the working directory.

Each run writes `ak135_tutorial.nd`, `library/`, `disp.npz`, `disp.png` and
`summary.json`. QSEIS2025 `--observables all` also writes strain and stress arrays
and figures. The arrays contain physical values for M0 = 10^15 N m. They record
the components, units, distance grid and (for dynamic outputs) each trace's time
axis. The default QSEIS2025 introduction crops displacement, strain and stress to 0–100 s inclusive
before saving NPZ arrays and figures; every plot uses a fixed 0–100 s x-axis.
At 0.5 s spacing, its saved displacement has shape `(3, 3, 201)` and its
saved tensor outputs have shape `(3, 6, 201)`. Figures use Matplotlib's
noninteractive Agg backend.

`--regional` selects 300, 600 and 900 km for QSEIS06/QSEIS2025, with 4 s
sampling, a 4092 s native window, a 64 s wavelet and the flat-Earth transformation
enabled. It saves 0–1020 s inclusive: displacement `(3, 3, 256)` and QSEIS2025
tensors `(3, 6, 256)`. These distances need a longer window because their main
arrivals extend beyond 100 s. The default output directories are
`examples/output/qseis06-regional` and `examples/output/qseis2025-regional`.

Regional QSEIS uses a custom, normalized 64 s sin-squared moment-rate function
(`wavelet_type=0`). `source_time_function.py` writes 1024 input nodes with the
numerical damping precompensated, so the effective physical pulse matches
SPGRN2020/QSSP rather than merely sharing their duration. The input samples
must not be normalized again: their area is about 0.964709, while the effective
physical pulse has unit area and a 32 s centroid. The scripts read rate
kernels and integrate them once before exporting displacement, strain or stress.
The helper archives the source samples and numerical checks in
`source_time_function.npz/json`; reuse verifies these and the generated input.

After the five dynamic workflows above, `compare_backends.py` reads their NPZ
files, aligns source-origin time and plots 0–500 s in
`examples/output/backend-comparison/`. It checks the paired QSEIS results and
records cross-backend differences without fitting time shifts or amplitudes.
Its path options accept independently calculated libraries. See the
[comparison guide](../docs/guides/backend-comparison.md) for physical limits.

Use `--reuse` only with the same script and observables as a successful earlier
run. It reads and plots the existing library without rebuilding it and saves
its verification report as `summary-reuse.json`. Run a different parameter set
in a fresh output directory; a fresh run rejects an existing library because binary components could belong
to earlier solver settings.

## Model and numerical choices

`common.py` writes the bundled `pygrnwang.ak135fc.s` elastic structure with
illustrative constant Qp = 600 and Qs = 300 appended to the numeric rows.
These constants are not the published AK135-F attenuation profile. Fluid-layer
Qs is unused where Vs is zero. Discontinuity labels are retained for TauP.

QSEIS and EDGRN use the first 24 numeric rows (down to 809.5 km) as a regional
layered half-space approximation; the spherical solvers use the full model.
The default QSEIS introductions use 0.5 s sampling and a native 127.5 s window
(256 samples). QSEIS2025 then crops its introductory outputs to 0–100 s
(201 samples); QSEIS06 retains all 256 samples. The regional options instead
use the longer native and exported windows described above.
The spherical tutorials use 4 s sampling, a 4092 s spectral window and a
1020 s output window, with a
64 s source duration and 0.0625 Hz cutoff. The longer spectral window reduces
periodic contamination in the shorter output window. SPGRN2020 uses
`max_slowness=0` to select its complete-wavefield branch. QSSP uses
`min_harmonic=2000, max_harmonic=8000`; raising the former to 4000 changed
these example waveforms by less than 0.008%. Both settings affect the
low-frequency content and spatial summation, so they need renewed convergence
checks for other source depths, distances or bands. SPGRN2020 plots the actual
native origin-time starts, including their integer-second rounding.

EDGRN requires at least two source depths, so the static library uses 10 and
11 km while the plotted query uses 10 km. Its distance grid spans 0–120 km,
with queries at 30, 60 and 90 km kept away from the table boundaries. EDCMP's
finite source corners can exceed a boundary even when the reference point
lies on it. The material lookup uses the generated four-column `library/noQ.nd`.

Each dynamic output is checked for three distances, the requested component
count, finite values and a nonzero signal at each distance. QSEIS2025's
default saved outputs contain exactly 201 samples over 0–100 s inclusive;
its regional outputs and the other dynamic tutorials retain 256 samples.
Static output is checked for shape `(3, 3)`, finite values and a nonzero result.
A successful tutorial checks installation and data flow; research calculations
also require physical validation and parameter convergence tests.
