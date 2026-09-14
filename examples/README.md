# Executable tutorials

Install `pygrnwang` into the Python environment first, including its native
executables. Run these commands from the repository root:

```sh
python examples/qseis2025.py
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025-all
python examples/qseis06.py
python examples/spgrn2012.py
python examples/spgrn2020.py
python examples/qssp2020.py
python examples/edgrn_edcmp.py
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
axis. Figures use Matplotlib's noninteractive Agg backend.

Use `--reuse` only with the same script and observables as a successful earlier
run. It reads and plots the existing library without rebuilding it and saves
its verification report as `summary-reuse.json`. Run a different parameter set
in a fresh output directory; existing binary components could otherwise belong
to earlier solver settings.

## Model and numerical choices

`common.py` writes the bundled `pygrnwang.ak135fc.s` elastic structure with
illustrative constant Qp = 600 and Qs = 300 appended to the numeric rows.
These constants are not the published AK135-F attenuation profile. Fluid-layer
Qs is unused where Vs is zero. Discontinuity labels are retained for TauP.

QSEIS and EDGRN use the first 24 numeric rows (down to 809.5 km) as a regional
layered half-space approximation; the spherical solvers use the full model.
QSEIS uses 0.5 s sampling and a 127.5 s output window. The spherical tutorials
use 4 s sampling, a 4092 s spectral window and a 1020 s output window, with a
64 s source duration and 0.0625 Hz cutoff. The longer spectral window reduces
periodic contamination in the shorter output window. QSSP's harmonic limit
is 800; these settings demonstrate the workflow and are not a convergence
study for arbitrary distances or near-field static offsets.

EDGRN requires at least two source depths, so the static library uses 10 and
11 km while the plotted query uses 10 km. Its distance grid spans 0–120 km,
with queries at 30, 60 and 90 km kept away from the table boundaries. EDCMP's
finite source corners can exceed a boundary even when the reference point
lies on it. The material lookup uses the generated four-column `library/noQ.nd`.

Each dynamic output is checked for three distances, the requested component
count, exactly 256 samples, finite values and a nonzero signal at each distance.
Static output is checked for shape `(3, 3)`, finite values and a nonzero result.
A successful tutorial checks installation and data flow; research calculations
also require physical validation and parameter convergence tests.
