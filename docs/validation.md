# Tutorial validation record

All six backend workflows were executed locally on 14 September 2026 using the
scripts included in this repository. This record describes the actual Windows
runs used to produce the tutorial figures. A separate GitHub Actions run also
built and executed all six workflows on Linux, Windows and macOS, as recorded
below. The current spherical examples were recalculated after the harmonic
cutoff audit, and additional QSEIS regional workflows were run at 300/600/900 km.

## Fresh 2 Hz displacement and stress comparison

On 15 September 2026, all five dynamic backends and two Gaussian smoothing
controls completed new calculations at source depth 10 km, receiver depth
1 km and distances 300/600/900 km. The shared AK135-FC model, 0.25 s sampling,
2 Hz maximum and 1.25 s sin² source differ from the lightweight tutorials
below. See the [fresh report](guides/backend-comparison.md) and
[public validation summary](_static/comparisons/2026-09-15/result-summary.json).

| Calculation | Complete processed displacement shape | Runtime | Recorded output size, decimal |
| --- | --- | ---: | ---: |
| QSEIS06, native 0.05 smoothing | `(3, 3, 2048)` | 971 s | 2.3 MB |
| QSEIS2025, 0.05 smoothing control | `(3, 3, 2048)` | 983 s | 5.0 MB |
| QSEIS06, isolated point-source control, eight frequency tasks | `(3, 3, 2048)` | 609 s | 24.9 MB |
| QSEIS2025, point source, eight frequency tasks | `(3, 3, 2048)` | 689 s | 52.8 MB |
| SPGRN2012 | `(3, 3, 2048)` | 1419 s | 1.207 GB |
| SPGRN2020 | `(3, 3, 2048)` | 1542 s | 2.777 GB |
| QSSP2020, new spectra and six basis-source syntheses | `(3, 3, 2048)` | 1908 s | 11.715 GB |

QSEIS2025 and QSSP2020 stress arrays have shape `(3, 6, 2048)` and units Pa.
The public comparison arrays contain 1601 samples on 0–400 s. Complete
native records, finite values, source area/centroid/spectrum, executable
hashes and QSSP's native rate/direct-output integration identity were checked.
All 636 published metric rows were independently recomputed from the public
arrays, with a maximum discrepancy of `1.42e-14` percentage points.

Matched-radius QSEIS versions were bit-identical in displacement at both
ratios 0.05 and 0. The latter QSEIS06 result requires an isolated control
build; it is not the released default. The isolated QSEIS builds increase
layer capacity to 2048, and the point calculations sum eight disjoint
frequency ranges with validated linear superposition. Their source diffs
and executable hashes are included in the public data record. No installed
solver or package default was changed by these controls.

These are **Windows 11** measurements using Python 3.12.13, NumPy 2.3.5,
SciPy 1.18.0, Matplotlib 3.11.0 and ObsPy 1.5.0. Some calculations overlapped,
so their elapsed times are not a controlled performance benchmark. Sizes
include each result directory's retained data and exclude later website
assets. These scientific runs were not added to per-PR tutorial CI and do
not extend the older cross-platform validation claim to this new setup.

The {download}`portable plotting script <../examples/plot_documented_comparison.py>`
replots the committed comparison arrays with NumPy and Matplotlib. It does
not recompute a native Green's library.

## Environment and results for the lightweight tutorials

The local environment was Windows 11 build 26200, 64-bit CPython 3.12.13, NumPy
2.3.5, SciPy 1.18.0, pandas 3.0.0, Matplotlib 3.11.0 and ObsPy 1.5.0. Native
executables were installed in the same Conda environment. Commands used
`conda run -n ffipy python examples/<script>.py`.

Times below measure preparation, solver execution, output reading, checks and
plotting after imports; they exclude interpreter startup and package installation.
Sizes include the retained solver ASCII files, binary files, spectra and figures.
They are measurements for these small examples, not performance guarantees.

| Workflow | Command argument | Validated output shape | Time | Output size | Machine-readable record |
|---|---|---|---:|---:|---|
| QSEIS2025 introduction | `--observables all` | displacement `(3, 3, 201)`; strain/stress `(3, 6, 201)` | 13.3 s | 1.21 MiB | [JSON](_static/examples/qseis2025.json) |
| QSEIS06 introduction | default | `(3, 3, 256)` | 13.0 s | 0.38 MiB | [JSON](_static/examples/qseis06.json) |
| SPGRN2012 | default, impulse + shared STF | `(3, 3, 256)` | 44.7 s | 203.41 MiB | [JSON](_static/examples/spgrn2012.json) |
| SPGRN2020 | default, complete wavefield | `(3, 3, 256)` | 62.5 s | 266.14 MiB | [JSON](_static/examples/spgrn2020.json) |
| QSSP2020 | default, harmonics 2000/8000 | `(3, 3, 256)` | 86.3 s | 657.98 MiB | [JSON](_static/examples/qssp2020.json) |
| EDGRN2 → EDCMP2 | default, including both solvers | `(3, 3)` | 3.1 s | 0.17 MiB | [JSON](_static/examples/edgrn_edcmp.json) |
| QSEIS2025 regional | `--regional --observables all` | displacement `(3, 3, 256)`; strain/stress `(3, 6, 256)` | 187.4 s | 3.66 MiB | [JSON](_static/examples/qseis2025-regional.json) |
| QSEIS06 regional | `--regional` | `(3, 3, 256)` | 186.3 s | 1.18 MiB | [JSON](_static/examples/qseis06-regional.json) |

The dynamic array axes are distance, component and sample. The static axes are
distance and component. The default displacement-only QSEIS2025 command was also
executed independently in a fresh output directory.

The QSEIS2025 run and figures were refreshed after cropping every exported
observable to 0–100 s inclusive. At 0.5 s spacing, the NPZ arrays contain
201 samples; each was checked against the first 201 samples of the previous
full waveform. Its native library still contains 256 samples over 127.5 s.
The JSON records the exported interval in `output_time_range_s` separately
from the native `time_window_s`. The refreshed displacement-only run took
11.4 s.

## Matched frequency band, mechanism and temporal source

All five regional scripts were run in fresh directories with 4 s sampling,
1024 native FFT samples, frequency spacing 1/4096 Hz and a 0.125 Hz maximum.
Native spectrum headers (SPGRN/QSSP) and input/output grids (QSEIS) were checked:
512 nonnegative bins are computed through 511/4096 Hz, and Nyquist is zero.
The three old spherical libraries were rejected because their native headers
retained only 257 bins, even before checking the declared JSON settings.

Every script uses the same strike/dip/rake (30/45/90 degrees), azimuth
(30 degrees), scalar moment (10^15 N m), source depth (10 km), receiver depth
(0 km), and effective normalized 64 s sin-squared moment-rate pulse.
SPGRN2012 now uses full-wavefield spectra and the native zero-duration impulse
branch. Its complete 1024-point velocity is convolved forward with the analytic
source at the native complex frequencies, then integrated once and cropped to
256 samples. Quadrature independently verified the full-band source transform
with relative L2 error 1.75e-15. Input, native-velocity and four NPZ hashes were
unchanged by compatible reuse. The old positive-duration library was rejected.

QSEIS retains its compensated custom source. Its independent spectral check
now covers 0–0.125 Hz. Both QSEIS native calculations completed normally;
their displacement samples were identical and their exported tensor arrays
were finite. QSEIS2025's original short example still exports 0–100 s.

The [earlier 64 s comparison](guides/backend-comparison-64s.md) uses the common native
origin-time intervals [3,500], [40,500] and [78,500] s, interpolated to 1 s without
amplitude or time-shift fitting. Relative ENU L2 differences from SPGRN2020 are:

| Backend | 300 km | 600 km | 900 km |
|---|---:|---:|---:|
| SPGRN2012 | 1.2266% | 0.0739% | 0.0534% |
| QSSP2020 | 1.2772% | 0.4085% | 2.3022% |
| QSEIS06 / QSEIS2025, default spatial smoothing | 12.7410% | 19.5740% | 22.7979% |
| QSEIS2025, point-source control | 11.4926% | 17.0758% | 19.7176% |

The point-source control only changes QSEIS2025's Gaussian spatial-source ratio
from 0.05 to zero; its model and temporal source samples are unchanged. That
additional Windows run took 647.470 s and retained 1,005,149 bytes in 27 files;
its `(3,3,256)` displacement was finite. See the
[point-source record](_static/examples/qseis2025-point-source.json).
It reduces only part of the residual and was not added to every CI calculation.
QSEIS06 has the ratio fixed at 0.05 in the solver. No Fortran algorithm was changed.

The fresh runtimes and output sizes above are preserved separately from the
short compatible rereads used to check the final scripts and metadata. Archived
JSON records identify those post-run checks explicitly. Older frequency-cutoff,
harmonic and STF comparisons remain in the
[historical audit](guides/backend-comparison-history.md).
The initial cross-platform archive below predates these updated settings.
## Initial cross-platform installation and execution

[GitHub Actions run 34801126315](https://github.com/Zhou-Jiangcheng/pygrnwang/actions/runs/34801126315)
built commit `0d38d44` before the QSEIS2025 100 s crop. These archived records
therefore retain its original 256-sample exported waveforms; the current
201-sample exports are verified by the refreshed local run above. The CI run
installed source in fresh Conda environments, using Python
3.12.14 on Linux x86-64, Windows x86-64 and macOS arm64. All six workflows
passed on every platform. The Windows job used
`gfortran=15.2.0=hf1b5d6d_19`; the installation guide records why that exact
compiler build is selected.

Each cell below gives calculation time / retained output size. The timing
boundary is the same as the local table; compiler installation and compilation
are excluded. QSEIS2025 includes displacement, strain and stress.

| Workflow | Linux | Windows | macOS |
|---|---:|---:|---:|
| QSEIS2025 | 15.4 s / 1.19 MiB | 23.2 s / 1.21 MiB | 9.1 s / 1.19 MiB |
| QSEIS06 | 14.2 s / 0.37 MiB | 19.9 s / 0.37 MiB | 8.3 s / 0.37 MiB |
| SPGRN2012 | 17.3 s / 12.54 MiB | 18.5 s / 12.54 MiB | 9.8 s / 12.54 MiB |
| SPGRN2020 | 17.4 s / 17.09 MiB | 19.2 s / 17.10 MiB | 9.7 s / 17.09 MiB |
| QSSP2020 | 25.4 s / 63.72 MiB | 27.9 s / 63.73 MiB | 19.1 s / 63.72 MiB |
| EDGRN2 → EDCMP2 | 0.3 s / 0.16 MiB | 1.6 s / 0.17 MiB | 0.3 s / 0.16 MiB |

The [preserved CI records](_static/examples/ci-validation.json) include all
18 run summaries, exact platform/Python/dependency versions, dimensions,
component names, units, amplitudes and output sizes. The 24 uploaded NPZ
arrays were reopened and checked for finite values, matching shapes,
component labels and units. The workflow also retains plots and selected
native inputs as downloadable artifacts for 14 days; the JSON records here
remain part of the documentation after those artifacts expire.

These checks exercise editable source installation. They do not claim a
standard-wheel validation, multi-node MPI validation or scientific convergence
for arbitrary models. They also do not exercise every Python version in the
package support range.

## Documentation checks

The isolated Python 3.12 documentation environment uses the locked requirements
in this repository. Sphinx's strict HTML build passed without warnings; the
explicit API checker covered 92 functions, classes and methods. Internal files
and anchors, MathJax formulas, figure loading, search, and desktop (1440 px)
and mobile (390 px) layouts were checked. All 25 modified Python source files
had identical ASTs after removing docstrings.

The focused TauP suite passed 12 tests, with one installed-wheel-only test
skipped in the source environment. A source-distribution archive was checked
to include documentation, figures and executable tutorial sources, excluding
generated HTML and calculation libraries.

## What was checked

- Every dynamic example returned three-component displacement at three distances,
  with finite values and a nonzero waveform. The QSEIS2025 introduction exports
  201 samples over 0–100 s; its regional variant and other dynamic examples
  export 256 samples. QSEIS2025 strain and stress have six components on the
  same exported distance/time grid as displacement in each mode.
- EDGRN generated the layered kernels; EDCMP used those kernels for all five
  mechanism bases at both source depths. ASCII-to-binary conversion produced
  the bulk shape `(2, 1, 5, 5, 3)`, and the three queried displacements were finite.
- Component labels and time axes were checked against the readers. Vector figures
  use east, north, up; QSEIS2025 tensors use EE, EN, EU, NN, NU, UU. SPGRN2012
  uses its distance-dependent native start time, and the revised SPGRN2020
  example uses the origin-time starts stored in its native binary headers.
- The static material lookup uses four-column `noQ.nd`, avoiding an invalid
  six-column material reshape. Its result is explicitly multiplied by seismic
  moment after EDCMP's unit-moment normalization.
- Figures were inspected for readable axes, legends, component labels and layout.
  The spherical examples use a longer spectral window than output window to
  reduce periodic end-of-window contamination.

## Reproduce and interpret the record

Run the commands in the [backend tutorials](backends/index.md). Each script writes
its own `summary.json`, with the environment, elapsed time, array dimensions,
maximum absolute value and output size. `--reuse` repeats the reader/plot checks
without native recomputation and writes `summary-reuse.json`.

The base runs check installation and workflow; the additional comparison records verified frequency grids, temporal sources and specific numerical checks. The tutorial model uses the bundled
AK135 elastic structure with explicit constant Qp = 600 and Qs = 300. Regional
examples truncate it at 809.5 km; spherical examples retain the full structure.
The figures are not a claim that these different discretizations, source time
functions and solver approximations produce interchangeable research results.
Check convergence in spectral window, bandwidth, spatial sampling, harmonic or
wavenumber settings and model resolution for the scientific problem of interest.

EDGRN requires two source depths. Static queries are kept within the distance
table because EDCMP finite source corners can cross a table edge. Multi-node MPI
and large production calculations were not exercised by these local runs.
