# Tutorial validation record

All six backend workflows were executed locally on 14 September 2026 using the
scripts included in this repository. This record describes the actual Windows
runs used to produce the tutorial figures. A separate GitHub Actions run also
built and executed all six workflows on Linux, Windows and macOS, as recorded
below.

## Environment and results

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
| QSEIS2025 | `--observables all` | displacement `(3, 3, 256)`; strain/stress `(3, 6, 256)` | 14.6 s | 1.23 MiB | [JSON](_static/examples/qseis2025.json) |
| QSEIS06 | default | `(3, 3, 256)` | 13.0 s | 0.38 MiB | [JSON](_static/examples/qseis06.json) |
| SPGRN2012 | default | `(3, 3, 256)` | 9.4 s | 12.54 MiB | [JSON](_static/examples/spgrn2012.json) |
| SPGRN2020 | default | `(3, 3, 256)` | 9.7 s | 17.10 MiB | [JSON](_static/examples/spgrn2020.json) |
| QSSP2020 | default, including spectra | `(3, 3, 256)` | 14.6 s | 63.74 MiB | [JSON](_static/examples/qssp2020.json) |
| EDGRN2 → EDCMP2 | default, including both solvers | `(3, 3)` | 3.1 s | 0.17 MiB | [JSON](_static/examples/edgrn_edcmp.json) |

The dynamic array axes are distance, component and sample. The static axes are
distance and component. The default displacement-only QSEIS2025 command was also
executed independently in a fresh output directory.

## Cross-platform installation and execution

[GitHub Actions run 34801126315](https://github.com/Zhou-Jiangcheng/pygrnwang/actions/runs/34801126315)
built commit `0d38d44` from source in fresh Conda environments, using Python
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

- Every dynamic reader returned three-component displacement at three distances,
  with 256 samples, finite values and a nonzero waveform. The QSEIS2025 extension
  returned six-component strain and stress with the same distance/time grid.
- EDGRN generated the layered kernels; EDCMP used those kernels for all five
  mechanism bases at both source depths. ASCII-to-binary conversion produced
  the bulk shape `(2, 1, 5, 5, 3)`, and the three queried displacements were finite.
- Component labels and time axes were checked against the readers. Vector figures
  use east, north, up; QSEIS2025 tensors use EE, EN, EU, NN, NU, UU. SPGRN2012
  uses its distance-dependent native start time, and SPGRN2020 plots time relative
  to its library P arrival.
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

These are installation and workflow checks. The tutorial model uses the bundled
AK135 elastic structure with explicit constant Qp = 600 and Qs = 300. Regional
examples truncate it at 809.5 km; spherical examples retain the full structure.
The figures are not a claim that these different discretizations, source time
functions and solver approximations produce interchangeable research results.
Check convergence in spectral window, bandwidth, spatial sampling, harmonic or
wavenumber settings and model resolution for the scientific problem of interest.

EDGRN requires two source depths. Static queries are kept within the distance
table because EDCMP finite source corners can cross a table edge. Multi-node MPI
and large production calculations were not exercised by these local runs.
