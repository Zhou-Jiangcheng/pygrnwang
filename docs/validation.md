# Tutorial validation record

All six backend workflows were executed locally on 14 September 2026 using the
scripts included in this repository. This record describes the actual Windows
runs used to produce the tutorial figures. A separate GitHub Actions run also
built and executed all six workflows on Linux, Windows and macOS, as recorded
below. The current spherical examples were recalculated after the harmonic
cutoff audit, and additional QSEIS regional workflows were run at 300/600/900 km.

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
| QSEIS2025 introduction | `--observables all` | displacement `(3, 3, 201)`; strain/stress `(3, 6, 201)` | 13.3 s | 1.21 MiB | [JSON](_static/examples/qseis2025.json) |
| QSEIS06 introduction | default | `(3, 3, 256)` | 13.0 s | 0.38 MiB | [JSON](_static/examples/qseis06.json) |
| SPGRN2012 | default | `(3, 3, 256)` | 9.4 s | 12.54 MiB | [JSON](_static/examples/spgrn2012.json) |
| SPGRN2020 | default, complete wavefield | `(3, 3, 256)` | 17.1 s | 116.02 MiB | [JSON](_static/examples/spgrn2020.json) |
| QSSP2020 | default, harmonics 2000/8000 | `(3, 3, 256)` | 38.2 s | 330.63 MiB | [JSON](_static/examples/qssp2020.json) |
| EDGRN2 → EDCMP2 | default, including both solvers | `(3, 3)` | 3.1 s | 0.17 MiB | [JSON](_static/examples/edgrn_edcmp.json) |
| QSEIS2025 regional | `--regional --observables all` | displacement `(3, 3, 256)`; strain/stress `(3, 6, 256)` | 248.4 s | 3.66 MiB | [JSON](_static/examples/qseis2025-regional.json) |
| QSEIS06 regional | `--regional` | `(3, 3, 256)` | 242.4 s | 1.18 MiB | [JSON](_static/examples/qseis06-regional.json) |

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

## Harmonic correction and regional comparison

The revised SPGRN2020 example uses `max_slowness=0`; QSSP2020 uses
`min_harmonic=2000, max_harmonic=8000`. Both current scripts were executed in
fresh directories to update their figures and JSON records. The
[controlled harmonic audit](guides/backend-comparison.md) includes nine
additional parameter-variation calculations: increasing QSSP's minimum from
2000 to 4000 at fixed maximum 8000 changed the complete saved waveform by
less than 0.008%. The original parameter sets had produced finite files but
had not established waveform convergence.

The QSEIS regional runs use 300/600/900 km, 4 s sampling, the 24-row model
and the flat-Earth transformation. Native arrays contain 1024 samples over
4092 s, with 256 samples over 0–1020 s exported. The first regional runs used
the native real-frequency wavelet and differed from SPGRN2020 by 10.82%,
17.52% and 20.74%; those measurements are retained in the
[before-STF record](_static/examples/backend-comparison-before-stf.json).

The current regional examples instead supply 1024 custom moment-rate nodes
with damping precompensated to match SPGRN2020's physical 64 s source.
The effective pulse has unit integral, a 32 s centroid, and time-domain and
spectral relative L2 errors below 1e-5 against the analytic target. These are
checked independently of the resulting seismogram amplitude; the example
never fits a waveform scale factor to improve agreement. It reads velocity,
strain rate and stress rate explicitly and integrates each once. See the
[STF record](_static/examples/source-time-function.json) and the updated
[backend comparison](guides/backend-comparison.md) for the final waveforms.

Both custom-source QSEIS runs completed normally. Displacement arrays again
matched sample for sample; the additional strain/stress outputs were finite.
Relative L2 differences from SPGRN2020 were 12.68%, 19.75% and 22.95%, so
matching the source did not explain or remove the earlier regional discrepancy.
The source-area excess in the old run had partly offset other differences.
These results are preserved without fitting a waveform scale or time shift.

The default near-distance QSEIS results were reread and compared with their
previous arrays: displacement and the QSEIS2025 tensor outputs were unchanged.
Trying to reuse a near-distance library in regional mode was rejected, as was
requesting missing tensors from a displacement-only library. A new tutorial
run also rejects an existing library directory to prevent old completion
markers from being mistaken for freshly calculated data.

`examples/compare_backends.py` saves three overlay figures and a
[comparison record](_static/examples/backend-comparison.json), using each
distance's common physical interval through 500 s. It checks the paired QSEIS
outputs and requires the revised QSSP/SPGRN2020 relative L2 difference to
remain below 5% for this example. Model, source-spectrum and integration
limitations are explained in the [comparison guide](guides/backend-comparison.md).
The example workflow now runs both regional QSEIS commands and this comparison
on each supported CI platform; the archived run below predates these additions.

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

The base runs check installation and workflow; the additional comparison records specific numerical checks. The tutorial model uses the bundled
AK135 elastic structure with explicit constant Qp = 600 and Qs = 300. Regional
examples truncate it at 809.5 km; spherical examples retain the full structure.
The figures are not a claim that these different discretizations, source time
functions and solver approximations produce interchangeable research results.
Check convergence in spectral window, bandwidth, spatial sampling, harmonic or
wavenumber settings and model resolution for the scientific problem of interest.

EDGRN requires two source depths. Static queries are kept within the distance
table because EDCMP finite source corners can cross a table edge. Multi-node MPI
and large production calculations were not exercised by these local runs.
