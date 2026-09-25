# Quickstart: QSEIS2025

This example computes a small Green's-function library, synthesizes displacement and writes a figure. It uses the AK135 model content included with pygrnwang, one source depth, one receiver depth and three distances. Computation is serial. The source is 10 km deep, receivers are at the surface and distances are 300, 600 and 900 km. The native library uses a 4 s sampling interval, a 0.125 Hz Nyquist limit and a 4092 s window containing 1024 samples, with the flat-Earth transformation enabled. After synthesis, the example saves 0–1020 s inclusive, giving 256 samples without changing the underlying library. The layered calculation uses the first 24 numeric model rows, down to 809.5 km. Constant Qp=600 and Qs=300 are illustrative tutorial choices, not the full AK135-F attenuation model. The source is a normalized 64 s squared half-sinusoid moment-rate pulse. The mechanism is strike 30°, dip 45°, rake 90°, scaled to M0 = 10^15 N m.

## 1. Prepare the environment

Follow [installation](installation.md), including cloning the repository to obtain `examples/`. Run the following commands from the repository root in the activated environment.

## 2. Calculate displacement

```bash
python examples/qseis2025.py --regional --output-dir examples/output/qseis2025-regional
```

On Windows with Conda, the equivalent non-interactive command is:

```powershell
conda run -n pygrnwang python examples/qseis2025.py --regional --output-dir examples/output/qseis2025-regional
```

The script prepares a local model, writes solver input, installs the 64 s source samples, runs QSEIS2025 and converts the library. It then reads the selected source mechanism with `output_type="disp"` and saves displacement curves. It checks that arrays have the expected components and contain finite values. All generated files stay below the selected output directory; no files from `test/` are required. The solver run takes a few minutes.

## 3. Inspect the result

Look for `disp.png`, `disp.npz`, `source_time_function.json` and `summary.json` in the output directory. The saved displacement array has shape `(3, 3, 256)`: three distances, three components and samples from 0 to 1020 s inclusive. The summary records the environment, calculation time, array dimensions and output size. The underlying library retains its full 1024-sample solver output, geometry metadata and converted binary arrays.

```{figure} _static/examples/qseis2025-regional.png
:alt: QSEIS2025 example displacement traces at 300, 600 and 900 km from a small AK135 Green's-function library.
:width: 100%

Validated QSEIS2025 displacement at 300, 600 and 900 km over 0–1020 s since source origin. Its model, mechanism and numerical choices are shown in the script below.
```

The vector reader uses **east, north, up** when `rotate=True`. Displacement is reported in metres for the moment specified by the script. The plotting time axis must be interpreted with the example's reduction and sampling settings; it is not automatically a P-relative axis. See [scientific conventions](conventions.md) before changing these settings.

The source uses `wavelet_type=0` with 1024 custom moment-rate samples. The helper precompensates QSEIS's numerical damping, so the effective pulse has unit area and a 32 s centroid; see the [regional tutorial](backends/qseis2025.md#regional-waveforms-at-300-600-and-900-km). The library stores rate kernels, and the reader integrates them once for `output_type="disp"`.

## 4. Add strain and stress

Use a separate directory because the output flags change the computed library:

```bash
python examples/qseis2025.py --regional --observables all --output-dir examples/output/qseis2025-regional-tensors
```

The script requests displacement, strain and stress directly from the reader and saves 0–1020 s of each. The saved strain/stress arrays have shape `(3, 6, 256)`. With geographic rotation enabled, symmetric tensors are stored as `[EE, EN, EU, NN, NU, UU]`; U is the same upward vertical component called Z elsewhere in the code. Strain is dimensionless and stress is in pascals for the chosen source moment. These curves are a small workflow example, not a convergence study.

## 5. Reuse a finished example

For the same model, geometry and output flags:

```bash
python examples/qseis2025.py --regional --output-dir examples/output/qseis2025-regional --reuse
```

Use a fresh output directory after changing calculation parameters. Reuse reloads the existing model and data and writes `summary-reuse.json`; it does not prove that existing files match newly selected scientific settings.

## The complete script

The executable script is included directly here, so the documentation and the tested example share one source.

```{literalinclude} ../examples/qseis2025.py
:language: python
:linenos:
```

## Shorter introductory calculation

Without `--regional`, the same script runs a lighter calculation at 30, 60 and 90 km with 0.5 s sampling. It exports 0–100 s, giving 201 samples, from a 256-sample native library and completes in seconds:

```console
python examples/qseis2025.py --output-dir examples/output/qseis2025
```

See the [QSEIS2025 tutorial](backends/qseis2025.md) for its figures and the [backend comparison](guides/backend-comparison.md) for cross-backend results and limitations.

## Continue

- [QSEIS2025 details](backends/qseis2025.md): time functions, observables and numerical settings.
- [Read and process a library](guides/reading.md): interpolation, rotation, filtering and resampling.
- [Choose another backend](backends/index.md): static deformation and spherical Earth calculations.
- [Validation record](validation.md): measured example results and platform coverage.
