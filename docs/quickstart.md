# Quickstart: QSEIS2025

This example computes a small Green's-function library, synthesizes displacement and writes a figure. It uses the AK135 model content included with pygrnwang, one source depth, one receiver depth and three distances. Computation is serial. The source is 10 km deep, receivers are at the surface and distances are 30, 60 and 90 km. The native library uses a 0.5 s sampling interval and a 127.5 s window containing 256 samples. After synthesis, the example crops every saved waveform and plot to 0–100 s inclusive, giving 201 samples without changing the underlying library. The layered calculation uses the first 24 numeric model rows, down to 809.5 km. Constant Qp=600 and Qs=300 are illustrative tutorial choices, not the full AK135-F attenuation model. The mechanism is strike 30°, dip 45°, rake 90°, scaled to M0 = 10^15 N m.

## 1. Prepare the environment

Follow [installation](installation.md), including cloning the repository to obtain `examples/`. Run the following commands from the repository root in the activated environment.

## 2. Calculate displacement

```bash
python examples/qseis2025.py --output-dir examples/output/qseis2025
```

On Windows with Conda, the equivalent non-interactive command is:

```powershell
conda run -n pygrnwang python examples/qseis2025.py --output-dir examples/output/qseis2025
```

The script prepares a local model, writes solver input, runs QSEIS2025, converts the library, reads the selected source mechanism and saves displacement curves. It checks that arrays have the expected components and contain finite values. All generated files stay below the selected output directory; no files from `test/` are required.

## 3. Inspect the result

Look for `disp.png`, `disp.npz` and `summary.json` in the output directory. The saved displacement array has shape `(3, 3, 201)`: three distances, three components and samples from 0 to 100 s inclusive. The summary records the environment, calculation time, array dimensions and output size. The underlying library retains its full 256-sample solver output, geometry metadata and converted binary arrays.

```{figure} _static/examples/qseis2025.png
:alt: QSEIS2025 example displacement traces from a small AK135 Green's-function library.
:width: 100%

Validated QSEIS2025 displacement, cropped to 0–100 s since source origin. Its model, mechanism and numerical choices are shown in the script below.
```

The vector reader uses **east, north, up** when `rotate=True`. Displacement is reported in metres for the moment specified by the script. The plotting time axis must be interpreted with the example's reduction and sampling settings; it is not automatically a P-relative axis. See [scientific conventions](conventions.md) before changing these settings.

## 4. Add strain and stress

Use a separate directory because the output flags change the computed library:

```bash
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025-all
```

The script reads displacement, strain and stress, then crops all three to 0–100 s inclusive before saving arrays and figures. The saved strain/stress arrays have shape `(3, 6, 201)`. With geographic rotation enabled, symmetric tensors are stored as `[EE, EN, EU, NN, NU, UU]`; U is the same upward vertical component called Z elsewhere in the code. Strain is dimensionless and stress is in pascals for the chosen source moment. These curves are a small workflow example, not a convergence study.

## 5. Reuse a finished example

For the same model, geometry and output flags:

```bash
python examples/qseis2025.py --output-dir examples/output/qseis2025 --reuse
```

Use a fresh output directory after changing calculation parameters. Reuse reloads the existing model and data and writes `summary-reuse.json`; it does not prove that existing files match newly selected scientific settings.

## The complete script

The executable script is included directly here, so the documentation and the tested example share one source.

```{literalinclude} ../examples/qseis2025.py
:language: python
:linenos:
```

## Longer-distance calculation

Keep the introductory 100 s output above for 30/60/90 km. For 300/600/900 km,
where major arrivals extend beyond that window, run:

```console
python examples/qseis2025.py --regional --observables all
```

This separate mode uses 4 s sampling, a 4092 s native window, a damping-compensated 64 s source
and the flat-Earth transformation. It exports 0–1020 s to
`examples/output/qseis2025-regional/`. See the
[regional tutorial](backends/qseis2025.md#regional-waveforms-at-300-600-and-900-km)
and [backend comparison](guides/backend-comparison.md) for results and limitations.

## Continue

- [QSEIS2025 details](backends/qseis2025.md): time functions, observables and numerical settings.
- [Read and process a library](guides/reading.md): interpolation, rotation, filtering and resampling.
- [Choose another backend](backends/index.md): static deformation and spherical Earth calculations.
- [Validation record](validation.md): measured example results and platform coverage.
