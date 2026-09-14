# QSEIS2025: displacement, strain and stress

QSEIS2025 calculates dynamic Green's functions in a layered half-space and
provides direct strain/stress and rotation kernels alongside vectors. Use
this tutorial for a complete first calculation and for the introductory
tensor workflow.

## Build and inspect displacement

From the repository root:

```console
python examples/qseis2025.py
```

The script prepares the model and input files, computes the library
sequentially, reads three distances and saves an ENU displacement figure.
E/N/U in example labels is the same east/north/up convention called ENZ
by the API.

The source depth is 10 km, receiver depth 0 km, and distances are 30, 60 and
90 km. The 0.5 s sample interval and 127.5 s window produce 256 native samples.
After synthesis, every saved example waveform and figure is cropped to
0–100 s inclusive (201 samples); the library keeps its full native window.
The model uses the first 24 numeric rows and disables the flat-Earth
transformation for this small half-space example. Moment is `10^15 N m`.

```{literalinclude} ../../examples/qseis2025.py
:language: python
:caption: Complete QSEIS2025 tutorial
```

```{figure} ../_static/examples/qseis2025.png
:alt: Three-component QSEIS2025 displacement at three epicentral distances.

Example displacement in metres with a fixed 0–100 s horizontal axis. Time is
measured from source origin because this example uses zero time reduction.
```

Results are under `examples/output/qseis2025/`. `disp.npz` contains the
three distance traces, component labels, seconds and units; its waveform
array has shape `(3, 3, 201)` and includes both 0 s and 100 s.
`summary.json` records checks and run details. The solver library still
contains 256 samples per trace.
Use `python examples/qseis2025.py --reuse` to read the completed library
again without launching the solver.

## Add tensor outputs

Use a separate output directory so the selected quantities match the library:

```console
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025_tensors
```

This enables displacement, strain and stress; it writes `strain.npz/png`
and `stress.npz/png` in addition to displacement. The tensor arrays have
shape `(3, 6, 201)` and components EE, EN, EU, NN, NU, UU. All three
observables are cropped to 0–100 s inclusive before saving arrays and plots.
Strain is dimensionless and stress is in Pa after the script's explicit
moment scaling.

```{figure} ../_static/examples/qseis2025-strain.png
:alt: Six ENU strain tensor components from the QSEIS2025 tensor tutorial.

Dimensionless tensor strain at 30, 60 and 90 km, cropped to 0–100 s since
source origin; shear entries are tensor strains rather than engineering shear.
```

```{figure} ../_static/examples/qseis2025-stress.png
:alt: Six ENU stress tensor components from the QSEIS2025 tensor tutorial.

Stress in Pa for the same moment and geometry, cropped to 0–100 s since source
origin. The surface receiver's traction components can be zero under the
free-surface boundary condition.
```

The five `output_observables` positions are:

| Zero-based position | Family | Reader requests |
| --- | --- | --- |
| 0 | Vector motion | `disp`, `velo`, `acce` |
| 1 | Fractional volume change | `volume` |
| 2 | Strain | `strain`, `strain_rate` |
| 3 | Stress | `stress`, `stress_rate` |
| 4 | Rotation | `rota`, `rota_rate` |

Thus `[1, 0, 1, 1, 0]` enables the demonstrated tensor calculation.
The physical stored quantity depends on `wavelet_type`. Here type 2 is
a tapered Heaviside, so the stored non-rate kernels include displacement,
strain and stress. `wavelet_duration=4` means **four samples**, or 2 s,
not four seconds.

## Parameters that control the calculation

`N_each_group` splits receiver distances into jobs; it is 3 in this example.
It controls job size and file layout, not the source model.
`wavenumber_sampling_rate`, `anti_alias`, `slowness_int_algorithm` and an
optional `slowness_window` control integration. QSEIS2025 additionally
exposes `eps_estimate_wavenumber` and `source_radius_ratio`.
Establish convergence when changing band, distance or depth.

`time_reduction_velo` is km/s; zero means no reduction. The raw
`free_surface` switch is the backend partial-solution switch, with 0
retaining free-surface effects. It is not the same Python boolean interface
as QSEIS06's `free_surface`.

The reader supports nearest or trilinear interpolation, ENZ rotation,
optional arrival adjustments and filtering/resampling. See
[reading](../guides/reading.md). An output family's source files must exist
before its derived rate/non-rate version can be requested.

## Outputs and limitations

Each depth pair contains distance-group directories with `grn.inp` and
basis files prefixed `ex`, `ss`, `ds` and `cl`. Components `e*` are strain
and `s*` are stress, despite swapped labels in some template comments.
The example retains ASCII and creates reader-compatible binary files.
Keep `green_lib_info.json` and travel-time tables with the data.

At `rotate=False`, six-component tensor output is a north-reference ENZ
basis, not the vector's RTZ convention. Use `rotate=True` and the documented
[component order](../conventions.md#six-component-tensors) for comparison.

The tutorial verifies executable output and reader behavior. A 30 km distance
grid and its integration defaults do not validate interpolation, high
frequencies, zero-distance behavior or interface stresses for your study.
