# Choose a backend

Choose the physical formulation and observable before choosing numerical
settings. Every tutorial below contains a complete build/read/plot script.

| Backend | Formulation and use | Main reader outputs |
| --- | --- | --- |
| [QSEIS2025](qseis2025.md) | Layered half-space; dynamic waveforms with direct tensor/rotation outputs | Displacement, velocity, acceleration, strain/stress and rates, rotation and rate, volume |
| [QSEIS06](qseis06.md) | Layered half-space; established vector waveform workflow | Displacement, velocity, acceleration; separate derivative-library tools |
| [SPGRN2012](spgrn2012.md) | Spherical layered Earth; reduced-time waveform library | Displacement, velocity, acceleration |
| [SPGRN2020](spgrn2020.md) | Spherical layered Earth; windows tied to direct P onset | Displacement, velocity, acceleration |
| [QSSP2020](qssp2020.md) | Spherical layered Earth with optional self-gravitation and broad observables | Vector, tensor, rotation and gravity families; see conversion limitation |
| [EDGRN → EDCMP](edgrn_edcmp.md) | Layered elastic half-space; static dislocation response | Displacement, strain, stress and tilt |

QSEIS2025 is the recommended first tutorial. Older backends retain their
own workflows and are not marked as deprecated. Selection depends on
physics and validated parameter choices, not simply the newest year.

## Shared tutorial setup

From a checkout root, install the package and plotting dependency, then run
a tutorial with the active Python environment:

```console
python -m pip install -e .
python -m pip install matplotlib
python examples/qseis2025.py
```

Use [installation](../installation.md) for wheels and platform requirements.
A source/editable installation requires gfortran. The scripts require the
compiled executables from the package and write only to their output
directories.

All scripts support `--output-dir PATH` and `--reuse`. The default is
`examples/output/<backend>`. A successful run writes `disp.npz`,
`disp.png`, `summary.json` and a `library/` directory; QSEIS2025's tensor
option adds strain/stress arrays and figures. Scripts check shape,
finiteness and nonzero response. The common helpers are documented in
[model preparation](../guides/models.md).

Tutorial physics is explicit: a 10 km source, surface receiver, 30° azimuth,
strike/dip/rake 30°/45°/90°, and moment `10^15 N m`. Elastic velocities and
density come from the bundled AK135 model; attenuation is the illustrative
constant `Qp=600, Qs=300`, not AK135-F attenuation. QSEIS/EDGRN use the first
24 numeric model rows (to 809.5 km); spherical examples use the full Earth.
EDGRN additionally builds 11 km depth to satisfy its two-depth minimum.

The small frequency bands and coarse grids demonstrate computation and
data interpretation. They are not numerical-convergence evidence for a
research application. Consult [scientific conventions](../conventions.md)
before comparing backend amplitudes, tensor components or time origins.

```{toctree}
:maxdepth: 1

qseis2025
qseis06
spgrn2012
spgrn2020
qssp2020
edgrn_edcmp
```
