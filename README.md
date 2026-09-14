# pygrnwang

**[Documentation](https://zhou-jiangcheng.github.io/pygrnwang/)** · [Quickstart](https://zhou-jiangcheng.github.io/pygrnwang/quickstart.html) · [中文入门](https://zhou-jiangcheng.github.io/pygrnwang/zh/index.html) · [Examples](https://github.com/Zhou-Jiangcheng/pygrnwang/tree/main/examples)

This Python package serves as the frontend for calculating and building a Green's function library for synthetic seismograms. The backend consists of Wang Rongjiang's program for calculating synthetic seismograms, including EDGRN/EDCMP, [QSEIS_STRESS](https://github.com/Zhou-Jiangcheng/QSEIS_2006_STRESS), SPGRN, and QSSP (Wang, 1999; Wang 2003; Wang and Wang 2007; Wang et al., 2017). The code includes two parallel modes: one using the multiprocessing library (single-node multi-process) and the other using MPI (multi-node).

The **QSEIS06** and **SPGRN2012** backends are **deprecated** in pygrnwang.
For new calculations, use [QSEIS2025](https://zhou-jiangcheng.github.io/pygrnwang/backends/qseis2025.html)
and [SPGRN2020](https://zhou-jiangcheng.github.io/pygrnwang/backends/spgrn2020.html), respectively.
Existing interfaces remain available; migration requires rebuilding libraries
and validating the replacement backend's settings and output conventions.

# Installation

Python 3.9 or later is supported. Building from source requires gfortran and
setuptools >=77; installing a compatible wheel does not require gfortran.

Travel times use TauP 2.6.1 through a Java subprocess when a JDK (`java` and
`javac` on PATH) is available, and otherwise fall back to ObsPy. JPype is not
required. The Java bridge is compiled on the first query and cached for the
current Python process; importing the library does not start Java. Travel-time
tables send all distances to one Java process. Individual queries start one
Java process per call, so use `create_tpts_table` for large distance grids.

Wheels include `TauP.jar` inside the package and install a second copy in the
target environment's `Scripts` (Windows) or `bin` directory. Resource lookup
prefers the package copy and falls back to the environment copy.

1. For user mode

```
pip install pygrnwang
```

On Linux and macOS, the current bulk runners require an editable source
installation because of an existing executable-path limitation. Follow the
[installation guide](https://zhou-jiangcheng.github.io/pygrnwang/installation.html)
for the complete platform instructions.

2. For developer mode

Prepare and activate the environment and platform-specific compiler described
in the [installation guide](https://zhou-jiangcheng.github.io/pygrnwang/installation.html), then:

```bash
git clone https://github.com/Zhou-Jiangcheng/pygrnwang.git
cd pygrnwang
pip install -e .
```

# Usage

After cloning the repository and installing the current source, run the small
QSEIS2025 example from the repository root:

```bash
python examples/qseis2025.py --output-dir examples/output/qseis2025
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025-all
```

The scripts prepare a model from bundled AK135 elastic parameters, compute three
distances, read displacement (and optionally strain/stress), and save NPZ files,
figures and a JSON run summary. Constant attenuation and small numerical grids
are explicit tutorial choices, not convergence recommendations for research.

Other complete workflows are listed in [the example scripts](https://github.com/Zhou-Jiangcheng/pygrnwang/tree/main/examples).
Read the [scientific conventions](https://zhou-jiangcheng.github.io/pygrnwang/conventions.html) before interpreting
amplitudes, component order or time origins. The [validation record](https://zhou-jiangcheng.github.io/pygrnwang/validation.html)
records measured results and platform coverage.

# Documentation development

Use an isolated Python 3.12 environment, then run:

```bash
python -m pip install -r docs/requirements.txt
python docs/api/check_coverage.py
python -m sphinx -b html -W --keep-going docs docs/_build/html
```

The build imports source directly and does not compile Fortran or run numerical
tutorials. See [development instructions](https://zhou-jiangcheng.github.io/pygrnwang/development.html) for previewing the
site, maintaining APIs and GitHub Pages deployment.

# Reference

Wang, R. (1999). A simple orthonormalization method for stable and efficient computation of Green’s functions.  *Bulletin of the Seismological Society of America* ,  *89* (3), 733–741. https://doi.org/10.1785/BSSA0890030733

Wang, R. (2003). Computation of deformation induced by earthquakes in a multi-layered elastic crust—FORTRAN programs EDGRN/EDCMP. Computers & Geosciences, 29(2), 195–207. https://doi.org/10.1016/S0098-3004(02)00111-5

Wang, R., & Wang, H. (2007). A fast converging and anti-aliasing algorithm for green’s functions in terms of spherical or cylindrical harmonics. Geophysical Journal International, 170(1), 239–248. https://doi.org/10.1111/j.1365-246X.2007.03385.x

Wang, R., Heimann, S., Zhang, Y., Wang, H., & Dahm, T. (2017). Complete synthetic seismograms based on a spherical self-gravitating earth model with an atmosphere–ocean–mantle–core structure. Geophysical Journal International, 210(3), 1739–1764. https://doi.org/10.1093/gji/ggx259

Zhou, J., Wang, R., & Zhang, Y. (2026). DynCFS: a program for modeling dynamic coulomb failure stress changes in layered elastic media. Geophysical Journal International, ggaf534. https://doi.org/10.1093/gji/ggaf534