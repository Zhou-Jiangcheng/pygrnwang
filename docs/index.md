# Green's functions, from model to waveform

```{raw} html
<p class="hero-kicker">pygrnwang · computational seismology</p>
<p class="hero-copy">Build Green's-function libraries with Wang's Fortran solvers, synthesize seismograms and static deformation, and keep model, source, time and coordinate conventions explicit.</p>
```

::::{grid} 1 1 3 3
:gutter: 3

:::{grid-item-card} Start with QSEIS2025
:link: quickstart
:link-type: doc

Run a small, complete example: prepare AK135, compute a library, read three-component displacement and save a figure.
:::

:::{grid-item-card} Choose a backend
:link: backends/index
:link-type: doc

Compare static and dynamic workflows, spherical and layered models, and supported observables.
:::

:::{grid-item-card} 中文入门
:link: zh/index
:link-type: doc

安装程序、运行同一套 QSEIS2025 示例，并理解输出文件和坐标约定。
:::
::::

## A complete workflow

1. **Prepare the model.** Define layered velocities, density, attenuation and the source/receiver geometry.
2. **Compute the library.** Preprocess the solver input, run the backend and convert its output where required.
3. **Synthesize observables.** Select a source mechanism, distance, sampling and output quantity.
4. **Interpret the result.** Check component order, amplitude normalization and the time origin before comparing observations.

The [scientific conventions](conventions.md) are part of the interface. In particular, vector and tensor orders, time reduction and source-duration units vary across backends.

## What is included

The package wraps **EDGRN/EDCMP**, **QSEIS06**, **QSEIS2025**, **SPGRN2012**, **SPGRN2020** and **QSSP2020**. Python provides preprocessing, serial and parallel orchestration, library lookup, source synthesis, coordinate transforms and signal processing. TauP travel times run through Java subprocesses or an ObsPy fallback.

Python **3.9 or newer** is supported. Binary wheels cover the platforms listed in the [installation guide](installation.md). This site follows `main`; the package version appears in the header.

```{toctree}
:maxdepth: 2
:caption: Getting started

installation
quickstart
zh/index
```

```{toctree}
:maxdepth: 2
:caption: User guide

backends/index
conventions
guides/index
cli
```

```{toctree}
:maxdepth: 2
:caption: Reference and development

api/index
validation
development
changelog
references
```
