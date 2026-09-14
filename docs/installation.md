# Installation

## Choose an installation route

| Route | Python | Fortran compiler | Java | MPI |
|---|---|---|---|---|
| Install a matching wheel | 3.9+ | Not needed for pygrnwang | Optional JDK for TauP | Optional, for multi-node workflows |
| Build from source | 3.9+ | `gfortran` on PATH | Optional JDK for TauP | Optional |
| Build this documentation | 3.12 | Not needed | Not needed | Not needed |

The release workflow targets **Linux x86-64**, **Windows x86-64** and **macOS Apple Silicon (arm64)**, with CPython 3.9–3.14 selected subject to interpreter and dependency availability. A platform in the build configuration does not guarantee that every release has every wheel. Intel macOS and other architectures need a source build unless a matching artifact is provided.

## Install from PyPI

Create and activate an isolated environment, then install:

```bash
python -m pip install --upgrade pip
python -m pip install pygrnwang
python -c "import pygrnwang; print(pygrnwang.__version__)"
```

To require a wheel for pygrnwang and fail immediately if none is available, add `--only-binary=pygrnwang`. Otherwise pip may attempt a source build, which needs `gfortran`.

The Python dependencies are NumPy, SciPy, pandas, ObsPy and tqdm. Matplotlib is used by the example plots. MPI needs both `mpi4py` and a working MPI runtime.

```{important}
On Unix, the current Python bulk runners look for `<solver>.bin` in the
environment's `bin` directory, while a standard wheel stores native binaries
inside `pygrnwang/exec`. The command-line wrappers can find the package copy,
but bulk creation may fail with a missing executable. Use the editable source
installation below for these calculation workflows. This existing packaging
limitation is documented rather than changing runtime behavior in a docs update.
```

## Run the examples from this repository

The example scripts live in the Git repository. Check out the revision matching the package you want to use. To use the current documentation and source together:

```bash
git clone https://github.com/Zhou-Jiangcheng/pygrnwang.git
cd pygrnwang
conda create -n pygrnwang -c conda-forge python=3.12 numpy scipy pandas obspy tqdm matplotlib gfortran
conda activate pygrnwang
python -m pip install -e .
```

`setuptools>=77` is required for source builds; pip installs the build requirements in its isolated build environment. The setup hook compiles seven Fortran executables. An editable installation makes Python source changes visible without reinstalling; changes to Fortran still require rebuilding.

On Windows, keep the Conda environment activated when running Python. For scripts and automation, use:

```powershell
conda run -n pygrnwang python examples/qseis2025.py
```

The active environment supplies compiler and numerical-library DLL directories. A bare path to an inactive Conda environment's `python.exe` does not provide this setup.

## Source-build prerequisites by platform

- **Linux:** install GCC/gfortran through the system package manager or Conda. The current setup requests static linking; the corresponding static runtime libraries must be available.
- **Windows:** use the activated Conda environment above or a compatible MinGW-w64 toolchain. Confirm `gfortran --version` works in the same terminal as the build command.
- **macOS:** install GCC/gfortran through Homebrew or Conda. Build for your machine's architecture and ensure the resulting Fortran runtime libraries are available. CI repairs runtime dependencies for published wheels.

Compilation flags can be extended through `PYGRNWANG_FFLAGS`. Prefer the release wheels when available; compiler versions and static-library availability affect source builds.

## Optional Java travel times

Install a **JDK** and make both `java` and `javac` available on PATH:

```bash
java -version
javac -version
```

The bundled TauP version is 2.6.1. Its small bridge is compiled on the first query and cached for the current Python process. JPype is not used. If the JDK or JAR cannot be found at module import, the general arrival-time functions use ObsPy. An explicitly requested `taup_time_java` call requires the Java backend.

Wheels retain `TauP.jar` in the package's `exec` directory and install another copy in the environment's `Scripts` directory on Windows or `bin` on Unix. Lookup prefers the package copy. See the [TauP guide](guides/taup.md) for custom models and batching.

## Optional MPI

For multi-node work, install `mpi4py` against your cluster's MPI implementation. Verify the same environment, solver binaries and shared paths on all nodes. Use the site-provided job launcher; a workstation installation of `mpi4py` alone does not configure a cluster. See the [parallel execution guide](guides/parallel.md).

## Confirm the installation

The `pygrnwang` command prints an installation message. It does not run the seven solvers. Use the [quickstart](quickstart.md) to validate an actual calculation and consult [troubleshooting](guides/troubleshooting.md) when a binary, compiler or model is missing.
