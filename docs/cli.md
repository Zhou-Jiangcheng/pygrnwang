# Command-line entry points

The package registers the following commands during installation:

| Command | Behavior |
|---|---|
| `pygrnwang` | Print an installation-success message |
| `edgrn2` | Launch the packaged EDGRN 2 executable |
| `edcmp2` | Launch the packaged EDCMP 2 executable |
| `qseis06` | Launch the packaged QSEIS06 executable |
| `qseis2025` | Launch the packaged QSEIS2025 executable |
| `spgrn2012` | Launch the packaged SPGRN2012 executable |
| `spgrn2020` | Launch the packaged SPGRN2020 executable |
| `qssp2020` | Launch the packaged QSSP2020 executable |

The solver commands are thin wrappers. They pass command-line arguments through, inherit the current working directory and standard input/output, and return the solver's exit status. They do not provide a common argparse interface or generate input files. Do not assume all of them implement `--help` or accept an input filename as a positional argument.

The Fortran programs prompt for their input filename on standard input. From the directory containing a prepared input file, launch the relevant command and enter that filename at its prompt. Relative paths inside the input file are interpreted by the solver from its working directory. The Python preprocessing and execution functions manage these paths for the supported workflows.

On Windows the wrappers add the active environment's `Library/bin` to PATH if present. On Unix they ensure the packaged binary is executable. A missing binary produces an error and exit status 1; verify a matching wheel or successful source compilation.

## Executable tutorials

The example scripts expose their own documented options:

```bash
python examples/qseis2025.py --help
python examples/qseis2025.py --output-dir examples/output/qseis2025
python examples/qseis2025.py --observables all --output-dir examples/output/qseis2025-all
```

All six backend scripts accept `--output-dir` and `--reuse`. `--observables` belongs to the QSEIS2025 tutorial. Start with these scripts or the [Python API](api/index.md) for reproducible library creation.
