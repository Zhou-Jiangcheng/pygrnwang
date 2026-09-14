# Development and documentation maintenance

## Work on the package

Use the [source installation](installation.md) to create an activated environment with gfortran, then run `python -m pip install -e .`. Fortran changes require a rebuild. Keep a clean output directory for each numerical configuration.

Run the focused TauP regression suite with:

```bash
python -m unittest discover -s tests -p test_pytaup_subprocess.py -v
```

The test suite checks Java subprocess queries, ObsPy fallback and installed-wheel JAR placement. It is not a numerical validation of every solver. The `examples/` scripts provide small end-to-end calculations; the [validation record](validation.md) describes what was actually run. Multi-node MPI requires a separate cluster test.

## Build the documentation

The documentation environment uses Python 3.12 independently of the package's Python 3.9 support floor. Install only the locked documentation dependencies; Sphinx imports this checkout directly and does not run setup.py.

```bash
python3.12 -m venv .venv
```

Activate `.venv` using the command appropriate to your shell, then:

```bash
python -m pip install -r docs/requirements.txt
python docs/api/check_coverage.py
python -m sphinx -b html -W --keep-going docs docs/_build/html
python -m http.server 8000 --directory docs/_build/html
```

On Windows, create the environment with an available Python 3.12 interpreter (for example `py -3.12 -m venv .venv`), and activate it with `.venv\Scripts\Activate.ps1`. Open `http://localhost:8000/` after starting the server.

HTML generation uses already checked-in example figures. It does not execute examples, start Java or compile the Fortran programs. Core modules are imported normally rather than replaced with mock objects.

`docs/requirements.in` describes the direct constraints. Refresh the lock deliberately using Python 3.12 resolution, review the changes, and rebuild before committing:

```bash
uv pip compile docs/requirements.in --python-version 3.12 --universal --no-header --output-file docs/requirements.txt
```

An optional external-link check is separate because publisher sites and other external servers can be temporarily unavailable:

```bash
python -m sphinx -b linkcheck docs docs/_build/linkcheck
```

## Add or change an interface

1. Update the English NumPy-style docstring: purpose, parameter meanings and units, defaults, output shape/component order, side effects and errors.
2. If it is a supported public interface, add it to the API inventory and the corresponding reference page. Internal helpers remain in the advanced reference.
3. Explain changed scientific behavior in the relevant guide and update a small executable example when it illustrates that behavior.
4. Rebuild with warnings treated as errors and run the relevant example/regression check. Update the English and Chinese quickstart together if its commands or output change.
5. Record the user-visible change in the changelog.

Do not infer units from a variable name alone. Verify Python synthesis/rotation, the input template and the Fortran output definition together. In particular, unrotated tensors and wavelet durations differ between backends. Existing numerical defects are recorded in [troubleshooting](guides/troubleshooting.md); documentation work must not silently change numerical behavior.

## Example conventions

Every backend script supports `--output-dir` and `--reuse`, includes a main guard and uses a model prepared from repository data. QSEIS2025 additionally accepts `--observables all`. The common helper validates arrays and writes machine-readable output and a run summary. Keep tutorial grids small. Generated solver libraries and logs are excluded from Git; the selected figures and validation report are reviewed artifacts.

## GitHub Pages

The documentation workflow builds on pull requests and uploads HTML for review. A successful build on `main` deploys to:

<https://zhou-jiangcheng.github.io/pygrnwang/>

Repository settings must select **Settings → Pages → Build and deployment → Source: GitHub Actions**. The deployment job uses the `github-pages` environment, `pages: write` and `id-token: write`. Fork pull requests build with read permissions and cannot deploy. See [GitHub's custom-workflow documentation](https://docs.github.com/en/pages/getting-started-with-github-pages/using-custom-workflows-with-github-pages).

The first site follows `main` and shows the package version. It does not archive historical versions. The package release workflow remains tag-triggered; publishing documentation does not create a PyPI release.

## Package release checklist

Before tagging, check the version in package metadata and `pygrnwang.__version__`, the changelog and dependency support matrix. Build source and wheel artifacts, inspect packaged executables/JARs, and test installation outside the checkout. The existing publishing workflow uses Trusted Publishing; the repository/workflow identity must match the PyPI configuration. Add release notes describing behavior and validation rather than assuming a successful import proves a solver works.
