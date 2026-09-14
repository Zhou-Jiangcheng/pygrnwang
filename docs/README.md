# Documentation development

The published site uses Sphinx, MyST Markdown and the PyData theme.

Use Python 3.12 and install `docs/requirements.txt` in an isolated environment.
From the repository root:

    python docs/api/check_coverage.py
    python -m sphinx -b html -W --keep-going docs docs/_build/html
    python -m http.server 8000 --directory docs/_build/html

Full maintenance instructions are in `development.md`. The package is imported
from source: building these pages does not invoke setup.py, compile Fortran,
start a JVM or run the examples. Tutorial figures are checked-in outputs from
separately validated example scripts.
