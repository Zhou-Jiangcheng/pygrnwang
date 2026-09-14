"""Build directly from source without invoking package build hooks."""
from pathlib import Path
import os
import sys
import tomllib

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
os.environ.setdefault("MPLBACKEND", "Agg")
metadata = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]
project = "pygrnwang"
author = "Jiangcheng Zhou"
copyright = "2026, Jiangcheng Zhou"
release = metadata["version"]
version = release
language = "en"
root_doc = "index"
extensions = [
    "myst_parser", "sphinx.ext.autodoc", "sphinx.ext.autosummary",
    "sphinx.ext.napoleon", "sphinx.ext.mathjax", "sphinx.ext.viewcode",
    "sphinx.ext.githubpages", "sphinx_design",
]
myst_enable_extensions = ["colon_fence", "dollarmath", "amsmath", "deflist"]
myst_heading_anchors = 3
exclude_patterns = ["_build", "requirements*.txt", "requirements.in", "README.md"]
autosummary_generate = True
autodoc_member_order = "bysource"
autodoc_typehints = "description"
napoleon_numpy_docstring = True
napoleon_google_docstring = False
napoleon_use_param = True
napoleon_use_rtype = False
html_theme = "pydata_sphinx_theme"
html_title = f"pygrnwang {release}"
html_short_title = "pygrnwang"
html_baseurl = "https://zhou-jiangcheng.github.io/pygrnwang/"
html_static_path = ["_static"]
html_css_files = ["custom.css"]
html_theme_options = {
    "github_url": "https://github.com/Zhou-Jiangcheng/pygrnwang",
    "show_nav_level": 1,
    "navigation_depth": 3,
    "show_toc_level": 2,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navbar_persistent": ["search-button"],
    "collapse_navigation": True,
    "announcement": f"Documentation for the main branch · package {release}",
}
html_context = {
    "github_user": "Zhou-Jiangcheng", "github_repo": "pygrnwang",
    "github_version": "main", "doc_path": "docs",
}
html_show_sourcelink = True
html_last_updated_fmt = None
linkcheck_anchors = False
linkcheck_timeout = 20
linkcheck_retries = 2
