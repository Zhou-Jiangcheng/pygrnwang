# Changelog

## Unreleased documentation

- Added an English user guide and API reference, with Chinese installation and QSEIS2025 quickstart pages.
- Added small executable tutorials for QSEIS2025, QSEIS06, SPGRN2012, SPGRN2020, QSSP2020 and EDGRN/EDCMP.
- Documented coordinate/component order, source normalization, time reduction, waveform units and solver-specific limitations.
- Added strict documentation builds, example validation and GitHub Pages deployment configuration.

This documentation change preserves calculation signatures, return types and numerical code. It records existing behavior and known limitations.

## Current 3.0.0 development baseline

The baseline used for this documentation declares Python 3.9+ and setuptools 77+. Travel-time calculations use a lazily compiled Java subprocess bridge, with ObsPy as the automatic alternative when a JDK is unavailable. JPype is no longer a dependency. Wheels contain a package-local TauP.jar and install a copy into the environment's scripts directory.

These entries describe the checked-out source baseline; they do not assert a new PyPI release date. For earlier releases, consult the [repository history](https://github.com/Zhou-Jiangcheng/pygrnwang/commits/main/).
