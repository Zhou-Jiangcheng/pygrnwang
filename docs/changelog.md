# Changelog

## Unreleased documentation

- Recalculated SPGRN2020 with its complete-wavefield option and QSSP2020 with
  harmonic settings 2000/8000 after identifying insufficient low-frequency
  cutoffs in the initial examples. Added convergence evidence, corrected
  harmonic-parameter documentation and plotted SPGRN2020 on its native
  source-origin time axis.
- Added QSEIS06/QSEIS2025 `--regional` calculations at 300/600/900 km, with a
  4092 s native window and 0–1020 s exports. The default QSEIS2025 introduction
  retains its 0–100 s crop. Added a reproducible five-backend comparison and
  included the regional calculations in tutorial CI.
- Matched regional QSEIS's effective source time function to SPGRN2020 using
  1024 custom moment-rate nodes with numerical damping precompensated. Kept
  the source's unit area and 32 s centroid, read rate kernels explicitly, and
  integrated them once in the examples before updating the comparison.
  Recorded that the remaining QSEIS/SPGRN2020 differences increased to about
  12.7–23.0%, ruling out the previous STF mismatch as their main explanation.

- Marked QSEIS06 and SPGRN2012 as deprecated in pygrnwang. New calculations should use QSEIS2025 and SPGRN2020 respectively; existing interfaces and tutorials remain available. Migration requires rebuilding libraries and validating the replacement backend's settings and output conventions.

- Added an English user guide and API reference, with Chinese installation and QSEIS2025 quickstart pages.
- Added small executable tutorials for QSEIS2025, QSEIS06, SPGRN2012, SPGRN2020, QSSP2020 and EDGRN/EDCMP.
- Cropped the default QSEIS2025 introductory displacement, strain and stress arrays and figures to 0–100 s inclusive (201 samples at 0.5 s). The native 127.5 s, 256-sample Green's library and other backend examples retain their original windows.
- Documented coordinate/component order, source normalization, time reduction, waveform units and solver-specific limitations.
- Added strict documentation builds, example validation and GitHub Pages deployment configuration.

This documentation change preserves calculation signatures, return types and numerical code. It records existing behavior and known limitations.

## Current 3.0.0 development baseline

The baseline used for this documentation declares Python 3.9+ and setuptools 77+. Travel-time calculations use a lazily compiled Java subprocess bridge, with ObsPy as the automatic alternative when a JDK is unavailable. JPype is no longer a dependency. Wheels contain a package-local TauP.jar and install a copy into the environment's scripts directory.

These entries describe the checked-out source baseline; they do not assert a new PyPI release date. For earlier releases, consult the [repository history](https://github.com/Zhou-Jiangcheng/pygrnwang/commits/main/).
