# Published comparison data: 2026-09-15

These are measured outputs of the fresh single-depth AK135-FC comparison described in the backend comparison guide. The downloads allow redrawing the documented waveforms and independently recomputing the reported metrics without solving Green functions again. They do not include the multi-gigabyte Green-function libraries.

## Files

- `displacement.png`, `displacement-zoom.png`, `displacement-unfiltered.png`: five dynamic backends with a common spatial point source. The QSEIS06 trace is an isolated point-source control, not the released solver's default.
- `stress-300km.png`, `stress-600km.png`, `stress-900km.png` and their `-zoom.png` counterparts: six ENU stress components from point-source QSEIS2025 and QSSP2020.
- `source-radius-control.png`: default Gaussian smoothing versus point-source controls for both QSEIS versions.
- `model-ak135fc.nd`: the exact six-column AK135-FC input, including Qp/Qs. The spherical solvers used 138 numeric rows; QSEIS used the first 40 rows through 1601.5 km before its Earth-flattening transform.
- `result-summary.json`: settings, source/model/executable checksums, full-window combined metrics, controls, environment, computational cost, and limitations.
- `baseline-metrics.csv/.json`, `point-source-metrics.csv/.json`, `source-radius-metrics.csv/.json`: every component and combined metric. JSON files omit local paths and large private provenance; numerical rows match their original source files.
- `SHA256SUMS.txt`: SHA-256 for every other published file in this directory.

## Array names, shapes, and units

All three NPZ files contain `time_s` (1601,) = 0, 0.25, ..., 400 seconds since earthquake origin and `distances_km` (3,) = [300, 600, 900]. Arrays are float64, loadable with `numpy.load(..., allow_pickle=False)`.

Waveform keys have the form `<processing>_<backend>_<quantity>`. Processing is `raw` or `lowpass_0p4Hz`. Displacement arrays have shape **(3 distances, 3 ENU components, 1601 samples)** in **metres**, component order **[East, North, Up]**. Stress arrays have shape **(3 distances, 6 components, 1601 samples)** in **Pa**, order **[ee, en, eu, nn, nu, uu]**, tension positive; shear entries are tensor components without an extra factor of two.

| Archive | Backend keys | Contents and source convention |
| --- | --- | --- |
| `baseline-waves.npz` | `qseis06`, `qseis2025`, `spgrn2012`, `spgrn2020`, `qssp2020` | `disp` for all; `stress` for QSEIS2025/QSSP2020. QSEIS06 retains its default 0.05 Gaussian smoothing; QSEIS2025 and all spherical backends are point sources. |
| `point-source-waves.npz` | `qseis06`, `qseis2025`, `spgrn2012`, `spgrn2020`, `qssp2020` | `disp` for all; QSEIS06 is the isolated point-source control. No stress keys. |
| `source-radius-waves.npz` | `qseis06_stock`, `qseis2025_gaussian`, `qseis06_point`, `qseis2025_point` | `disp` for both Gaussian and point-source controls. |

The QSEIS06 key intentionally refers to different physical controls in the baseline and point-source archives. Read the table before comparing arrays across archives.

## Processing and metrics

`raw` means no additional presentation filter: the common 1.25 s sin-squared moment-rate STF and continuous damped-Fourier rate integral have already been applied. The QSSP exported-rate indexing has also been deterministically restored using its zero-Nyquist condition. These archives therefore are not untouched native displacement files.

`lowpass_0p4Hz` applies a fourth-order Butterworth filter forward and backward to complete native records before the common grid and comparison window are selected. The two-pass amplitude at the 0.4 Hz critical frequency is one half (-6.02 dB). Do not apply this filter again when redrawing a low-pass key.

The displayed units are micrometres for displacement and Pa for stress. No individual trace normalization, fitted amplitude, or fitted time shift is used in any published plot or metric. True native time origins are retained; only linear interpolation puts traces on the same 0.25 s grid. Zoom windows are 75–115 s, 155–220 s, and 240–310 s for 300, 600, and 900 km. The main reported metrics use 0–400 s; additional `view=zoom` rows in the point/source-radius metric files explicitly identify their shorter intervals.

Relative L2 (%) = 100 * sqrt(sum((u-ref)^2) / sum(ref^2)). This combines timing, amplitude, and waveform differences; it is not an amplitude-error percentage. The printed paper Eq. (14) has no square root, so `paper_energy_misfit = (relative_l2_percent/100)^2`. Correlation is zero-shift Pearson correlation. `all_components` concatenates component arrays without normalization. The reference column must be respected: SPGRN2020 for displacement, QSSP2020 for stress, and the specifically named QSEIS control in source-radius comparisons. The reference is not asserted to be an exact solution.

## Isolated native controls

`qseis2025-sharded-source.diff` and `qseis06-point-sharded-source.diff` show the exact reviewed source changes with repository-relative paths. Both raise the sublayer capacity to 2048, partition the frequency loop without changing the common FFT grid, and preserve more output digits for summing partitions. The QSEIS06 control additionally changes its hard-coded Gaussian radius ratio from 0.05 to 0. These diffs were used in isolated executables only; they were not applied to installed/released solvers or repository Fortran sources. They are provenance for this calculation, not instructions to change the package's default behavior. They require prepared native inputs and explicit frequency-shard settings; they are not a standalone rebuild recipe.

The full 3×3×2048 displacement arrays are bit-identical between QSEIS06 and QSEIS2025 with the same Gaussian setting, and also between their point-source controls. These checks apply to this particular experiment. Residual QSEIS/spherical differences remain, and geometry, deep truncation, discretization, and integration convergence have not been isolated as individual causes.

The comparison follows the presentation of [Zhou et al. (2026), DynCFS](https://doi.org/10.1093/gji/ggaf534), Fig. 3, extending its 10 km near-field stress example to regional distances and dynamic displacement. It is not a reproduction of that figure. The validated platform was Windows 11; the environment and timing limits are recorded in the summary.
