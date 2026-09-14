# Reading, synthesis and caching

Most dynamic readers share geometry, sampling, filtering and time-window arguments. Their tensor layouts and native time origins differ. Read each function's Notes before combining results.

With `only_seismograms=False`, the dynamic `seek_*` functions return **seven values**. `first_p` and `first_s` remain `None` unless `shift=True`; nearest-grid metadata can differ from the geometry used for interpolation. EDCMP readers return static arrays and normalize the mechanism magnitude.

## read_syn

```{eval-rst}
.. autofunction:: pygrnwang.read_syn.read_syn
```

## read_qseis2025

```{eval-rst}
.. autofunction:: pygrnwang.read_qseis2025.seek_qseis2025
```

## read_qseis06

**Deprecated backend: QSEIS06.** Use the QSEIS2025 reader with a newly built QSEIS2025 library for new calculations; revalidate parameter choices and results.

```{eval-rst}
.. autofunction:: pygrnwang.read_qseis06.seek_qseis06
```

## read_qseis06_diff

**Deprecated backend: QSEIS06.** Its strain-rate and stress-rate readers are deprecated; use QSEIS2025 direct tensor outputs with a newly built and validated library for new calculations. The generic convert_strain2stress utility below is not deprecated.

This is an advanced historical workflow. It requires a dedicated perturbation
library with ASCII files retained. The docstrings record current binary-reader
and azimuth-rotation limitations; verify these conventions for your application.
The QSEIS2025 tutorial demonstrates directly computed strain and stress.

```{eval-rst}
.. autofunction:: pygrnwang.read_qseis06_diff.seek_qseis06_strain_rate_diff
```

```{eval-rst}
.. autofunction:: pygrnwang.read_qseis06_diff.seek_qseis06_stress_rate_diff
```

```{eval-rst}
.. autofunction:: pygrnwang.read_qseis06_diff.convert_strain2stress
```

## read_spgrn2012

**Deprecated backend: SPGRN2012.** Use the SPGRN2020 reader with a newly built SPGRN2020 library for new calculations; revalidate parameter choices, time alignment and results.

```{eval-rst}
.. autofunction:: pygrnwang.read_spgrn2012.seek_spgrn2012
```

## read_spgrn2020

```{eval-rst}
.. autofunction:: pygrnwang.read_spgrn2020.seek_spgrn2020
```

```{eval-rst}
.. autoclass:: pygrnwang.read_spgrn2020.GridGFCache
```

```{eval-rst}
.. automethod:: pygrnwang.read_spgrn2020.GridGFCache.get_block
```

```{eval-rst}
.. automethod:: pygrnwang.read_spgrn2020.GridGFCache.time_series
```

```{eval-rst}
.. automethod:: pygrnwang.read_spgrn2020.GridGFCache.tpts_table
```

```{eval-rst}
.. autofunction:: pygrnwang.read_spgrn2020.synthesize_from_cache
```

## read_qssp2020

```{eval-rst}
.. autofunction:: pygrnwang.read_qssp2020.seek_qssp2020
```

## read_by_qssp

```{eval-rst}
.. autofunction:: pygrnwang.read_by_qssp.read_by_qssp
```

## read_edcmp

```{eval-rst}
.. autofunction:: pygrnwang.read_edcmp.seek_edcmp2
```

```{eval-rst}
.. autofunction:: pygrnwang.read_edcmp.seek_edcmp2_bulk
```

## read_green_info_spgrn

```{eval-rst}
.. autofunction:: pygrnwang.read_green_info_spgrn.read_green_info_spgrn
```

## read_tpts_table

```{eval-rst}
.. autofunction:: pygrnwang.read_tpts_table.read_tpts_table
```

## utils

```{eval-rst}
.. autofunction:: pygrnwang.utils.read_tpts_table
```
