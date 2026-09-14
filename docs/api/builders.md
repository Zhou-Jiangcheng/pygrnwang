# Library preparation and execution

Prepare inputs first, run a matching executor, then read the outputs. All paths passed to executors should be absolute. Use the backend tutorials to select numerical parameters and to understand model and spectrum reuse.

The native programs can change the process working directory. The Python execution functions report progress; inspect result files and logs to establish calculation success. MPI entries require `mpi4py` and an MPI launch with the prepared group width.

## create_edgrn_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_edgrn_bulk.pre_process_edgrn2
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edgrn_bulk.create_grnlib_edgrn2_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edgrn_bulk.create_grnlib_edgrn2_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edgrn_bulk.create_grnlib_edgrn2_parallel_multi_nodes
```

## create_edcmp_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_edcmp_bulk.pre_process_edcmp2
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edcmp_bulk.create_grnlib_edcmp2_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edcmp_bulk.create_grnlib_edcmp2_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edcmp_bulk.create_grnlib_edcmp2_parallel_multi_nodes
```

```{eval-rst}
.. autofunction:: pygrnwang.create_edcmp_bulk.convert_pd2bin_edcmp2_all
```

## create_qseis06_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis06_bulk.pre_process_qseis06
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis06_bulk.pre_process_qseis06_strain_rate
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis06_bulk.create_grnlib_qseis06_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis06_bulk.create_grnlib_qseis06_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis06_bulk.create_grnlib_qseis06_parallel_multi_nodes
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis06_bulk.convert_pd2bin_qseis06_all
```

## create_qseis2025_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis2025_bulk.pre_process_qseis2025
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis2025_bulk.create_grnlib_qseis2025_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis2025_bulk.create_grnlib_qseis2025_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis2025_bulk.create_grnlib_qseis2025_parallel_multi_nodes
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qseis2025_bulk.convert_pd2bin_qseis2025_all
```

## create_spgrn2012_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2012_bulk.pre_process_spgrn2012
```

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2012_bulk.create_grnlib_spgrn2012_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2012_bulk.create_grnlib_spgrn2012_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2012_bulk.create_grnlib_spgrn2012_parallel_multi_nodes
```

## create_spgrn2020_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2020_bulk.pre_process_spgrn2020
```

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2020_bulk.create_grnlib_spgrn2020_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2020_bulk.create_grnlib_spgrn2020_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_spgrn2020_bulk.create_grnlib_spgrn2020_parallel_multi_nodes
```

## create_qssp2020_bulk

```{eval-rst}
.. autofunction:: pygrnwang.create_qssp2020_bulk.pre_process_qssp2020
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qssp2020_bulk.create_grnlib_qssp2020_sequential
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qssp2020_bulk.create_grnlib_qssp2020_parallel
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qssp2020_bulk.create_grnlib_qssp2020_spec_parallel_multi_nodes
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qssp2020_bulk.create_grnlib_qssp2020_func_parallel_multi_nodes
```

```{eval-rst}
.. autofunction:: pygrnwang.create_qssp2020_bulk.convert_pd2bin_qssp2020_all
```
