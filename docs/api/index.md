# API reference

This reference defines the documented support surface for pygrnwang 3.0. Import objects from their modules; the package root does not re-export the whole API.

**Deprecated backends:** QSEIS06 and SPGRN2012. Use QSEIS2025 and SPGRN2020, respectively, for new calculations. Deprecated entry points remain documented for existing workflows. Migration requires rebuilding Green libraries and revalidating numerical settings, source and time conventions, and results; replacement backends are not guaranteed to accept existing libraries or reproduce results without validation.

```python
from pygrnwang.create_qseis2025_bulk import pre_process_qseis2025
from pygrnwang.read_qseis2025 import seek_qseis2025
```

The explicit [public API manifest](public-api.json) lists 92 functions, classes and methods. It is reviewed when interfaces are added; documentation builds verify that every entry exists and every argument is documented. Names without a leading underscore are not automatically part of this surface.

```{eval-rst}
.. autosummary::
   :nosignatures:

   pygrnwang.create_qseis2025_bulk.pre_process_qseis2025
   pygrnwang.read_syn.read_syn
   pygrnwang.read_qseis2025.seek_qseis2025
   pygrnwang.read_spgrn2020.GridGFCache
   pygrnwang.read_edcmp.seek_edcmp2
   pygrnwang.pytaup.cal_first_p_s
   pygrnwang.focal_mechanism.check_convert_fm
```

```{toctree}
:maxdepth: 2

builders
readers
traveltimes
science
advanced
```

Numerical choices, units and complete runnable workflows are explained in the user guides and backend tutorials. Reference pages describe actual behavior, including historical parameter defaults that require an explicit model path in practice.
