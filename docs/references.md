# References and citation

Cite the numerical method/backend actually used, and identify the
`pygrnwang` version or Git commit in the software and data description.
Record the Earth model, source convention, numerical parameters and
processing steps needed to reproduce the calculation. A package version
alone does not describe a Green's library.

## Backend methods

**QSEIS and stable Green's-function computation.**
Wang, R. (1999). A simple orthonormalization method for stable and efficient
computation of Green's functions. *Bulletin of the Seismological Society of
America*, **89**(3), 733–741.
[DOI: 10.1785/BSSA0890030733](https://doi.org/10.1785/BSSA0890030733).

**Static layered deformation with EDGRN/EDCMP.**
Wang, R. (2003). Computation of deformation induced by earthquakes in a
multi-layered elastic crust—FORTRAN programs EDGRN/EDCMP.
*Computers & Geosciences*, **29**(2), 195–207.
[DOI: 10.1016/S0098-3004(02)00111-5](https://doi.org/10.1016/S0098-3004(02)00111-5).

**Spherical/cylindrical harmonic calculation.**
Wang, R., & Wang, H. (2007). A fast converging and anti-aliasing algorithm
for Green's functions in terms of spherical or cylindrical harmonics.
*Geophysical Journal International*, **170**(1), 239–248.
[DOI: 10.1111/j.1365-246X.2007.03385.x](https://doi.org/10.1111/j.1365-246X.2007.03385.x).

**Spherical self-gravitating synthetic seismograms.**
Wang, R., Heimann, S., Zhang, Y., Wang, H., & Dahm, T. (2017).
Complete synthetic seismograms based on a spherical self-gravitating
earth model with an atmosphere–ocean–mantle–core structure.
*Geophysical Journal International*, **210**(3), 1739–1764.
[DOI: 10.1093/gji/ggx259](https://doi.org/10.1093/gji/ggx259).

**Dynamic Coulomb failure stress.**
Zhou, J., Wang, R., & Zhang, Y. (2026). DynCFS: a program for modeling
dynamic Coulomb failure stress changes in layered elastic media.
*Geophysical Journal International*, article ggaf534.
[DOI: 10.1093/gji/ggaf534](https://doi.org/10.1093/gji/ggaf534).

These references are the backend/method bibliography maintained by the
project. A paper's full formulation can be broader than the subset of
physics or observables enabled in one Python tutorial.

## Travel times

Crotwell, H. P., Owens, T. J., & Ritsema, J. (1999). The TauP Toolkit:
Flexible seismic travel-time and ray-path utilities.
*Seismological Research Letters*, **70**, 154–160.
This is the citation requested by the
[TauP project](https://www.seis.sc.edu/TauP/).
Report the version used: this package bundles TauP 2.6.1 for its Java
bridge and uses the installed ObsPy version for its alternative backend.

If results use ObsPy processing or its travel-time backend, follow the
[ObsPy citation guidance](https://docs.obspy.org/citations.html) for the
components and version involved.

## Software and data

The software source is
[Zhou-Jiangcheng/pygrnwang](https://github.com/Zhou-Jiangcheng/pygrnwang).
Record the installed version and a commit for an unreleased checkout.
When archiving results, include the generated input files, complete model,
`green_lib_info.json` and a small executable driver.

The tutorials use bundled AK135 elastic velocities/density with deliberately
chosen constant `Qp=600` and `Qs=300`. Cite the actual model source and
describe these attenuation choices when reusing the tutorial model in a
study; do not label it as an unmodified AK135-F attenuation model.
The separate [15 September comparison](guides/backend-comparison.md) instead
uses the published model file listed in its downloads, with depth-dependent
AK135-FC attenuation. Its source setup and component overlays refer to
DynCFS Fig. 3, extended from the paper's 10 km distance to 300/600/900 km.
The exact model hash and processing choices distinguish this result from
the lightweight tutorial examples.

For a publication, replace demonstration settings with a documented
model and convergence study suited to the intended observations.

## Source audit used for this documentation

The scientific conventions were checked against the Python moment conversion,
reader synthesis/rotation and material normalization, together with:

- QSEIS2025 `qswvint.f` and `qsgetinp.f` for strain/stress kernels and file names.
- QSEIS input templates for axes, reduction velocity and wavelet duration.
- QSSP2020 `qpfftinv.f` for output selections and the input source/receiver conventions.
- SPGRN2020 `qpdocu.f` for native arrival table units.
- EDGRN/EDCMP input writers for distance conversion, basis source area/slip and output layout.

The audit records existing inconsistencies in
[known limitations](guides/troubleshooting.md#known-implementation-limitations).
The documentation changes do not claim to replace numerical method
validation or alter solver behavior.
