# **Introduction**

This Python package serves as the frontend for calculating and building a Green's function library for synthetic seismograms. The backend consists of Wang Rongjiang's program for calculating synthetic seismograms, including EDGRN/EDCMP, QSEIS, SPGRN, and QSSP (Wang, 1999; Wang 2003; Wang and Wang 2007; Wang et al., 2017). The code includes two parallel modes: one using the multiprocessing library (single-node multi-process) and the other using MPI (multi-node).

Wang, R. (1999). A simple orthonormalization method for stable and efficient computation of Green’s functions.  *Bulletin of the Seismological Society of America* ,  *89* (3), 733–741. [https://doi.org/10.1785/BSSA0890030733](https://doi.org/10.1785/BSSA0890030733)

Wang, R. (2003). Computation of deformation induced by earthquakes in a multi-layered elastic crust—FORTRAN programs EDGRN/EDCMP. Computers & Geosciences, 29(2), 195–207. https://doi.org/10.1016/S0098-3004(02)00111-5

Wang, R., & Wang, H. (2007). A fast converging and anti-aliasing algorithm for green’s functions in terms of spherical or cylindrical harmonics. Geophysical Journal International, 170(1), 239–248. https://doi.org/10.1111/j.1365-246X.2007.03385.x

Wang, R., Heimann, S., Zhang, Y., Wang, H., & Dahm, T. (2017). Complete synthetic seismograms based on a spherical self-gravitating earth model with an atmosphere–ocean–mantle–core structure. Geophysical Journal International, 210(3), 1739–1764. https://doi.org/10.1093/gji/ggx259

# Installation

1. Install the requirments.
```
conda create -n pygrnwang python=3.11
conda install openjdk jpype1 gfortran numpy scipy pandas tqdm mpi4py -c conda-forge
```
2. Install by pip.
```
git clone https://github.com/Zhou-Jiangcheng/pygrnwang.git
cd pygrnwang
pip install .
```
