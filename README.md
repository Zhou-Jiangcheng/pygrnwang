# **Introduction**

This Python package serves as the frontend for calculating and building a Green's function library for synthetic seismograms. The backend consists of Wang Rongjiang's program for calculating synthetic seismograms, including QSEIS, SPGRN, and QSSP (Wang, 1999; Wang et al., 2017). The code includes two parallel modes: one using the multiprocessing library (single-node multi-process) and the other using MPI (multi-node).

Wang, R. (1999). A simple orthonormalization method for stable and efficient computation of Green’s functions.  *Bulletin of the Seismological Society of America* ,  *89* (3), 733–741. [https://doi.org/10.1785/BSSA0890030733](https://doi.org/10.1785/BSSA0890030733)

Wang, R., Heimann, S., Zhang, Y., Wang, H., & Dahm, T. (2017). Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmosphere–ocean–mantle–core structure.  *Geophysical Journal International* ,  *210* (3), 1739–1764.

# Installation

1. Install the requirments. (Debian 12, Python 3.11)

```
sudo apt install gfortran
conda install numpy scipy pandas mpi4py -c conda-forge
```

2. Add the folder path containing this Python library to the environment.

For example, if the folder path is `/home/pygrnwang` and you have already created a Python 3.11 virtual environment named `pygrnwang` using Anaconda, add `/home` to the `/path_anaconda/envs/pygrnwang/lib/python3.11/site-packages/custom.pth` file.

3. Compile the corresponding Fortran source files (in `qseis06_src`, `qssp2020_src`, and `spgrn2020_src`, respectively).

```
cd ./qseis06_src
gfortran ./*.f -O3 -o qseis06.bin
cp qseis06.bin ../
cd ../spgrn2020_src
gfortran ./*.f -O3 -o spgrn2020.bin
cp spgrn2020.bin ../
cd ../qssp2020_src
gfortran ./*.f -O3 -o qssp2020.bin
cp qssp2020.bin ../
```
