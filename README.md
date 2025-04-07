# **Introduction**

This Python package serves as the frontend for calculating and building a Green's function library for synthetic seismograms. The backend consists of Wang Rongjiang's program for calculating synthetic seismograms, including QSEIS, SPGRN, and QSSP (Wang, 1999; Wang et al., 2017). The code includes two parallel modes: one using the multiprocessing library (single-node multi-process) and the other using MPI (multi-node).

Wang, R. (1999). A simple orthonormalization method for stable and efficient computation of Green’s functions.  *Bulletin of the Seismological Society of America* ,  *89* (3), 733–741. [https://doi.org/10.1785/BSSA0890030733](https://doi.org/10.1785/BSSA0890030733)

Wang, R., Heimann, S., Zhang, Y., Wang, H., & Dahm, T. (2017). Complete synthetic seismograms based on a spherical self-gravitating Earth model with an atmosphere–ocean–mantle–core structure.  *Geophysical Journal International* ,  *210* (3), 1739–1764.

# Installation

1. Install the requirments. (Debian 12, Python 3.11)

```
sudo apt install default-jdk
conda create -n pygrnwang python=3.9
conda install gfortran numpy scipy pandas mpi4py -c conda-forge
git clone https://github.com/Zhou-Jiangcheng/pygrnwang.git
cd pygrnwang
pip install .
```
