# **Introduction**

This Python package serves as the frontend for calculating and building a Green's function library for synthetic seismograms. The backend consists of Wang Rongjiang's program for calculating synthetic seismograms, including EDGRN/EDCMP, [QSEIS_STRESS](https://github.com/Zhou-Jiangcheng/QSEIS_2006_STRESS), SPGRN, and QSSP (Wang, 1999; Wang 2003; Wang and Wang 2007; Wang et al., 2017). The code includes two parallel modes: one using the multiprocessing library (single-node multi-process) and the other using MPI (multi-node).

# Installation

1. For user mode (with Python 3.12)

```
pip install pygrnwang
```

2. For developer mode

```
conda create -n pygrnwang python=3.12
conda install gfortran obspy numpy scipy pandas matplotlib tqdm mpi4py -c conda-forge
git clone https://github.com/Zhou-Jiangcheng/pygrnwang.git
cd pygrnwang
pip install -e . --no-build-isolation
```

# Usage
1. An example for creating a Green's function library with qssp2020

```
from pygrnwang.create_qssp2020_bulk import *

if __name__ == '__main__':
    wavelet_duration = 0
    sampling_interval = 1
    time_window = 4096 - sampling_interval
    path_green = r'path\grns_qssp2020\ak135fc'
    os.makedirs(path_green, exist_ok=True)
    output_observables = [0 for _ in range(11)]
    output_observables[0] = 1
    output_observables[1] = 1
    output_observables[2] = 1
    pre_process_qssp2020(
        processes_num=24,
        path_green=path_green,
        event_depth_list=[h for h in range(1, 41, 2)],
        receiver_depth_list=[0],
        dist_range=[3000, 12000],
        delta_dist=10,
        spec_time_window=time_window,
        sampling_interval=sampling_interval,
        max_frequency=0.2,
        max_slowness=0.4,
        anti_alias=0.01,
        turning_point_filter=0,
        turning_point_d1=0,
        turning_point_d2=0,
        free_surface_filter=1,
        gravity_fc=0,
        gravity_harmonic=0,
        cal_sph=1,
        cal_tor=1,
        min_harmonic=4000,
        max_harmonic=10000,
        source_radius=0,
        source_duration=wavelet_duration * sampling_interval,
        output_observables=output_observables,
        time_window=time_window,
        time_reduction=-20,
        path_nd=r'path\ak135fc.nd',
        earth_model_layer_num=None,
        physical_dispersion=0,
        check_finished_tpts_table=False
    )
    create_grnlib_qssp2020_parallel(
        path_green=path_green, check_finished=False, cal_spec=False
    )

```

2. An example for reading from a Green's function library created by qssp2020
```
from pygrnwang.read_qssp2020 import seek_qssp2020


if __name__ == "__main__":
    seismograms, tpts_table, first_p, first_s, grn_dep, grn_receiver, green_dist = (
        seek_qssp2020(
            path_green="/e/grns_test/test_qssp",
            event_depth_km=10,
            receiver_depth_km=0,
            az_deg=60,
            dist_km=5000,
            focal_mechanism=[30, 40, 50],
            srate=1,
            before_p=20,
            pad_zeros=False,
            shift=False,
            rotate=True,
            only_seismograms=False,
            output_type='disp',
            model_name=r"path\ak135fc.nd",
        )
    )

    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(nrows=3, ncols=1)
    axs[0].plot(seismograms[0])
    axs[1].plot(seismograms[1])
    axs[2].plot(seismograms[2])
    plt.show()

```

# Reference

Wang, R. (1999). A simple orthonormalization method for stable and efficient computation of Green’s functions.  *Bulletin of the Seismological Society of America* ,  *89* (3), 733–741. https://doi.org/10.1785/BSSA0890030733

Wang, R. (2003). Computation of deformation induced by earthquakes in a multi-layered elastic crust—FORTRAN programs EDGRN/EDCMP. Computers & Geosciences, 29(2), 195–207. https://doi.org/10.1016/S0098-3004(02)00111-5

Wang, R., & Wang, H. (2007). A fast converging and anti-aliasing algorithm for green’s functions in terms of spherical or cylindrical harmonics. Geophysical Journal International, 170(1), 239–248. https://doi.org/10.1111/j.1365-246X.2007.03385.x

Wang, R., Heimann, S., Zhang, Y., Wang, H., & Dahm, T. (2017). Complete synthetic seismograms based on a spherical self-gravitating earth model with an atmosphere–ocean–mantle–core structure. Geophysical Journal International, 210(3), 1739–1764. https://doi.org/10.1093/gji/ggx259

Zhou, J., Wang, R., & Zhang, Y. (2026). DynCFS: a program for modeling dynamic coulomb failure stress changes in layered elastic media. Geophysical Journal International, ggaf534. https://doi.org/10.1093/gji/ggaf534