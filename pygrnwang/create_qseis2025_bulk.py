import os
import math
import pickle
import json
import datetime

from tqdm import tqdm

from .create_qseis2025 import (
    QSEIS2025_NRMAX,
    QSEIS2025_MEMORY_PER_JOB_GB,
    create_dir_qseis2025,
    create_inp_qseis2025,
    call_qseis2025,
    check_output_qseis2025,
    convert_pd2bin_qseis2025,
    read_output_observables_qseis2025,
)
from .pytaup import create_tpts_table
from .utils import (
    group,
    convert_earth_model_nd2nd_without_Q,
    cal_grid,
    check_path_lengths,
    run_checked_job,
    run_jobs_sequential,
    run_jobs_parallel,
    run_jobs_mpi,
    run_until_complete,
    finish_library,
)


def _load_green_info(path_green):
    with open(
        os.path.join(path_green, "green_lib_info.json"), "r", encoding="utf-8"
    ) as fr:
        return json.load(fr)


def _job_dir(path_green, event_depth, receiver_depth, n_group):
    return str(
        os.path.join(
            path_green, "%.2f" % event_depth, "%.2f" % receiver_depth, "%d_0" % n_group
        )
    )


def _check_job(path_green, green_info, event_depth, receiver_depth, n_group,
               check_values=False):
    sub_sub_dir = _job_dir(path_green, event_depth, receiver_depth, n_group)
    path_inp = os.path.join(sub_sub_dir, "grn.inp")
    if not os.path.exists(path_inp):
        return ["%s is missing" % path_inp]
    try:
        output_observables = read_output_observables_qseis2025(path_inp)
    except (OSError, ValueError):
        return ["%s cannot be read" % path_inp]
    N_each_group = green_info["N_each_group"]
    return check_output_qseis2025(
        sub_sub_dir,
        output_observables,
        green_info["sampling_num"],
        min(N_each_group, green_info["N_dist"] - n_group * N_each_group),
        check_values,
    )


def _run_job(task):
    """Run one job unless check_finished finds it complete; return its problems."""
    event_depth, receiver_depth, n_group, path_green, check_finished = task
    green_info = _load_green_info(path_green)
    return run_checked_job(
        _job_dir(path_green, event_depth, receiver_depth, n_group),
        check_finished,
        lambda: call_qseis2025(event_depth, receiver_depth, n_group, path_green),
        lambda: _check_job(path_green, green_info, event_depth, receiver_depth, n_group),
        # outputs of an earlier run would pass the check if this run fails
        stale=["grn_*.bin*", "ex.*", "ss.*", "ds.*", "cl.*"],
    )


def _tasks(path_green, check_finished):
    with open(os.path.join(path_green, "group_list.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    return [[tuple(job + [path_green, check_finished]) for job in grp] for grp in group_list]


def _finish(path_green, run_problems, convert_pd2bin, remove_pd):
    """Convert the jobs, check the whole library and raise if incomplete."""
    if convert_pd2bin:
        convert_pd2bin_qseis2025_all(path_green, remove_pd)
    finish_library(path_green, run_problems, lambda: check_grnlib_qseis2025(path_green))


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_qseis2025(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    dist_range,
    delta_dist,
    N_each_group,
    time_window,
    sampling_interval,
    output_observables,
    slowness_int_algorithm=0,
    eps_estimate_wavenumber=1e-6,
    source_radius_ratio=0.05,
    slowness_window=None,
    time_reduction_velo=0,
    wavenumber_sampling_rate=12,
    anti_alias=0.01,
    free_surface=0,
    wavelet_duration=0,
    wavelet_type=1,
    flat_earth_transform=True,
    path_nd=None,
    earth_model_layer_num=None,
    check_finished_tpts_table=False,
):
    """Prepare the qseis2025 library grid, input files and job groups.

    Parameters
    ----------
    processes_num : int
        Positive worker count used to group jobs; MPI rank count must match the prepared group width.
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    event_depth_list : list of float
        Source depth nodes in km, positive down; supply a nonempty sorted list.
    receiver_depth_list : list of float
        Receiver depth nodes in km, positive down; supply a nonempty sorted list.
    dist_range : list of float
        Minimum and maximum epicentral distances in km.
    delta_dist : float
        Positive regular distance increment in km. The last grid point can exceed the requested maximum by less than one increment.
    N_each_group : int
        Positive maximum number of distances in each backend input file; at most 101, the distance limit (nrmax) of the bundled qseis2025 build.
    time_window : float
        Output time-window duration in seconds.
    sampling_interval : float
        Time step in seconds; choose it consistently with the highest modeled frequency.
    output_observables : list of int
        Five 0/1 flags in displacement (or velocity), volume, strain, stress, rotation order; rate versus integrated quantities depend on wavelet_type.
    slowness_int_algorithm : int, optional
        QSEIS integration selector: 0 for the full wavefield; 1 or 2 for narrow tapered slowness windows. Default: 0.
    eps_estimate_wavenumber : float, optional
        Full-wavefield wavenumber truncation tolerance; smaller values increase accuracy and computation. Default: 1e-06.
    source_radius_ratio : float, optional
        Dimensionless Gaussian spatial-smoothing ratio. At each frequency and
        receiver, the native radius is ``source_radius_ratio *
        min(sqrt(r**2 + (zs-zr)**2), Vp_source/(f+df))``, using epicentral
        distance r, source/receiver depths zs/zr, source-layer P-wave speed,
        frequency f and FFT increment df in consistent units. Coordinates and
        speed include any selected Earth flattening. The kernel is multiplied
        by ``exp(-(k*radius)**2/2)`` at wavenumber k; this is not a fixed-radius
        source disk. Zero disables smoothing. Larger values generally reduce
        the automatically estimated wavenumber cutoff and computation time.
        Default: 0.05.
    slowness_window : list of float or None, optional
        Four ordered slowness taper corners in s/km; None writes zeros for backend automatic limits. Default: None.
    time_reduction_velo : float, optional
        Reduction velocity in km/s; nonzero starts each trace at distance/velocity seconds, while zero disables reduction. Default: 0.
    wavenumber_sampling_rate : float, optional
        Dimensionless spatial Nyquist oversampling factor for wavenumber integration. Default: 12.
    anti_alias : float, optional
        Dimensionless time-domain alias suppression factor; use a small positive value below 1. Default: 0.01.
    free_surface : bool or int, optional
        Backend free-surface selection; see Notes for the backend-specific encoding. Default: 0.
    wavelet_duration : int, optional
        Wavelet duration in native time samples, not seconds. Nonpositive values request the backend default of two samples. Default: 0.
    wavelet_type : int, optional
        1 selects the normalized squared half-sinusoid; 2 selects its tapered Heaviside integral. Custom type 0 requires manually supplying wavelet samples; readers treat them as a moment-rate function, like type 1. Default: 1.
    flat_earth_transform : bool, optional
        Apply the backend flat-Earth transformation and receiver-radius distance correction. Default: True.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.
    check_finished_tpts_table : bool, optional
        Reuse existing P/S table files without validating their model or grid provenance. Default: False.

    Returns
    -------
    group_list : list
        Jobs grouped by processes_num; the same groups are saved as a pickle file.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    ValueError
        N_each_group exceeds the distance limit of the qseis2025 executable, or a job input path is longer than the 160 characters it can read.

    Notes
    -----
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. free_surface=0 includes the free surface; 1 removes it; 2 removes it with amplitude correction for surface receivers.
    """
    if N_each_group > QSEIS2025_NRMAX:
        # qseis2025 would stop with "nr > nrmax" and write no Green's functions
        raise ValueError(
            "N_each_group=%d exceeds %d, the number of distances qseis2025 can "
            "compute in one job" % (N_each_group, QSEIS2025_NRMAX)
        )
    N_dist = len(cal_grid(dist_range[0], dist_range[1], delta_dist))
    n_group_last = math.ceil(N_dist / N_each_group) - 1
    check_path_lengths(
        os.path.join(_job_dir(path_green, event_dep, receiver_dep, n_group_last), "grn.inp")
        for event_dep in event_depth_list
        for receiver_dep in receiver_depth_list
    )
    print("Preprocessing")
    os.makedirs(path_green, exist_ok=True)

    N_dist, N_dist_group = None, None
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            # print("creating green func dir, info, inp for event_depth=%.2f receiver_depth=%.2f" %
            #       (event_depth, receiver_depth))
            N_dist, N_dist_group = create_dir_qseis2025(
                path_green=path_green,
                event_depth=event_depth,
                receiver_depth=receiver_depth,
                dist_range=dist_range,
                delta_dist=delta_dist,
                N_each_group=N_each_group,
            )
            create_inp_qseis2025(
                path_green=path_green,
                event_depth=event_depth,
                receiver_depth=receiver_depth,
                dist_range=dist_range,
                delta_dist=delta_dist,
                N_dist=N_dist,
                N_dist_group=N_dist_group,
                N_each_group=N_each_group,
                time_window=time_window,
                sampling_interval=sampling_interval,
                output_observables=output_observables,
                slowness_int_algorithm=slowness_int_algorithm,
                eps_estimate_wavenumber=eps_estimate_wavenumber,
                source_radius_ratio=source_radius_ratio,
                slowness_window=slowness_window,
                time_reduction_velo=time_reduction_velo,
                wavenumber_sampling_rate=wavenumber_sampling_rate,
                anti_alias=anti_alias,
                free_surface=free_surface,
                wavelet_duration=wavelet_duration,
                wavelet_type=wavelet_type,
                flat_earth_transform=flat_earth_transform,
                path_nd=path_nd,
                earth_model_layer_num=earth_model_layer_num,
            )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    # creating tp and ts tables
    dist_kms = cal_grid(dist_range[0], dist_range[1], delta_dist)
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                path_green,
                event_depth,
                receiver_depth,
                dist_kms,
                path_nd_without_Q,
                check_finished_tpts_table,
            )

    green_info = {
        "processes_num": processes_num,
        "event_depth_list": event_depth_list,
        "receiver_depth_list": receiver_depth_list,
        "grn_dist_range": dist_range,
        "grn_delta_dist": delta_dist,
        "N_dist": N_dist,
        "N_dist_group": N_dist_group,
        "N_each_group": N_each_group,
        "time_window": time_window,
        "sampling_interval": sampling_interval,
        "sampling_num": round(time_window / sampling_interval + 1),
        "slowness_int_algorithm": slowness_int_algorithm,
        "slowness_window": slowness_window,
        "time_reduction_velo": time_reduction_velo,
        "wavenumber_sampling_rate": wavenumber_sampling_rate,
        "anti_alias": anti_alias,
        "free_surface": free_surface,
        "wavelet_duration": wavelet_duration,
        "wavelet_type": wavelet_type,
        "flat_earth_transform": flat_earth_transform,
        "path_nd": path_nd,
        "path_nd_without_Q": path_nd_without_Q,
        "earth_model_layer_num": earth_model_layer_num,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    inp_list = []
    for event_dep in event_depth_list:
        for receiver_dep in receiver_depth_list:
            for nn in range(N_dist_group):
                inp_list.append([event_dep, receiver_dep, nn])
    # inp_list_sorted = sorted(inp_list, key=lambda x: abs(x[0] - x[1]))
    group_list = group(inp_list, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore
    return group_list


def create_grnlib_qseis2025_sequential(
    path_green, check_finished=False, convert_pd2bin=True, remove_pd=True, max_retries=2
):
    """Compute the prepared qseis2025 library sequentially.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    convert_pd2bin : bool, optional
        Convert completed ASCII waveforms to the compact float32 reader format. Default: True.
    remove_pd : bool, optional
        Delete original ASCII output; retain it while validating a new calculation. Default: True.
    max_retries : int, optional
        Extra passes over the jobs that did not complete, for example because they ran out of memory; each pass recomputes only those jobs, with the same worker count. Default: 2.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    RuntimeError
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_qseis2025).

    Notes
    -----
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    problems = run_until_complete(
        lambda tasks: run_jobs_sequential(_run_job, tasks, desc="Computing dynamic stress"),
        sum(_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems, convert_pd2bin, remove_pd)


def create_grnlib_qseis2025_parallel(
    path_green,
    check_finished=False,
    convert_pd2bin=True,
    remove_pd=True,
    memory_per_job_gb=QSEIS2025_MEMORY_PER_JOB_GB,
    max_retries=2,
):
    """Compute the prepared qseis2025 library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    convert_pd2bin : bool, optional
        Convert completed ASCII waveforms to the compact float32 reader format. Default: True.
    remove_pd : bool, optional
        Delete original ASCII output; retain it while validating a new calculation. Default: True.
    memory_per_job_gb : float or None, optional
        Peak memory of one qseis2025 process in GiB. A RuntimeWarning is issued when min(processes_num, number of jobs) such processes may not fit in the currently available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again; None skips the warning. Default: 1.4, the measured commit of the bundled executable.
    max_retries : int, optional
        Extra passes over the jobs that did not complete, for example because they ran out of memory; each pass recomputes only those jobs, with the same worker count. Default: 2.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    RuntimeError
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_qseis2025).

    Notes
    -----
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A job fails when the executable cannot start, exits with an error (for example after being killed for lack of memory) or leaves incomplete output files; it never gets a .finished marker. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs. On Windows call under an if __name__ == "__main__" guard.
    """
    processes = _load_green_info(path_green).get("processes_num", None)
    problems = run_until_complete(
        lambda tasks: run_jobs_parallel(
            _run_job, tasks, processes, memory_per_job_gb, desc="Compute QSEIS2025 Green's library"
        ),
        sum(_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems, convert_pd2bin, remove_pd)


def convert_pd2bin_qseis2025_all(path_green, remove=False):
    """Convert all completed qseis2025 outputs to float32 binary libraries.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    remove : bool, optional
        Delete source ASCII files after converting them; keep False when inspecting backend output. Default: False.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.

    Notes
    -----
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units. Jobs with incomplete output are skipped and keep their files, so they can be inspected and recomputed.
    """
    print("Converting ascii files to bytes files")
    green_info = _load_green_info(path_green)
    for event_dep in green_info["event_depth_list"]:
        for receiver_dep in green_info["receiver_depth_list"]:
            for n_group in range(green_info["N_dist_group"]):
                if _check_job(path_green, green_info, event_dep, receiver_dep, n_group):
                    continue
                convert_pd2bin_qseis2025(
                    _job_dir(path_green, event_dep, receiver_dep, n_group),
                    remove=remove,
                )


def check_grnlib_qseis2025(path_green, check_values=False):
    """Check that a qseis2025 library holds every file the readers need.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_values : bool, optional
        Also read every binary Green's function file and report NaN or infinite values; this reads the whole library. Default: False.

    Returns
    -------
    problems : list of str
        One line per missing or incomplete file; empty when the library is complete.

    Raises
    ------
    OSError
        green_lib_info.json cannot be read.

    Notes
    -----
    For every source depth, receiver depth and distance group, each observable selected in the job's grn.inp must exist as a float32 binary file of the size the readers expect, or as complete ASCII files with one row per sample and one column per distance. The P and S travel-time tables must hold one value per distance. Rerun the create_grnlib function with check_finished=True to recompute only incomplete jobs.
    """
    green_info = _load_green_info(path_green)
    problems = []
    for event_dep in green_info["event_depth_list"]:
        for receiver_dep in green_info["receiver_depth_list"]:
            sub_dir = os.path.join(path_green, "%.2f" % event_dep, "%.2f" % receiver_dep)
            for name in ["tp_table.bin", "ts_table.bin"]:
                path_table = os.path.join(sub_dir, name)
                if not os.path.exists(path_table):
                    problems.append("%s is missing" % path_table)
                elif os.path.getsize(path_table) != 4 * green_info["N_dist"]:
                    problems.append(
                        "%s has %d bytes, expected %d"
                        % (path_table, os.path.getsize(path_table), 4 * green_info["N_dist"])
                    )
            for n_group in range(green_info["N_dist_group"]):
                problems += _check_job(
                    path_green, green_info, event_dep, receiver_dep, n_group, check_values
                )
    return problems


def create_grnlib_qseis2025_parallel_multi_nodes(
    path_green,
    check_finished=False,
    convert_pd2bin=True,
    remove_pd=True,
    memory_per_job_gb=QSEIS2025_MEMORY_PER_JOB_GB,
    max_retries=2,
):
    """Compute the prepared qseis2025 library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    convert_pd2bin : bool, optional
        Convert completed ASCII waveforms to the compact float32 reader format. Default: True.
    remove_pd : bool, optional
        Delete original ASCII output; retain it while validating a new calculation. Default: True.
    memory_per_job_gb : float or None, optional
        Peak memory of one qseis2025 process in GiB. A RuntimeWarning is issued when the ranks on a node may not fit in its available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again; None skips the warning. Default: 1.4, the measured commit of the bundled executable.
    max_retries : int, optional
        Extra passes over the jobs that did not complete, for example because they ran out of memory; each pass recomputes only those jobs, with the same worker count. Default: 2.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    RuntimeError
        mpi4py is unavailable, or (on rank 0) jobs still fail after the retries or the finished library is incomplete.
    ValueError
        MPI rank count does not match the prepared group width.

    Notes
    -----
    See the qseis2025 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A failed job does not stop the others: after all ranks finish, the jobs that did not complete are shared among the ranks and computed again, up to max_retries times. Rank 0 checks the library after all ranks finish; after an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    rank, problems = run_jobs_mpi(
        MPI, _tasks(path_green, check_finished), _run_job, memory_per_job_gb,
        max_retries=max_retries,
    )
    if rank == 0:
        _finish(path_green, problems, convert_pd2bin, remove_pd)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


if __name__ == "__main__":
    pass
