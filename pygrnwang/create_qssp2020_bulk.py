import os
import glob
import pickle
import json
import datetime

import numpy as np
from tqdm import tqdm

from .create_qssp2020 import (
    mt_com_list,
    create_inp_qssp2020,
    create_dir_qssp2020,
    call_qssp2020,
    check_spec_qssp2020,
    check_func_qssp2020,
    convert_pd2bin_qssp2020,
)
from .utils import (
    group,
    convert_earth_model_nd2nd_without_Q,
    cal_grid,
    check_path_lengths,
    check_file_size,
    run_checked_job,
    run_jobs_sequential,
    run_jobs_parallel,
    run_jobs_mpi,
    run_until_complete,
    finish_library,
    raise_if_incomplete,
)
from .pytaup import taup_create_npz_file, create_tpts_table


def _load_green_info(path_green):
    with open(
        os.path.join(path_green, "green_lib_info.json"), "r", encoding="utf-8"
    ) as fr:
        return json.load(fr)


def _job_dir(path_green, event_depth, receiver_depth, mt_com):
    if mt_com == "spec":
        return str(
            os.path.join(
                path_green, "GreenSpec", "%.2f" % event_depth, "%.2f" % receiver_depth
            )
        )
    return str(
        os.path.join(
            path_green,
            "GreenFunc",
            "%.2f" % event_depth,
            "%.2f" % receiver_depth,
            mt_com,
        )
    )


def _check_job(path_green, green_info, event_depth, receiver_depth, mt_com,
               check_values=False):
    if mt_com == "spec":
        return check_spec_qssp2020(path_green, event_depth, receiver_depth)
    dist_range = green_info["grn_dist_range"]
    return check_func_qssp2020(
        _job_dir(path_green, event_depth, receiver_depth, mt_com),
        green_info["output_observables"],
        green_info["sampling_num"],
        len(cal_grid(dist_range[0], dist_range[1], green_info["grn_delta_dist"])),
        check_values,
    )


def _spectra_problem(path_green, event_depth, receiver_depth, mt_com):
    """Return why a time-domain job cannot use its spectra, or None."""
    path_spec = _job_dir(path_green, event_depth, receiver_depth, "spec")
    if check_spec_qssp2020(path_green, event_depth, receiver_depth) or os.path.exists(
        os.path.join(path_spec, ".failed")
    ):
        return "%s failed: the spectra in %s are incomplete" % (
            _job_dir(path_green, event_depth, receiver_depth, mt_com),
            path_spec,
        )
    return None


def _run_job(task):
    """Run one job unless check_finished finds it complete; return its problems."""
    event_depth, receiver_depth, mt_com, path_green, check_finished = task
    if mt_com != "spec":
        problem = _spectra_problem(path_green, event_depth, receiver_depth, mt_com)
        if problem:
            return [problem]
    green_info = _load_green_info(path_green)
    return run_checked_job(
        _job_dir(path_green, event_depth, receiver_depth, mt_com),
        check_finished,
        lambda: call_qssp2020(event_depth, receiver_depth, mt_com, path_green),
        lambda: _check_job(path_green, green_info, event_depth, receiver_depth, mt_com),
        # outputs of an earlier run would pass the check if this run fails
        stale=["?_Green_*"] if mt_com == "spec" else ["_*.bin*", "_*.dat"],
    )


def _tasks(path_green, stage, check_finished):
    with open(os.path.join(path_green, "group_list_%s.pkl" % stage), "rb") as fr:
        group_list = pickle.load(fr)
    return [[tuple(job + [path_green, check_finished]) for job in grp] for grp in group_list]


def _run_stages(path_green, cal_spec, check_finished, max_retries, run_pass):
    """Run the spectral, then the time-domain jobs, recomputing failed ones.

    run_pass(tasks, desc) runs jobs and returns the failed ones. Time-domain
    jobs whose spectra are incomplete are not run. Returns the problems of the
    jobs still failing.
    """
    problems = []
    if cal_spec:
        problems += run_until_complete(
            lambda tasks: run_pass(tasks, "in the transformed domain"),
            sum(_tasks(path_green, "spec", check_finished), []),
            max_retries,
        )
    func_tasks = []
    for task in sum(_tasks(path_green, "func", check_finished), []):
        problem = _spectra_problem(path_green, *task[:3])
        if problem:
            problems.append(problem)
        else:
            func_tasks.append(task)
    problems += run_until_complete(
        lambda tasks: run_pass(tasks, "in the time domain"), func_tasks, max_retries
    )
    return problems


def _finish(path_green, run_problems, convert_pd2bin, remove_pd):
    """Convert the jobs, check the whole library and raise if incomplete."""
    if convert_pd2bin:
        convert_pd2bin_qssp2020_all(path_green)
    if remove_pd:
        remove_dat_files(path_green)
    finish_library(path_green, run_problems, lambda: check_grnlib_qssp2020(path_green))


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_spec(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    turning_point_filter,
    turning_point_d1,
    turning_point_d2,
    free_surface_filter,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    min_harmonic,
    max_harmonic,
    source_radius,
    source_duration,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    item_list_spec = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            create_dir_qssp2020(event_depth, receiver_depth, path_green)
            create_inp_qssp2020(
                "spec",
                path_green,
                event_depth,
                receiver_depth,
                spec_time_window,
                sampling_interval,
                max_frequency,
                max_slowness,
                anti_alias,
                turning_point_filter,
                turning_point_d1,
                turning_point_d2,
                free_surface_filter,
                gravity_fc,
                gravity_harmonic,
                cal_sph,
                cal_tor,
                min_harmonic,
                max_harmonic,
                source_radius,
                1,
                source_duration,
                [0 for _ in range(11)],
                time_window,
                time_reduction,
                dist_range,
                delta_dist,
                path_nd,
                earth_model_layer_num,
                physical_dispersion,
            )
            item_list_spec.append([event_depth, receiver_depth, "spec"])

    group_list_spec = group(item_list_spec, processes_num)
    with open(os.path.join(path_green, "group_list_spec.pkl"), "wb") as fw:
        pickle.dump(group_list_spec, fw)  # type: ignore
    return group_list_spec


def pre_process_func(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    turning_point_filter,
    turning_point_d1,
    turning_point_d2,
    free_surface_filter,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    min_harmonic,
    max_harmonic,
    source_radius,
    source_duration,
    output_observables: list,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    item_list_func = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            for mt_com in mt_com_list:
                create_inp_qssp2020(
                    mt_com,
                    path_green,
                    event_depth,
                    receiver_depth,
                    spec_time_window,
                    sampling_interval,
                    max_frequency,
                    max_slowness,
                    anti_alias,
                    turning_point_filter,
                    turning_point_d1,
                    turning_point_d2,
                    free_surface_filter,
                    gravity_fc,
                    gravity_harmonic,
                    cal_sph,
                    cal_tor,
                    min_harmonic,
                    max_harmonic,
                    source_radius,
                    0,
                    source_duration,
                    output_observables,
                    time_window,
                    time_reduction,
                    dist_range,
                    delta_dist,
                    path_nd,
                    earth_model_layer_num,
                    physical_dispersion,
                )
                item_list_func.append([event_depth, receiver_depth, mt_com])

    group_list_func = group(item_list_func, processes_num)
    with open(os.path.join(path_green, "group_list_func.pkl"), "wb") as fw:
        pickle.dump(group_list_func, fw)  # type: ignore
    return group_list_func


def pre_process_qssp2020(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    turning_point_filter,
    turning_point_d1,
    turning_point_d2,
    free_surface_filter,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    min_harmonic,
    max_harmonic,
    source_radius,
    source_duration,
    output_observables: list,
    time_window,
    time_reduction,
    dist_range,
    delta_dist,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
    check_finished_tpts_table=False,
):
    """Prepare the qssp2020 library grid, input files and job groups.

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
    spec_time_window : float
        Spectral calculation duration in seconds; it must cover the requested output window.
    sampling_interval : float
        Time step in seconds; choose it consistently with the highest modeled frequency.
    max_frequency : float
        Highest modeled frequency in Hz; must be compatible with the time step.
    max_slowness : float
        Maximum modeled slowness in s/km; in SPGRN, nonpositive requests the complete wavefield.
    anti_alias : float
        Dimensionless time-domain alias suppression factor; use a small positive value below 1.
    turning_point_filter : int
        1 selects the QSSP turning-depth filter; 0 disables it.
    turning_point_d1 : float
        Minimum allowed turning depth in km when the turning-point filter is enabled.
    turning_point_d2 : float
        Maximum allowed turning depth in km when the turning-point filter is enabled.
    free_surface_filter : int
        QSSP switch: 1 includes free-surface reflection, 0 removes it.
    gravity_fc : float
        Critical frequency in Hz below which self-gravity is included together with the harmonic cutoff.
    gravity_harmonic : int
        Critical spherical harmonic degree for self-gravity.
    cal_sph : int
        1 enables spheroidal (P-SV) modes; 0 disables them.
    cal_tor : int
        1 enables toroidal (SH) modes; 0 disables them.
    min_harmonic : int
        Control for estimating the low-frequency baseline of the
        frequency-dependent upper harmonic cutoff. This is not the lowest
        retained degree: the spectral sum still starts at degree zero.
        Converge it together with max_harmonic for the requested observables.
    max_harmonic : int
        Upper harmonic cutoff cap. It also influences the spatial
        differential-transform order through the internal maximum degree,
        so changing it can affect synthesis even when the stored spectral
        cutoff is unchanged. Converge it together with min_harmonic.
    source_radius : float
        Source patch radius in km.
    source_duration : float
        Squared half-sinusoid source-time-function duration in seconds; zero requests the backend minimum.
    output_observables : list of int
        Eleven 0/1 flags: displacement, velocity, acceleration, strain, strain rate, stress, stress rate, rotation, rotation rate, gravitation, gravimeter.
    time_window : float
        Output time-window duration in seconds.
    time_reduction : float
        QSSP trace start time in seconds relative to source origin, not a velocity.
    dist_range : list of float
        Minimum and maximum epicentral distances in km.
    delta_dist : float
        Positive regular distance increment in km. The last grid point can exceed the requested maximum by less than one increment.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.
    physical_dispersion : int, optional
        0 disables, 1 enables the backend physical-dispersion correction associated with attenuation. Default: 0.
    check_finished_tpts_table : bool, optional
        Reuse existing P/S table files without validating their model or grid provenance. Default: False.

    Returns
    -------
    None
        Writes backend inputs, metadata or output files to the library.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    ValueError
        An input, spectrum or output path is longer than the 160 characters qssp2020 can hold.

    Notes
    -----
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. A first computation must include both spectral and time-domain stages. See the :doc:`QSSP2020 tutorial </backends/qssp2020>` for harmonic-cutoff convergence and comparison settings.
    """
    paths = []
    for event_dep in event_depth_list:
        for receiver_dep in receiver_depth_list:
            path_spec = _job_dir(path_green, event_dep, receiver_dep, "spec")
            path_func = _job_dir(path_green, event_dep, receiver_dep, "mrr")
            paths += [
                os.path.join(path_spec, "spec.inp"),
                os.path.join(path_spec, "U_Green_%.2fkm" % event_dep),
                os.path.join(path_func, "mrr.inp"),
                # the longest output file name
                os.path.join(path_func, "_stress_rate_ee.dat"),
            ]
    check_path_lengths(paths)
    print("Preprocessing")
    pre_process_spec(
        processes_num,
        path_green,
        event_depth_list,
        receiver_depth_list,
        spec_time_window,
        sampling_interval,
        max_frequency,
        max_slowness,
        anti_alias,
        turning_point_filter,
        turning_point_d1,
        turning_point_d2,
        free_surface_filter,
        gravity_fc,
        gravity_harmonic,
        cal_sph,
        cal_tor,
        min_harmonic,
        max_harmonic,
        source_radius,
        source_duration,
        time_window,
        time_reduction,
        dist_range,
        delta_dist,
        path_nd,
        earth_model_layer_num,
        physical_dispersion,
    )
    pre_process_func(
        processes_num,
        path_green,
        event_depth_list,
        receiver_depth_list,
        spec_time_window,
        sampling_interval,
        max_frequency,
        max_slowness,
        anti_alias,
        turning_point_filter,
        turning_point_d1,
        turning_point_d2,
        free_surface_filter,
        gravity_fc,
        gravity_harmonic,
        cal_sph,
        cal_tor,
        min_harmonic,
        max_harmonic,
        source_radius,
        source_duration,
        output_observables,
        time_window,
        time_reduction,
        dist_range,
        delta_dist,
        path_nd,
        earth_model_layer_num,
        physical_dispersion,
    )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    # creating tp and ts tables
    npz_file = taup_create_npz_file(nd_file=path_nd_without_Q)
    dist_kms = cal_grid(dist_range[0], dist_range[1], delta_dist)
    for event_depth in tqdm(event_depth_list, desc="Creating travel time tables"):
        for receiver_depth in receiver_depth_list:
            create_tpts_table(
                os.path.join(path_green, "GreenFunc"),
                event_depth,
                receiver_depth,
                dist_kms,
                npz_file,
                check_finished_tpts_table,
            )

    green_info = {
        "processes_num": processes_num,
        "event_depth_list": event_depth_list,
        "receiver_depth_list": receiver_depth_list,
        "spec_time_window": spec_time_window,
        "sampling_interval": sampling_interval,
        "max_frequency": max_frequency,
        "max_slowness": max_slowness,
        "anti_alias": anti_alias,
        "turning_point_filter": turning_point_filter,
        "turning_point_d1": turning_point_d1,
        "turning_point_d2": turning_point_d2,
        "free_surface_filter": free_surface_filter,
        "gravity_fc": gravity_fc,
        "gravity_harmonic": gravity_harmonic,
        "cal_sph": cal_sph,
        "cal_tor": cal_tor,
        "min_harmonic": min_harmonic,
        "max_harmonic": max_harmonic,
        "source_radius": source_radius,
        "source_duration": source_duration,
        "output_observables": output_observables,
        "time_window": time_window,
        "sampling_num": round(time_window / sampling_interval) + 1,
        "time_reduction": time_reduction,
        "grn_dist_range": dist_range,
        "grn_delta_dist": delta_dist,
        "path_nd": path_nd,
        "path_nd_without_Q": path_nd_without_Q,
        "earth_model_layer_num": earth_model_layer_num,
        "physical_dispersion": physical_dispersion,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)


def create_grnlib_qssp2020_sequential(
    path_green, cal_spec=True, check_finished=False, convert_pd2bin=True, remove_pd=True, max_retries=2
):
    """Compute the prepared qssp2020 library sequentially.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    cal_spec : bool, optional
        Compute spectra before time-domain synthesis. Keep True for a new QSSP library; False requires compatible existing spectra. Default: True.
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
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_qssp2020).

    Notes
    -----
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    problems = _run_stages(
        path_green,
        cal_spec,
        check_finished,
        max_retries,
        lambda tasks, desc: run_jobs_sequential(
            _run_job, tasks, desc="Compute the Green's function library %s." % desc
        ),
    )
    _finish(path_green, problems, convert_pd2bin, remove_pd)


def create_grnlib_qssp2020_parallel(
    path_green,
    cal_spec=True,
    check_finished=False,
    convert_pd2bin=True,
    remove_pd=True,
    memory_per_job_gb=None,
    max_retries=2,
):
    """Compute the prepared qssp2020 library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    cal_spec : bool, optional
        Compute spectra before time-domain synthesis. Keep True for a new QSSP library; False requires compatible existing spectra. Default: True.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    convert_pd2bin : bool, optional
        Convert completed ASCII waveforms to the compact float32 reader format. Default: True.
    remove_pd : bool, optional
        Delete original ASCII output; retain it while validating a new calculation. Default: True.
    memory_per_job_gb : float or None, optional
        Peak memory of one qssp2020 process in GiB, for the larger of the spectral and time-domain stages. A RuntimeWarning is issued when min(processes_num, number of a stage's jobs) such processes may not fit in the currently available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again. qssp2020 allocates its arrays from the input, so there is no fixed value; None skips the warning. Default: None.
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
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_qssp2020).

    Notes
    -----
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A job fails when the executable cannot start, exits with an error (for example after being killed or refused memory) or leaves incomplete output files; it never gets a .finished marker. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs. On Windows call under an if __name__ == "__main__" guard.
    """
    processes = _load_green_info(path_green).get("processes_num", None)
    problems = _run_stages(
        path_green,
        cal_spec,
        check_finished,
        max_retries,
        lambda tasks, desc: run_jobs_parallel(
            _run_job,
            tasks,
            processes,
            memory_per_job_gb,
            desc="Compute QSSP2020 Green's library %s." % desc,
        ),
    )
    _finish(path_green, problems, convert_pd2bin, remove_pd)


def create_grnlib_qssp2020_spec_parallel_multi_nodes(
    path_green, check_finished=False, memory_per_job_gb=None, max_retries=2
):
    """Compute the prepared qssp2020 library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one qssp2020 spectral process in GiB. A RuntimeWarning is issued when the ranks on a node may not fit in its available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again. qssp2020 allocates its arrays from the input, so there is no fixed value; None skips the warning. Default: None.
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
        mpi4py is unavailable, or (on rank 0) spectral jobs still fail after the retries.
    ValueError
        There are fewer MPI ranks than the prepared group width.

    Notes
    -----
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. This computes the spectral stage only; run create_grnlib_qssp2020_func_parallel_multi_nodes afterwards. Ranks beyond the group width stay idle. A failed job does not stop the others: after all ranks finish, the jobs that did not complete are shared among the ranks and computed again, up to max_retries times. After an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    rank, problems = run_jobs_mpi(
        MPI,
        _tasks(path_green, "spec", check_finished),
        _run_job,
        memory_per_job_gb,
        exact_ranks=False,
        max_retries=max_retries,
    )
    if rank == 0:
        raise_if_incomplete(problems, path_green)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_qssp2020_func_parallel_multi_nodes(
    path_green, check_finished=False, memory_per_job_gb=None, max_retries=2
):
    """Compute the prepared qssp2020 library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one qssp2020 time-domain process in GiB. A RuntimeWarning is issued when the ranks on a node may not fit in its available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again. qssp2020 allocates its arrays from the input, so there is no fixed value; None skips the warning. Default: None.
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
        There are fewer MPI ranks than the prepared group width.

    Notes
    -----
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. Run it after the spectral stage has finished. The output stays in ASCII; convert it with convert_pd2bin_qssp2020_all. Ranks beyond the group width stay idle. A failed job does not stop the others: after all ranks finish, the jobs that did not complete are shared among the ranks and computed again, up to max_retries times. Rank 0 checks the library after all ranks finish; after an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    rank, problems = run_jobs_mpi(
        MPI,
        _tasks(path_green, "func", check_finished),
        _run_job,
        memory_per_job_gb,
        exact_ranks=False,
        max_retries=max_retries,
    )
    if rank == 0:
        _finish(path_green, problems, False, False)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def convert_pd2bin_qssp2020_all(path_green):
    """Convert all completed qssp2020 outputs to float32 binary libraries.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.

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
    See the qssp2020 tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units. Source/receiver depth pairs with incomplete output are skipped and keep their files, so they can be inspected and recomputed.
    """
    print("converting ascii files to byte files")
    green_info = _load_green_info(path_green)
    output_observables = np.nonzero(np.array(green_info["output_observables"]))[0]
    for event_dep in green_info["event_depth_list"]:
        for receiver_dep in green_info["receiver_depth_list"]:
            if any(
                _check_job(path_green, green_info, event_dep, receiver_dep, mt_com)
                for mt_com in mt_com_list
            ):
                continue
            for output_type_ind in output_observables:
                convert_pd2bin_qssp2020(
                    path_green, event_dep, receiver_dep, int(output_type_ind)
                )


def remove_dat_files(path_green):
    print("removing dat files")
    path_func = os.path.join(path_green, "GreenFunc")
    for root, dirs, files in os.walk(path_func):
        for file in glob.glob(os.path.join(root, "*.dat")):
            # keep output that was not converted
            if os.path.exists(file[:-4] + ".bin"):
                os.remove(file)


def check_grnlib_qssp2020(path_green, check_values=False):
    """Check that a qssp2020 library holds every file the readers need.

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
    For every source depth, receiver depth and moment-tensor component, each observable selected by output_observables must exist as a float32 .bin file of the size the readers expect, or as a complete .dat file with one row per sample and one column per distance. The P and S travel-time tables must hold one value per distance. Spectra under GreenSpec are not checked. Rerun the create_grnlib function with check_finished=True to recompute only incomplete jobs.
    """
    green_info = _load_green_info(path_green)
    dist_range = green_info["grn_dist_range"]
    n_dist = len(cal_grid(dist_range[0], dist_range[1], green_info["grn_delta_dist"]))
    problems = []
    for event_dep in green_info["event_depth_list"]:
        for receiver_dep in green_info["receiver_depth_list"]:
            sub_dir = os.path.join(
                path_green, "GreenFunc", "%.2f" % event_dep, "%.2f" % receiver_dep
            )
            for name in ["tp_table.bin", "ts_table.bin"]:
                problem = check_file_size(os.path.join(sub_dir, name), 4 * n_dist)
                problems += [problem] if problem else []
            for mt_com in mt_com_list:
                problems += _check_job(
                    path_green, green_info, event_dep, receiver_dep, mt_com, check_values
                )
    return problems


if __name__ == "__main__":
    pass
