import os
import pickle
import json
import datetime

from tqdm import tqdm

from .create_spgrn2020 import create_dir_spgrn, check_output_spgrn
from .create_spgrn2020_bulk import (
    record_spgrn_distances,
    spgrn_job_dir,
    check_spgrn_paths,
    spgrn_tasks,
    check_spgrn_library,
)
from .create_spgrn2012 import create_inp_spgrn2012, call_spgrn2012
from .utils import (
    group,
    convert_earth_model_nd2nd_without_Q,
    check_file_size,
    run_checked_job,
    run_jobs_sequential,
    run_jobs_parallel,
    run_jobs_mpi,
    run_until_complete,
    finish_library,
)
from .pytaup import taup_create_npz_file, create_tpts_table


def _check_job(path_func, event_depth, expected_info=None, check_values=False):
    problems = check_output_spgrn(
        path_func, event_depth, expected_info=expected_info, check_values=check_values
    )
    if expected_info is not None:
        # the Python travel-time tables, one float32 per distance
        for name in ["tp_table.bin", "ts_table.bin"]:
            problem = check_file_size(
                os.path.join(path_func, name), 4 * len(expected_info["dist_list"])
            )
            problems += [problem] if problem else []
    return problems


def _run_job(task):
    """Run one job unless check_finished finds it complete; return its problems."""
    event_depth, receiver_depth, path_green, check_finished = task
    path_func = spgrn_job_dir(path_green, event_depth, receiver_depth)
    return run_checked_job(
        path_func,
        check_finished,
        lambda: call_spgrn2012(event_depth, receiver_depth, path_green),
        lambda: _check_job(path_func, event_depth),
        # outputs of an earlier run would pass the check if this run fails
        stale=["grn_d*", "GreenInfo*"],
    )


def _finish(path_green, run_problems):
    """Record the distances, write travel-time tables, check the library."""
    green_info = record_spgrn_distances(path_green)
    if green_info is not None:
        npz_file = taup_create_npz_file(nd_file=green_info["path_nd_without_Q"])
        for event_depth in tqdm(
            green_info["event_depth_list"], desc="Creating travel time tables"
        ):
            for receiver_depth in green_info["receiver_depth_list"]:
                create_tpts_table(
                    os.path.join(path_green, "GreenFunc"),
                    event_depth,
                    receiver_depth,
                    green_info["dist_list"],
                    npz_file,
                    False,
                )
    finish_library(path_green, run_problems, lambda: check_grnlib_spgrn2012(path_green))


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_spgrn2012(
    processes_num,
    path_green,
    event_depth_list,
    receiver_depth_list,
    spec_time_window,
    sampling_interval,
    max_frequency,
    max_slowness,
    anti_alias,
    gravity_fc,
    gravity_harmonic,
    cal_sph,
    cal_tor,
    source_radius,
    cal_gf,
    time_window,
    t0,
    v0,
    source_duration,
    dist_range,
    delta_dist_range,
    path_nd=None,
    earth_model_layer_num=None,
    physical_dispersion=0,
):
    """Prepare the spgrn2012 library grid, input files and job groups.

    .. warning::

        SPGRN2012 is deprecated. Use SPGRN2020 for new calculations.
        Rebuild the Green library with the replacement backend and revalidate
        numerical settings, source and time conventions, and results; existing
        libraries and parameter choices are not guaranteed to be interchangeable.

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
    gravity_fc : float
        Critical frequency in Hz below which self-gravity is included together with the harmonic cutoff.
    gravity_harmonic : int
        Critical spherical harmonic degree for self-gravity.
    cal_sph : int
        1 enables spheroidal (P-SV) modes; 0 disables them.
    cal_tor : int
        1 enables toroidal (SH) modes; 0 disables them.
    source_radius : float
        Source patch radius in km.
    cal_gf : int
        1 updates source spectra; 0 reuses spectra with identical model and spectral parameters.
    time_window : float
        Output time-window duration in seconds.
    t0 : float
        Reference start time in seconds in the SPGRN2012 rule t_start = t0 + distance/v0.
    v0 : float
        SPGRN2012 reduction velocity in km/s; use a nonzero value consistent with the desired phase window.
    source_duration : float
        Squared half-sinusoid source-time-function duration in seconds; zero requests the backend minimum.
    dist_range : list of float
        Minimum and maximum epicentral distances in km.
    delta_dist_range : list of float
        Smallest and largest distance increments in km at the near and far limits; equal values produce a regular grid.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.
    physical_dispersion : int, optional
        0 disables, 1 enables the backend physical-dispersion correction associated with attenuation. Default: 0.

    Returns
    -------
    group_list : list
        Jobs grouped by processes_num; the same groups are saved as a pickle file.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    ValueError
        An input or output path is longer than the 160 characters spgrn2012 can hold.

    Notes
    -----
    See the spgrn2012 tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions.
    """
    check_spgrn_paths(path_green, event_depth_list, receiver_depth_list)
    item_list = []
    for event_depth in event_depth_list:
        for receiver_depth in receiver_depth_list:
            create_dir_spgrn(event_depth, receiver_depth, path_green)
            create_inp_spgrn2012(
                path_green,
                event_depth,
                receiver_depth,
                spec_time_window,
                sampling_interval,
                max_frequency,
                max_slowness,
                anti_alias,
                gravity_fc,
                gravity_harmonic,
                cal_sph,
                cal_tor,
                source_radius,
                cal_gf,
                time_window,
                t0,
                v0,
                source_duration,
                dist_range,
                delta_dist_range,
                path_nd,
                earth_model_layer_num,
                physical_dispersion,
            )
            item_list.append([event_depth, receiver_depth])
    group_list = group(item_list, processes_num)
    with open(os.path.join(path_green, "group_list.pkl"), "wb") as fw:
        pickle.dump(group_list, fw)  # type: ignore

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    params = {
        "processes_num": processes_num,
        "path_green": path_green,
        "event_depth_list": event_depth_list,
        "receiver_depth_list": receiver_depth_list,
        "spec_time_window": spec_time_window,
        "sampling_interval": sampling_interval,
        "max_frequency": max_frequency,
        "max_slowness": max_slowness,
        "anti_alias": anti_alias,
        "gravity_fc": gravity_fc,
        "gravity_harmonic": gravity_harmonic,
        "cal_sph": cal_sph,
        "cal_tor": cal_tor,
        "source_radius": source_radius,
        "cal_gf": cal_gf,
        "time_window": time_window,
        # "sampling_num": round(time_window / sampling_interval + 1),
        "t0": t0,
        "v0": v0,
        "source_duration": source_duration,
        "dist_range": dist_range,
        "delta_dist_range": delta_dist_range,
        "path_nd": path_nd,
        "earth_model_layer_num": earth_model_layer_num,
        "physical_dispersion": physical_dispersion,
        "path_nd_without_Q": path_nd_without_Q,
    }
    json_str = json.dumps(params, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    return group_list


def create_grnlib_spgrn2012_sequential(path_green, check_finished=False, max_retries=2):
    """Compute the prepared spgrn2012 library sequentially.

    .. warning::

        SPGRN2012 is deprecated. Use SPGRN2020 for new calculations.
        Rebuild the Green library with the replacement backend and revalidate
        numerical settings, source and time conventions, and results; existing
        libraries and parameter choices are not guaranteed to be interchangeable.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
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
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_spgrn2012).

    Notes
    -----
    See the spgrn2012 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. The travel-time tables are written after the run. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    problems = run_until_complete(
        lambda tasks: run_jobs_sequential(_run_job, tasks, desc="Computing SPGRN2012 library"),
        sum(spgrn_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems)
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))


def create_grnlib_spgrn2012_parallel(
    path_green, check_finished=False, memory_per_job_gb=None, max_retries=2
):
    """Compute the prepared spgrn2012 library with local worker processes.

    .. warning::

        SPGRN2012 is deprecated. Use SPGRN2020 for new calculations.
        Rebuild the Green library with the replacement backend and revalidate
        numerical settings, source and time conventions, and results; existing
        libraries and parameter choices are not guaranteed to be interchangeable.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one spgrn2012 process in GiB. A RuntimeWarning is issued when min(processes_num, number of jobs) such processes may not fit in the currently available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again. spgrn2012 allocates most arrays from the input, so there is no fixed value; None skips the warning. Default: None.
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
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_spgrn2012).

    Notes
    -----
    See the spgrn2012 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A job fails when the executable cannot start, exits with an error (for example after being killed or refused memory) or leaves incomplete output files; it never gets a .finished marker. The travel-time tables are written after the run. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs. On Windows call under an if __name__ == "__main__" guard.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        processes = json.load(fr).get("processes_num", None)
    problems = run_until_complete(
        lambda tasks: run_jobs_parallel(
            _run_job, tasks, processes, memory_per_job_gb, desc="Computing SPGRN2012 library"
        ),
        sum(spgrn_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def create_grnlib_spgrn2012_parallel_multi_nodes(
    path_green, check_finished=False, memory_per_job_gb=None, max_retries=2
):
    """Compute the prepared spgrn2012 library with MPI.

    .. warning::

        SPGRN2012 is deprecated. Use SPGRN2020 for new calculations.
        Rebuild the Green library with the replacement backend and revalidate
        numerical settings, source and time conventions, and results; existing
        libraries and parameter choices are not guaranteed to be interchangeable.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one spgrn2012 process in GiB. A RuntimeWarning is issued when the ranks on a node may not fit in its available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again. spgrn2012 allocates most arrays from the input, so there is no fixed value; None skips the warning. Default: None.
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
    See the spgrn2012 tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A failed job does not stop the others: after all ranks finish, the jobs that did not complete are shared among the ranks and computed again, up to max_retries times. After all ranks finish, rank 0 records the distances in green_lib_info.json, writes the travel-time tables and checks the library; after an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    rank, problems = run_jobs_mpi(
        MPI, spgrn_tasks(path_green, check_finished), _run_job, memory_per_job_gb,
        max_retries=max_retries,
    )
    if rank == 0:
        _finish(path_green, problems)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def check_grnlib_spgrn2012(path_green, check_values=False):
    """Check that a spgrn2012 library holds every file the readers need.

    .. warning::

        SPGRN2012 is deprecated. Use SPGRN2020 for new calculations.
        Rebuild the Green library with the replacement backend and revalidate
        numerical settings, source and time conventions, and results; existing
        libraries and parameter choices are not guaranteed to be interchangeable.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_values : bool, optional
        Also read every Green's function file and report NaN or infinite values; this reads the whole library. Default: False.

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
    For every source/receiver depth pair, GreenInfo must list the same distances and samples as green_lib_info.json, grn_d must hold every distance, and the P and S travel-time tables tp_table.bin and ts_table.bin must hold one value per distance. Spectra under GreenSpec are not checked. Rerun the create_grnlib function with check_finished=True to recompute only incomplete jobs.
    """
    return check_spgrn_library(path_green, _check_job, check_values)


if __name__ == "__main__":
    pass
