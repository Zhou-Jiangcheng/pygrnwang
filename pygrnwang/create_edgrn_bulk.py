import os
import pickle
import json
import datetime

from .create_edgrn import (
    EDGRN2_NRMAX,
    EDGRN2_NZSMAX,
    EDGRN2_MEMORY_PER_JOB_GB,
    create_inp_edgrn2,
    call_edgrn2,
    check_output_edgrn2,
)
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


def _job_dir(path_green, obs_depth):
    return str(os.path.join(path_green, "edgrn2", "%.2f" % obs_depth))


def _run_job(task):
    """Run one job unless check_finished finds it complete; return its problems."""
    obs_depth, path_green, check_finished = task
    path_obs_dep = _job_dir(path_green, obs_depth)
    return run_checked_job(
        path_obs_dep,
        check_finished,
        lambda: call_edgrn2(obs_depth, path_green),
        lambda: check_output_edgrn2(path_obs_dep),
        # outputs of an earlier run would pass the check if this run fails
        stale=["edgrn.ss", "edgrn.ds", "edgrn.cl"],
    )


def _tasks(path_green, check_finished):
    with open(os.path.join(path_green, "group_list_edgrn.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    return [[(obs_depth, path_green, check_finished) for obs_depth in grp]
            for grp in group_list]


def _finish(path_green, run_problems):
    finish_library(path_green, run_problems, lambda: check_grnlib_edgrn2(path_green))


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_edgrn2(
    processes_num,
    path_green,
    grn_source_depth_range,
    grn_source_delta_depth,
    grn_dist_range,
    grn_delta_dist,
    obs_depth_list,
    wavenumber_sampling_rate=12,
    path_nd=None,
    earth_model_layer_num=None,
):
    # print("preprocessing edgrn2")

    """Prepare the edgrn library grid, input files and job groups.

    Parameters
    ----------
    processes_num : int
        Positive worker count used to group jobs; MPI rank count must match the prepared group width.
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    grn_source_depth_range : list of float
        Minimum and maximum source depths in km.
    grn_source_delta_depth : float
        Positive source-depth grid increment in km.
    grn_dist_range : list of float
        Minimum and maximum epicentral distances in km.
    grn_delta_dist : float
        Positive regular epicentral-distance increment in km.
    obs_depth_list : list of float
        Nonempty receiver depth list in km, positive down.
    wavenumber_sampling_rate : float, optional
        Dimensionless spatial Nyquist oversampling factor for wavenumber integration. Default: 12.
    path_nd : str or None, optional
        Six-column named-discontinuity model path: depth (km), Vp/Vs (km/s), density (g/cm3), Qp/Qs. Bulk preprocessing requires a real path even though the signature default is None. Default: None.
    earth_model_layer_num : int or None, optional
        Number of numeric model rows retained, not the number of discontinuities; None retains all. Default: None.

    Returns
    -------
    group_list : list
        Jobs grouped by processes_num; the same groups are saved as a pickle file.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    ValueError
        The grid has fewer than 2 or more than 10000 distances, fewer than 2 or more than 401 source depths (the limits of edgrn2 and edcmp2), or an output path is longer than the 160 characters edgrn2 can hold.

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions.
    """
    n_dist = len(cal_grid(grn_dist_range[0], grn_dist_range[1], grn_delta_dist))
    n_source_depth = len(
        cal_grid(
            grn_source_depth_range[0], grn_source_depth_range[1], grn_source_delta_depth
        )
    )
    if not 2 <= n_dist <= EDGRN2_NRMAX:
        raise ValueError(
            "The grid has %d distances; edgrn2 needs 2 to %d" % (n_dist, EDGRN2_NRMAX)
        )
    if not 2 <= n_source_depth <= EDGRN2_NZSMAX:
        raise ValueError(
            "The grid has %d source depths; edgrn2 needs 2 to %d"
            % (n_source_depth, EDGRN2_NZSMAX)
        )
    check_path_lengths(
        os.path.join(_job_dir(path_green, obs_depth), "edgrn.ss")
        for obs_depth in obs_depth_list
    )
    for obs_depth in obs_depth_list:
        sub_sub_dir = str(os.path.join(path_green, "edgrn2", "%.2f" % obs_depth))
        os.makedirs(sub_sub_dir, exist_ok=True)
        create_inp_edgrn2(
            path_green,
            obs_depth,
            grn_dist_range,
            grn_delta_dist,
            grn_source_depth_range,
            grn_source_delta_depth,
            wavenumber_sampling_rate,
            path_nd,
            earth_model_layer_num,
        )

    path_nd_without_Q = os.path.join(path_green, "noQ.nd")
    convert_earth_model_nd2nd_without_Q(path_nd, path_nd_without_Q)

    green_info = {
        "processes_num": processes_num,
        "grn_source_depth_range": grn_source_depth_range,
        "grn_source_delta_depth": grn_source_delta_depth,
        "grn_dist_range": grn_dist_range,
        "grn_delta_dist": grn_delta_dist,
        "obs_depth_list": obs_depth_list,
        "wavenumber_sampling_rate": wavenumber_sampling_rate,
        "path_nd": path_nd,
        "path_nd_without_Q": path_nd_without_Q,
        "earth_model_layer_num": earth_model_layer_num,
    }
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    group_list_edgrn = group(obs_depth_list, processes_num)
    with open(os.path.join(path_green, "group_list_edgrn.pkl"), "wb") as fw:
        pickle.dump(group_list_edgrn, fw)  # type: ignore
    return group_list_edgrn


def create_grnlib_edgrn2_sequential(path_green, check_finished=False, max_retries=2):
    """Compute the prepared edgrn library sequentially.

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
    elapsed : datetime.timedelta
        Wall-clock duration of the computation loop.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    RuntimeError
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_edgrn2).

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    problems = run_until_complete(
        lambda tasks: run_jobs_sequential(_run_job, tasks, desc="Computing Green's function library"),
        sum(_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems)
    e = datetime.datetime.now()
    print("run time:%s" % str(e - s))
    return e - s


def create_grnlib_edgrn2_parallel(
    path_green, check_finished=False, memory_per_job_gb=EDGRN2_MEMORY_PER_JOB_GB, max_retries=2
):
    """Compute the prepared edgrn library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one edgrn2 process in GiB. A RuntimeWarning is issued when min(processes_num, number of jobs) such processes may not fit in the currently available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again; None skips the warning. Default: 0.05, above the measured commit of the bundled executable.
    max_retries : int, optional
        Extra passes over the jobs that did not complete, for example because they ran out of memory; each pass recomputes only those jobs, with the same worker count. Default: 2.

    Returns
    -------
    elapsed : datetime.timedelta
        Wall-clock duration of the computation loop.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    RuntimeError
        Jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_edgrn2).

    Notes
    -----
    See the edgrn tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A job fails when the executable cannot start, exits with an error (for example after being killed for lack of memory) or leaves incomplete output files; it never gets a .finished marker. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs. On Windows call under an if __name__ == "__main__" guard.
    """
    s = datetime.datetime.now()
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        processes = json.load(fr).get("processes_num", None)
    problems = run_until_complete(
        lambda tasks: run_jobs_parallel(
            _run_job, tasks, processes, memory_per_job_gb, desc="Computing Green's function library"
        ),
        sum(_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems)
    e = datetime.datetime.now()
    return e - s


def create_grnlib_edgrn2_parallel_multi_nodes(
    path_green, check_finished=False, memory_per_job_gb=EDGRN2_MEMORY_PER_JOB_GB, max_retries=2
):
    """Compute the prepared edgrn library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one edgrn2 process in GiB. A RuntimeWarning is issued when the ranks on a node may not fit in its available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again; None skips the warning. Default: 0.05, above the measured commit of the bundled executable.
    max_retries : int, optional
        Extra passes over the jobs that did not complete, for example because they ran out of memory; each pass recomputes only those jobs, with the same worker count. Default: 2.

    Returns
    -------
    elapsed : datetime.timedelta
        Wall-clock duration of the computation loop.

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
    See the edgrn tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A failed job does not stop the others: after all ranks finish, the jobs that did not complete are shared among the ranks and computed again, up to max_retries times. Rank 0 checks the library after all ranks finish; after an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    rank, problems = run_jobs_mpi(
        MPI, _tasks(path_green, check_finished), _run_job, memory_per_job_gb,
        max_retries=max_retries,
    )
    if rank == 0:
        _finish(path_green, problems)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))
    return e - s


def check_grnlib_edgrn2(path_green):
    """Check that an edgrn library holds complete Green's function tables.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.

    Returns
    -------
    problems : list of str
        One line per missing or incomplete table; empty when the library is complete.

    Raises
    ------
    OSError
        green_lib_info.json cannot be read.

    Notes
    -----
    For every receiver depth, edgrn.ss, edgrn.ds and edgrn.cl must hold one row per distance and source depth listed in their own parameter line. Rerun the create_grnlib function with check_finished=True to recompute only incomplete jobs.
    """
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        obs_depth_list = json.load(fr)["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]
    problems = []
    for obs_depth in obs_depth_list:
        problems += check_output_edgrn2(_job_dir(path_green, obs_depth))
    return problems
