import os
import pickle
import json
import datetime
import warnings

import numpy as np

from .create_edcmp import (
    EDCMP2_NRECMAX,
    EDCMP2_MEMORY_PER_JOB_GB,
    output_name_list,
    output_cha_num,
    create_inp_edcmp2,
    call_edcmp2,
    check_output_edcmp2,
    convert_edcmp2,
)
from .create_edgrn_bulk import check_grnlib_edgrn2
from .utils import (
    group,
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
    write_bin_atomic,
)


def _load_green_info(path_green):
    with open(
        os.path.join(path_green, "green_lib_info.json"), "r", encoding="utf-8"
    ) as fr:
        return json.load(fr)


def _job_dir(path_green, event_depth, obs_depth, mt_ind):
    return str(
        os.path.join(
            path_green, "edcmp2", "%.2f" % event_depth, "%.2f" % obs_depth, "%d" % mt_ind
        )
    )


def _grid(green_info):
    """Source depths, receiver depths and the number of distances of the library."""
    depth_range = green_info["grn_source_depth_range"]
    event_depth_list = cal_grid(
        depth_range[0], depth_range[1], green_info["grn_source_delta_depth"]
    )
    obs_depth_list = green_info["obs_depth_list"]
    if not isinstance(obs_depth_list, list):
        obs_depth_list = [obs_depth_list]
    dist_range = green_info["grn_dist_range"]
    n_dist = len(cal_grid(dist_range[0], dist_range[1], green_info["grn_delta_dist"]))
    return event_depth_list, obs_depth_list, n_dist


def _check_job(path_green, green_info, event_depth, obs_depth, mt_ind,
               check_values=False):
    return check_output_edcmp2(
        _job_dir(path_green, event_depth, obs_depth, mt_ind),
        green_info["output_observables"],
        _grid(green_info)[2],
        check_values,
    )


def _run_job(task):
    """Run one job unless check_finished finds it complete; return its problems."""
    event_depth, obs_depth, mt_ind, path_green, check_finished = task
    green_info = _load_green_info(path_green)
    return run_checked_job(
        _job_dir(path_green, event_depth, obs_depth, mt_ind),
        check_finished,
        lambda: call_edcmp2(event_depth, obs_depth, mt_ind, path_green),
        lambda: _check_job(path_green, green_info, event_depth, obs_depth, mt_ind),
        # outputs of an earlier run would pass the check if this run fails
        stale=["hs.*", "*.bin*"],
    )


def _tasks(path_green, check_finished):
    """Job groups, after checking that a layered model has its edgrn2 tables."""
    if _load_green_info(path_green).get("layered", True):
        raise_if_incomplete(check_grnlib_edgrn2(path_green), path_green)
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "rb") as fr:
        group_list = pickle.load(fr)
    return [[tuple(job + [path_green, check_finished]) for job in grp] for grp in group_list]


def _finish(path_green, run_problems, convert_bulk, remove):
    """Convert the jobs, check the whole library and raise if incomplete."""
    # the combined files need every job
    if convert_bulk and not run_problems:
        convert_pd2bin_edcmp2_all(path_green, remove=remove)
    finish_library(path_green, run_problems, lambda: check_grnlib_edcmp2(path_green))


def _get_mpi():
    try:
        from mpi4py import MPI
    except ImportError as exc:
        raise RuntimeError("mpi4py is required for multi-node MPI mode") from exc
    return MPI


def pre_process_edcmp2(
    processes_num: int,
    path_green: str,
    grn_source_depth_range,
    grn_source_delta_depth,
    grn_dist_range,
    grn_delta_dist,
    obs_depth_list,
    output_observables=(1, 0, 0, 0),
    layered=True,
    lam=30516224000,
    mu=33701888000,
):
    """Prepare the edcmp library grid, input files and job groups.

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
    output_observables : sequence of int, optional
        Four 0/1 flags in displacement, strain, stress, tilt order. Default: (1, 0, 0, 0).
    layered : bool, optional
        True uses the EDGRN layered-medium library; False selects the homogeneous half-space formula. Default: True.
    lam : float, optional
        First Lame parameter in Pa for the homogeneous half-space. Default: 30516224000.
    mu : float, optional
        Shear modulus in Pa for the homogeneous half-space. Default: 33701888000.

    Returns
    -------
    group_list : list
        Jobs grouped by processes_num; the same groups are saved as a pickle file.

    Raises
    ------
    OSError
        Required files are missing or output paths cannot be read or written.
    ValueError
        The grid has more than 500000 distances, the receiver limit of edcmp2, or an input, output or Green's function path is longer than the 160 characters edcmp2 can hold.

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Preprocessing writes inputs and travel-time/model metadata; run the matching create_grnlib function to calculate Green functions. Run EDGRN preparation first: this function reads and updates its metadata, even for homogeneous-half-space mode.
    """
    # Load the current Green's function library information.
    with open(os.path.join(path_green, "green_lib_info.json"), "r") as fr:
        green_info = json.load(fr)

    item_list = []
    event_depth_list = cal_grid(
        grn_source_depth_range[0],
        grn_source_depth_range[1],
        grn_source_delta_depth,
    )
    n_dist = len(cal_grid(grn_dist_range[0], grn_dist_range[1], grn_delta_dist))
    if n_dist > EDCMP2_NRECMAX:
        raise ValueError(
            "The grid has %d distances; edcmp2 computes at most %d receivers"
            % (n_dist, EDCMP2_NRECMAX)
        )
    paths = []
    for event_depth in event_depth_list:
        for obs_depth in obs_depth_list:
            path_sub_dir = _job_dir(path_green, event_depth, obs_depth, 4)
            paths += [
                os.path.join(path_sub_dir, "grn.inp"),
                os.path.join(path_sub_dir, "hs.strain"),
            ]
            if layered:
                paths.append(
                    os.path.join(path_green, "edgrn2", "%.2f" % obs_depth, "edgrn.ss")
                )
    check_path_lengths(paths)
    for event_depth in event_depth_list:
        for obs_depth in obs_depth_list:
            for mt_ind in range(5):
                create_inp_edcmp2(
                    path_green=path_green,
                    event_depth=event_depth,
                    obs_depth=obs_depth,
                    dist_range=grn_dist_range,
                    delta_dist=grn_delta_dist,
                    mt_ind=mt_ind,
                    output_observables=output_observables,
                    layered=layered,
                    lam=lam,
                    mu=mu,
                )
                item_list.append([event_depth, obs_depth, mt_ind])

    # Update the green_info dictionary with new observation and model parameters.
    # seek_edcmp2 indexes the edcmp2 output on the grid used here, so record it;
    # pre_process_edgrn2 wrote the edgrn2 grid into the same keys, and the two
    # disagreeing means the caller passed different parameters to the two steps.
    grid_keys = {
        "grn_source_depth_range": list(grn_source_depth_range),
        "grn_source_delta_depth": grn_source_delta_depth,
        "grn_dist_range": list(grn_dist_range),
        "grn_delta_dist": grn_delta_dist,
        "obs_depth_list": list(obs_depth_list),
    }
    for key, value in grid_keys.items():
        old = green_info.get(key, None)
        if old is not None and old != value:
            warnings.warn(
                "pre_process_edcmp2 got %s=%r but pre_process_edgrn2 recorded %r; "
                "the edcmp2 grid is the one seek_edcmp2 will use." % (key, value, old)
            )
        green_info[key] = value
    green_info["layered"] = layered
    if layered:
        green_info["lam"] = None
        green_info["mu"] = None
    else:
        green_info["lam"] = lam
        green_info["mu"] = mu
    green_info["output_observables"] = output_observables
    json_str = json.dumps(green_info, indent=4, ensure_ascii=False)
    with open(
        os.path.join(path_green, "green_lib_info.json"), "w", encoding="utf-8"
    ) as file:
        file.write(json_str)

    # Group the observation depths into subgroups for parallel processing.
    group_list_edcmp = group(item_list, processes_num)
    with open(os.path.join(path_green, "group_list_edcmp.pkl"), "wb") as fw:
        pickle.dump(group_list_edcmp, fw)  # type: ignore

    return group_list_edcmp


def create_grnlib_edcmp2_sequential(path_green, check_finished=False, max_retries=2):
    """Compute the prepared edcmp library sequentially.

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
        The layered model's edgrn library is incomplete, jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_edcmp2).

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. The output stays in ASCII; convert it with convert_pd2bin_edcmp2_all. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    problems = run_until_complete(
        lambda tasks: run_jobs_sequential(_run_job, tasks, desc="Computing static stress"),
        sum(_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems, False, False)


def create_grnlib_edcmp2_parallel(
    path_green,
    check_finished=False,
    convert_bulk=True,
    remove=False,
    memory_per_job_gb=EDCMP2_MEMORY_PER_JOB_GB,
    max_retries=2,
):
    """Compute the prepared edcmp library with local worker processes.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    convert_bulk : bool, optional
        Generate combined EDCMP float32 libraries after jobs finish. Default: True.
    remove : bool, optional
        Delete source ASCII files after converting them; keep False when inspecting backend output. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one edcmp2 process in GiB. A RuntimeWarning is issued when min(processes_num, number of jobs) such processes may not fit in the currently available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again; None skips the warning. Default: 1.0, the measured commit of the bundled executable.
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
        The layered model's edgrn library is incomplete, jobs still fail after the retries, or the finished library is incomplete (see check_grnlib_edcmp2).

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. A job fails when the executable cannot start, exits with an error (for example after being killed for lack of memory) or leaves incomplete output files; it never gets a .finished marker. A failed job does not stop the others: afterwards the jobs that did not complete are computed again, up to max_retries times, and the library is checked. Ctrl+C stops the run at once and kills the running backends. After an error, rerun with check_finished=True to compute only unfinished jobs. On Windows call under an if __name__ == "__main__" guard.
    """
    s = datetime.datetime.now()
    processes = _load_green_info(path_green).get("processes_num", None)
    problems = run_until_complete(
        lambda tasks: run_jobs_parallel(
            _run_job, tasks, processes, memory_per_job_gb, desc="Computing static lib"
        ),
        sum(_tasks(path_green, check_finished), []),
        max_retries,
    )
    _finish(path_green, problems, convert_bulk, remove)
    e = datetime.datetime.now()
    return e - s


def create_grnlib_edcmp2_parallel_multi_nodes(
    path_green, check_finished=False, memory_per_job_gb=EDCMP2_MEMORY_PER_JOB_GB, max_retries=2
):
    """Compute the prepared edcmp library with MPI.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_finished : bool, optional
        Reuse jobs marked finished whose output is complete; recompute the others. Markers do not verify that inputs are unchanged. Default: False.
    memory_per_job_gb : float or None, optional
        Peak memory of one edcmp2 process in GiB. A RuntimeWarning is issued when the ranks on a node may not fit in its available memory (including Slurm/container cgroup limits); the run goes on and jobs that run out of memory are computed again; None skips the warning. Default: 1.0, the measured commit of the bundled executable.
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
        mpi4py is unavailable, the layered model's edgrn library is incomplete, or (on rank 0) jobs still fail after the retries or the finished library is incomplete.
    ValueError
        MPI rank count does not match the prepared group width.

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Prepare jobs first. Backend runners can change the process working directory; use absolute paths and restore the caller directory if needed. The output stays in ASCII; convert it with convert_pd2bin_edcmp2_all. A failed job does not stop the others: after all ranks finish, the jobs that did not complete are shared among the ranks and computed again, up to max_retries times. Rank 0 checks the library after all ranks finish; after an error, rerun with check_finished=True to compute only unfinished jobs.
    """
    s = datetime.datetime.now()
    MPI = _get_mpi()
    rank, problems = run_jobs_mpi(
        MPI, _tasks(path_green, check_finished), _run_job, memory_per_job_gb,
        max_retries=max_retries,
    )
    if rank == 0:
        _finish(path_green, problems, False, False)
    e = datetime.datetime.now()
    print("run time:" + str(e - s))


def convert_pd2bin_edcmp2_all(path_green, remove=False):
    """Convert all completed edcmp outputs to float32 binary libraries.

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
    RuntimeError
        A job's output is incomplete; nothing is converted.

    Notes
    -----
    See the edcmp tutorial for a complete prepare, run and read workflow. Conversion is a storage operation; it does not resample or change physical units.
    """
    print("converting ascii files to binary float32 files")
    green_info = _load_green_info(path_green)
    event_depth_list, obs_depth_list, n_dist = _grid(green_info)
    # the combined files need every job
    problems = []
    for event_depth in event_depth_list:
        for obs_depth in obs_depth_list:
            for k in range(5):
                problems += _check_job(path_green, green_info, event_depth, obs_depth, k)
    raise_if_incomplete(problems, path_green)

    output_observables = np.nonzero(np.array(green_info["output_observables"]))[0]
    n_dep = len(event_depth_list)
    n_obs = len(obs_depth_list)

    # Pre-allocate one bulk array per active output_type:
    # shape (n_dep, n_obs, 5, n_dist, cha_num)
    bulk_data = {}
    for o in output_observables:
        ot = output_name_list[int(o)]
        bulk_data[ot] = np.zeros(
            (n_dep, n_obs, 5, n_dist, output_cha_num[ot]), dtype=np.float32
        )

    for i, event_depth in enumerate(event_depth_list):
        for j, obs_depth in enumerate(obs_depth_list):
            for k in range(5):
                path_sub_dir = _job_dir(path_green, event_depth, obs_depth, k)
                for o in output_observables:
                    v_ijko = convert_edcmp2(
                        path_sub_dir=path_sub_dir,
                        output_type_ind=int(o),
                        remove=remove,
                    )
                    bulk_data[output_name_list[int(o)]][i, j, k] = v_ijko

    for ot, arr in bulk_data.items():
        out_path = os.path.join(path_green, "edcmp2_%s.bin" % ot)
        write_bin_atomic(arr, out_path)
        print("Saved %s, shape: %s" % (os.path.basename(out_path), str(arr.shape)))



def check_grnlib_edcmp2(path_green, check_values=False):
    """Check that an edcmp library holds every file the readers need.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    check_values : bool, optional
        Also read every binary file and report NaN or infinite values; this reads the whole library. Default: False.

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
    For every source depth, receiver depth and base mechanism, each output selected by output_observables must exist as <output>.bin with one row per distance, or as a complete hs.<output> file. Combined edcmp2_<output>.bin files, when present, must hold every job. The edgrn tables are checked by check_grnlib_edgrn2. Rerun the create_grnlib function with check_finished=True to recompute only incomplete jobs.
    """
    green_info = _load_green_info(path_green)
    event_depth_list, obs_depth_list, n_dist = _grid(green_info)
    problems = []
    for event_depth in event_depth_list:
        for obs_depth in obs_depth_list:
            for k in range(5):
                problems += _check_job(
                    path_green, green_info, event_depth, obs_depth, k, check_values
                )
    for ind, selected in enumerate(green_info["output_observables"]):
        path_bulk = os.path.join(
            path_green, "edcmp2_%s.bin" % output_name_list[ind]
        )
        if selected and os.path.exists(path_bulk):
            problem = check_file_size(
                path_bulk,
                len(event_depth_list) * len(obs_depth_list) * 5 * n_dist * 4
                * output_cha_num[output_name_list[ind]],
                check_values,
            )
            problems += [problem] if problem else []
    return problems
