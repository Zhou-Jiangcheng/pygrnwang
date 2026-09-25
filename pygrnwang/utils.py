import os
import sys
import math
import platform
import subprocess

import numpy as np
import pandas as pd

from .signal_process import linear_interp
from .pytaup import cal_first_p_s
from . import geo


# These three live in geo; they used to be duplicated here verbatim.
# Re-exported so the historical pygrnwang.utils import paths keep working.
cal_max_dist_from_2d_points = geo.cal_max_dist_from_2d_points
create_rotate_z_mat = geo.create_rotate_z_mat
rotate_symmetric_tensor_series = geo.rotate_symmetric_tensor_series


def read_source_array(source_inds, path_input, shift2corner=False, source_shapes=None):
    source_array = None
    for ind_src in range(len(source_inds)):
        source_plane = pd.read_csv(
            str(os.path.join(path_input, "source_plane%d.csv" % source_inds[ind_src])),
            index_col=False,
            header=None,
        ).to_numpy()
        if shift2corner:
            mu_strike = (
                source_plane[source_shapes[ind_src][1], :3] - source_plane[0, :3]
            )
            mu_dip = source_plane[1, :3] - source_plane[0, :3]
            source_plane[:, :3] = source_plane[:, :3] - mu_strike / 2 - mu_dip / 2
        if ind_src == 0:
            source_array = source_plane.copy()
        else:
            source_array = np.concatenate([source_array, source_plane.copy()], axis=0)
    return source_array


def cal_grid(v_min, v_max, delta):
    """Construct the regular grid shared by writers and readers.

    Parameters
    ----------
    v_min : float
        First grid value in any consistent unit.
    v_max : float
        Minimum required terminal value, in the same unit as v_min.
    delta : float
        Positive grid spacing in the same unit as v_min.

    Returns
    -------
    grid : numpy.ndarray

    Raises
    ------
    ZeroDivisionError
        delta is zero.
    ValueError
        Non-finite grid parameters prevent calculation of the sample count.

    Notes
    -----
    """
    n = math.ceil((v_max - v_min) / delta) + 1
    return v_min + np.arange(n) * delta


def group(inp_list, num_in_each_group):
    group_list = []
    for i in range(len(inp_list) // num_in_each_group):
        group_list.append(inp_list[i * num_in_each_group : (i + 1) * num_in_each_group])
    rest = len(inp_list) % num_in_each_group
    if rest != 0:
        group_list.append(inp_list[-rest:])
    return group_list


def shift_green2real_tpts(
    seismograms,
    tpts_table,
    green_before_p,
    srate,
    event_depth_km,
    dist_in_km,
    receiver_depth_km=0,
    model_name="ak135",
):
    first_p, first_s = cal_first_p_s(
        event_depth_km=event_depth_km,
        dist_km=dist_in_km,
        receiver_depth_km=receiver_depth_km,
        model_name=model_name,
    )
    p_count = round(green_before_p * srate)
    s_count = round(
        (tpts_table["s_onset"] - tpts_table["p_onset"] + green_before_p) * srate
    )
    p_count_new = round((first_p - tpts_table["p_onset"] + green_before_p) * srate)
    s_count_new = min(
        len(seismograms[0]),
        round((first_s - tpts_table["p_onset"] + green_before_p) * srate),
    )
    if s_count == p_count or s_count_new == p_count_new:
        return seismograms, first_p, first_s

    n_samples = seismograms.shape[1]
    for i in range(seismograms.shape[0]):
        # own local name: green_before_p must stay the scalar parameter
        before_p_part = seismograms[i][:p_count]
        p_s = linear_interp(seismograms[i][p_count:s_count], s_count_new - p_count_new)
        after_s = seismograms[i][s_count:]
        if len(after_s) > 0:
            after_s = linear_interp(
                after_s, max(0, n_samples - len(before_p_part) - len(p_s))
            )
            row = np.concatenate([before_p_part, p_s, after_s])
        else:
            row = np.concatenate([before_p_part, p_s])
        # the pieces do not always add up to the original length, e.g. when after_s
        # is empty and the real S is earlier than the one in the library
        if len(row) < n_samples:
            row = np.concatenate([row, np.zeros(n_samples - len(row))])
        seismograms[i] = row[:n_samples]

    return seismograms, first_p, first_s


def convert_earth_model_nd2inp(path_nd, path_output):
    """Convert ND numeric rows to numbered backend model input lines.

    Parameters
    ----------
    path_nd : str
        Path to a named-discontinuity text model with either four or six numeric columns as required by the operation.
    path_output : str
        Destination path for the converted model.

    Returns
    -------
    lines : list of str
        Numeric rows prefixed by a one-based row number and terminated by newlines.

    Raises
    ------
    OSError
        The input model cannot be read.

    Notes
    -----
    path_output is retained for API compatibility but is not written by this function. Discontinuity labels are removed.
    """
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            lines_new.append(temp)
    for i in range(len(lines_new)):
        # print(lines_new[i])
        lines_new[i] = "  ".join([str(int(i + 1))] + lines_new[i]) + "\n"  # type: ignore
    # with open(path_output, "w") as fw:
    #     fw.writelines(lines_new)
    return lines_new


def convert_earth_model_nd2nd_without_Q(path_nd, path_output):
    """Write a four-column ND model by removing the final two Q columns.

    Parameters
    ----------
    path_nd : str
        Path to a named-discontinuity text model with either four or six numeric columns as required by the operation.
    path_output : str
        Destination path for the converted model.

    Returns
    -------
    lines : list of str
        The exact converted lines written to path_output, with labels retained.

    Raises
    ------
    OSError
        Input or output cannot be accessed.

    Notes
    -----
    The input must have six numeric columns. Passing a four-column file removes real physical columns and is invalid.
    """
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            lines_new.append([str(float(_)) for _ in temp[:-2]])
            lines_new[i] = "  ".join(lines_new[i]) + "\n"
        else:
            lines_new.append(lines[i].strip() + "\n")
    with open(path_output, "w") as fw:
        fw.writelines(lines_new)
    return lines_new


def read_nd(path_nd, with_Q=False):
    """Read numeric rows from a named-discontinuity model.

    Parameters
    ----------
    path_nd : str
        Path to a named-discontinuity text model with either four or six numeric columns as required by the operation.
    with_Q : bool, optional
        True expects exactly six numeric columns; False expects four. The flag describes the file, it does not remove Q columns. Default: False.

    Returns
    -------
    model : numpy.ndarray
        Shape (N, 4) or (N, 6): depth km, Vp/Vs km/s, density g/cm3, and
        optional dimensionless Qp/Qs. Single-token discontinuity labels are skipped.

    Raises
    ------
    OSError
        The model cannot be read.
    ValueError
        Numeric rows do not match the requested column layout.
    """
    with open(path_nd, "r") as fr:
        lines = fr.readlines()
    lines_new = []
    for i in range(len(lines)):
        temp = lines[i].split()
        if len(temp) > 1:
            for j in range(len(temp)):
                lines_new.append(float(temp[j]))
    if with_Q:
        nd_model = np.array(lines_new).reshape(-1, 6)
    else:
        nd_model = np.array(lines_new).reshape(-1, 4)
    return nd_model


def read_material_nd(model_name, depth):
    """Select the first material row at or below the requested depth.

    Parameters
    ----------
    model_name : str
        TauP built-in model name or path to a custom model. Use a model consistent with the Green library.
    depth : float
        Material lookup depth in km, positive down.

    Returns
    -------
    material : numpy.ndarray
        [depth km, Vp km/s, Vs km/s, density g/cm3]; below the model bottom,
        the final row is returned.

    Raises
    ------
    FileNotFoundError
        model_name is neither ak135fc nor an existing model file.

    Notes
    -----
    This is a row selection, not interpolation. model_name accepts only the built-in ak135fc or a four-column no-Q ND file path; ak135 is not a material-lookup alias.
    """
    if model_name == "ak135fc":
        from .ak135fc import s as str_nd

        lines = str_nd.split("\n")
        lines_new = []
        for i in range(len(lines)):
            temp = lines[i].split()
            if len(temp) > 1:
                for j in range(len(temp)):
                    lines_new.append(float(temp[j]))
        nd_model = np.array(lines_new).reshape(-1, 4)
    else:
        if not os.path.isfile(model_name):
            raise FileNotFoundError(
                "model_name must be the built-in 'ak135fc' or a path to an nd file, "
                "got %r" % (model_name,)
            )
        nd_model = read_nd(model_name)
    inds = np.argwhere((nd_model[:, 0] - depth) >= 0)
    # below the bottom of the model: use its deepest layer
    ind = inds[0][0] if len(inds) > 0 else len(nd_model) - 1
    return nd_model[ind]


def read_layerd_material(path_layerd_dat, depth_in_km):
    # thickness, rho, vp, vs, qp, qs
    """Read the material layer containing a depth from a thickness table.

    Parameters
    ----------
    path_layerd_dat : str
        Text table of thickness (m), density, Vp, Vs, Qp, Qs; returned material values retain the file units.
    depth_in_km : float
        Material lookup depth in km.

    Returns
    -------
    material : numpy.ndarray
        One row [thickness, density, Vp, Vs, Qp, Qs] in the file units;
        below the model bottom the final row is returned.

    Raises
    ------
    OSError
        The table cannot be read.
    """
    depth_in_m = depth_in_km * 1e3
    dat = np.loadtxt(path_layerd_dat)
    inds = np.argwhere((np.cumsum(dat[:, 0]) - depth_in_m) >= 0)
    # below the bottom of the model: use its deepest layer
    ind = inds[0][0] if len(inds) > 0 else len(dat) - 1
    return dat[ind]


def create_stf(tau, srate):
    """Sample a normalized squared half-sine source-rate function.

    Parameters
    ----------
    tau : float
        Positive source duration in seconds.
    srate : float
        Positive output sampling rate in Hz.

    Returns
    -------
    stf : numpy.ndarray
        Shape (round(tau*srate)+1,), samples of 2/tau * sin(pi*t/tau)^2,
        in 1/s. Its continuous integral over [0, tau] is one.

    Notes
    -----
    The discrete integral depends on sampling; normalize explicitly when exact discrete convolution normalization is required.
    """
    t = np.linspace(0, tau, round(tau * srate) + 1, endpoint=True)
    stf = (2 / tau) * (np.sin(np.pi * t / tau)) ** 2
    return stf


def group_planes(strike_array):
    """
    It is necessary to ensure that the sub faults on
    the same fault plane have the same strikes!!!
    :param strike_array: numpy array

    Returns:
    np.array: An array containing the lengths of each group.
    """
    # Find the indices where the value changes
    # a[1:] != a[:-1] produces a boolean array that's True
    # at positions where a value differs from its predecessor.
    change_indices = np.where(strike_array[1:] != strike_array[:-1])[0] + 1

    # Include the start and end indices to get boundaries for each group.
    boundaries = np.concatenate(([0], change_indices, [len(strike_array)]))

    # The difference between consecutive boundaries gives the group lengths.
    lengths = np.diff(boundaries)
    return lengths


def reshape_sub_faults(sub_faults, num_strike, num_dip):
    mu_strike = sub_faults[num_dip] - sub_faults[0]
    mu_dip = sub_faults[1] - sub_faults[0]
    sub_faults = sub_faults - mu_strike / 2 - mu_dip / 2
    X: np.ndarray = sub_faults[:, 0]
    Y: np.ndarray = sub_faults[:, 1]
    Z: np.ndarray = sub_faults[:, 2]

    X = X.reshape(num_strike, num_dip)
    Y = Y.reshape(num_strike, num_dip)
    Z = Z.reshape(num_strike, num_dip)

    X = np.concatenate([X, np.array([X[:, -1] + mu_dip[0]]).T], axis=1)
    Y = np.concatenate([Y, np.array([Y[:, -1] + mu_dip[1]]).T], axis=1)
    Z = np.concatenate([Z, np.array([Z[:, -1] + mu_dip[2]]).T], axis=1)

    X = np.concatenate([X, np.array([X[-1, :] + mu_strike[0]])], axis=0)
    Y = np.concatenate([Y, np.array([Y[-1, :] + mu_strike[1]])], axis=0)
    Z = np.concatenate([Z, np.array([Z[-1, :] + mu_strike[2]])], axis=0)
    return X, Y, Z


def _fortran_stopped_with_error(stderr_text):
    # gfortran exits with code 0 after STOP 'message', and every backend uses a
    # STOP message only for errors; runtime errors print a recognizable prefix.
    for line in stderr_text.splitlines():
        line = line.strip()
        if line.startswith(("STOP", "ERROR STOP")) or "runtime error" in line.lower():
            return True
    return False


# set in worker processes of run_jobs_parallel; set() when the user stops the run
_stop_event = None


def _communicate(proc, data):
    """proc.communicate(data), killing proc when the run is interrupted.

    Ctrl+C (KeyboardInterrupt) or a stop request from run_jobs_parallel kills
    the backend process, so an interrupted run leaves no backend running.
    """
    try:
        while True:
            try:
                return proc.communicate(data, timeout=1)
            except subprocess.TimeoutExpired:
                data = None  # sent with the first call
                if _stop_event is not None and _stop_event.is_set():
                    raise KeyboardInterrupt("the run was stopped")
    except BaseException:
        if platform.system() == "Windows":
            # Scripts\<name>.exe is a launcher that starts Python, which starts
            # the Fortran executable; kill the whole tree, not only the launcher
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(proc.pid)], capture_output=True
            )
        proc.kill()
        proc.communicate()
        raise


def call_exe(path_inp, path_finished, name):
    """Run one backend executable and record whether it succeeded.

    Parameters
    ----------
    path_inp : str
        Absolute backend input file path, passed on standard input.
    path_finished : str
        Success marker path. The executable log is written here only on success;
        otherwise it is written to ``.failed`` in the same directory.
    name : str
        Executable name without the platform suffix.

    Returns
    -------
    ok : bool
        False when the executable could not start, exited with a nonzero code
        (for example after being killed for lack of memory) or stopped with a
        Fortran error message.
    """
    if platform.system() == "Windows":
        name_exe = "%s.exe" % name
        path_exe = os.path.join(sys.exec_prefix, "Scripts", name_exe)
    else:
        name_exe = "%s.bin" % name
        path_exe = os.path.join(sys.exec_prefix, "bin", name_exe)
    path_failed = os.path.join(os.path.dirname(path_finished), ".failed")
    # drop the status of an earlier run so a failure can never keep an old marker
    for path in (path_finished, path_failed):
        if os.path.exists(path):
            os.remove(path)
    try:
        proc = subprocess.Popen(
            [path_exe],
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
        )
    except OSError as exc:
        # e.g. Windows refuses to start a process when the commit limit is reached
        with open(path_failed, "w", encoding="utf-8") as fw:
            fw.write("%s could not run: %s\n" % (path_exe, exc))
        return False
    stdout_bytes, stderr_bytes = _communicate(proc, str.encode(path_inp))
    stdout_text = stdout_bytes.decode(errors="ignore")
    stderr_text = stderr_bytes.decode(errors="ignore")
    output = stdout_text + stderr_text
    ok = proc.returncode == 0 and not _fortran_stopped_with_error(stderr_text)
    with open(path_finished if ok else path_failed, "w", encoding="utf-8") as fw:
        fw.writelines(output)
        if not ok:
            fw.write("\nexit code: %d\n" % proc.returncode)
    return ok


def read_failure_log(path_dir, max_chars=300):
    """Return the error lines of the ``.failed`` log in path_dir, or an empty string."""
    path_failed = os.path.join(path_dir, ".failed")
    if not os.path.exists(path_failed):
        return ""
    with open(path_failed, "r", encoding="utf-8", errors="ignore") as fr:
        lines = [" ".join(line.split()) for line in fr if line.strip()]
    # the error message and exit code, not the banner or the backtrace
    errors = [
        line
        for line in lines
        if any(key in line for key in ("STOP", "rror", "exit code", "could not run"))
        and "Backtrace" not in line
        and not line.startswith("#")
    ]
    return " | ".join(errors or lines)[-max_chars:]


def _cgroup_free_memory_bytes():
    # Slurm and containers enforce memory with cgroups; /proc/meminfo shows the node.
    try:
        with open("/proc/self/cgroup", "r") as fr:
            entries = [line.rstrip("\n").split(":", 2) for line in fr]
    except OSError:
        return None
    for entry in entries:
        if len(entry) != 3:
            continue
        _, controllers, path = entry
        if controllers == "":  # cgroup v2
            base = "/sys/fs/cgroup" + path
            limit_name, usage_name = "memory.max", "memory.current"
        elif "memory" in controllers.split(","):  # cgroup v1
            base = "/sys/fs/cgroup/memory" + path
            limit_name, usage_name = "memory.limit_in_bytes", "memory.usage_in_bytes"
        else:
            continue
        try:
            with open(os.path.join(base, limit_name), "r") as fr:
                limit = fr.read().strip()
            with open(os.path.join(base, usage_name), "r") as fr:
                usage = int(fr.read().strip())
        except (OSError, ValueError):
            continue
        if limit == "max" or int(limit) >= 1 << 60:
            return None
        return max(0, int(limit) - usage)
    return None


def available_memory_bytes():
    """Estimate the memory new backend processes can use now.

    Returns
    -------
    available : int or None
        Bytes, or None when the platform offers no estimate. On Windows this is
        the smaller of free physical memory and free commit charge. On Linux it
        is MemAvailable, further limited by a cgroup (Slurm, container) limit.
    """
    if platform.system() == "Windows":
        import ctypes

        class MEMORYSTATUSEX(ctypes.Structure):
            _fields_ = [
                ("dwLength", ctypes.c_ulong),
                ("dwMemoryLoad", ctypes.c_ulong),
                ("ullTotalPhys", ctypes.c_ulonglong),
                ("ullAvailPhys", ctypes.c_ulonglong),
                ("ullTotalPageFile", ctypes.c_ulonglong),
                ("ullAvailPageFile", ctypes.c_ulonglong),
                ("ullTotalVirtual", ctypes.c_ulonglong),
                ("ullAvailVirtual", ctypes.c_ulonglong),
                ("ullAvailExtendedVirtual", ctypes.c_ulonglong),
            ]

        stat = MEMORYSTATUSEX()
        stat.dwLength = ctypes.sizeof(MEMORYSTATUSEX)
        if not ctypes.windll.kernel32.GlobalMemoryStatusEx(ctypes.byref(stat)):
            return None
        return min(stat.ullAvailPhys, stat.ullAvailPageFile)

    available = None
    try:
        with open("/proc/meminfo", "r") as fr:
            for line in fr:
                if line.startswith("MemAvailable:"):
                    available = int(line.split()[1]) * 1024
                    break
    except (OSError, ValueError):
        pass
    if available is None:
        try:
            available = os.sysconf("SC_AVPHYS_PAGES") * os.sysconf("SC_PAGE_SIZE")
        except (ValueError, OSError, AttributeError):
            pass
    cgroup_free = _cgroup_free_memory_bytes()
    if cgroup_free is not None:
        available = cgroup_free if available is None else min(available, cgroup_free)
    return available


def warn_if_memory_short(processes, memory_per_job_gb, where=""):
    """Warn when the backend processes may not fit in the available memory.

    The run goes on: jobs that fail for lack of memory are computed again
    after the others (see run_until_complete).

    Parameters
    ----------
    processes : int
        Number of backend processes that will run at the same time.
    memory_per_job_gb : float or None
        Peak memory of one backend process in GiB; None skips the check.
    where : str, optional
        Prefix of the warning, e.g. the node name. Default: "".

    Returns
    -------
    None
    """
    if not memory_per_job_gb:
        return
    available = available_memory_bytes()
    if available is None:
        return
    required = processes * memory_per_job_gb * 2**30
    if required > available:
        import warnings

        warnings.warn(
            "%s%d backend processes need about %.1f GiB (%.1f GiB each) but only "
            "%.1f GiB is available; jobs that run out of memory will be computed "
            "again after the others. At most %d processes fit."
            % (
                where,
                processes,
                required / 2**30,
                memory_per_job_gb,
                available / 2**30,
                int(available / (memory_per_job_gb * 2**30)),
            ),
            RuntimeWarning,
            stacklevel=2,
        )


def _init_worker(stop_event):
    global _stop_event
    _stop_event = stop_event
    # Ctrl+C reaches the main process, which stops the workers through
    # stop_event; a KeyboardInterrupt inside a worker would only break the pool
    import signal

    signal.signal(signal.SIGINT, signal.SIG_IGN)


def _run_unless_stopped(run_job, task):
    # the executor hands a few jobs to workers in advance; they must not start
    # after the user stopped the run
    if _stop_event.is_set():
        return []
    return run_job(task)


def _report_failure(problems):
    from tqdm import tqdm

    tqdm.write("A job failed; continuing with the other jobs:\n  %s" % problems[0])


def run_jobs_parallel(run_job, tasks, processes, memory_per_job_gb=None, desc=""):
    """Run backend jobs in worker processes; a failed job does not stop the others.

    Parameters
    ----------
    run_job : callable
        Module-level function taking one task and returning a list of problem
        strings; an empty list means the job completed.
    tasks : list of tuple
        Picklable, hashable job arguments.
    processes : int or None
        Worker count; None uses the CPU count.
    memory_per_job_gb : float or None, optional
        Peak memory of one backend process in GiB; a RuntimeWarning is issued
        when the workers may not fit in the available memory, and the run goes
        on. None skips the check. Default: None.
    desc : str, optional
        Progress-bar label.

    Returns
    -------
    failed : dict
        Task to problem list for every job that did not complete.

    Raises
    ------
    KeyboardInterrupt
        The user pressed Ctrl+C; no new job starts and running backends are killed.
    """
    import multiprocessing
    from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
    from tqdm import tqdm

    workers = min(processes or os.cpu_count() or 1, len(tasks))
    if workers == 0:
        return {}
    warn_if_memory_short(workers, memory_per_job_gb)
    print("Running %d jobs with %d worker processes" % (len(tasks), workers))
    failed = {}
    stop_event = multiprocessing.Event()
    # a worker killed by the OS breaks the executor instead of hanging it; its
    # unfinished jobs then report BrokenProcessPool and count as failed
    executor = ProcessPoolExecutor(
        max_workers=workers, initializer=_init_worker, initargs=(stop_event,)
    )
    try:
        futures = {
            executor.submit(_run_unless_stopped, run_job, task): task for task in tasks
        }
        pending = set(futures)
        with tqdm(total=len(futures), desc=desc) as bar:
            while pending:
                # a timed wait, because Ctrl+C cannot interrupt an endless one on Windows
                done, pending = wait(pending, timeout=1, return_when=FIRST_COMPLETED)
                for future in done:
                    bar.update(1)
                    try:
                        problems = future.result()
                    except Exception as exc:
                        problems = [
                            "job %s: %s: %s"
                            % (futures[future], type(exc).__name__, exc)
                        ]
                    if problems:
                        _report_failure(problems)
                        failed[futures[future]] = problems
    except KeyboardInterrupt:
        print("\nInterrupted; stopping the running jobs", flush=True)
        stop_event.set()
        executor.shutdown(wait=True, cancel_futures=True)
        raise
    executor.shutdown(wait=True)
    return failed


def run_jobs_sequential(run_job, tasks, desc=""):
    """Run backend jobs one after another; a failed job does not stop the others.

    Returns the task to problem list for every job that did not complete.
    """
    from tqdm import tqdm

    failed = {}
    for task in tqdm(tasks, desc=desc):
        problems = run_job(task)
        if problems:
            _report_failure(problems)
            failed[task] = problems
    return failed


def run_until_complete(run_pass, tasks, max_retries=2):
    """Run all jobs, then run the failed ones again until they complete.

    Parameters
    ----------
    run_pass : callable
        Runs a list of tasks and returns the task to problem list of the failed ones.
    tasks : list of tuple
        All jobs.
    max_retries : int, optional
        Number of extra passes over the jobs that did not complete, for
        example because they ran out of memory. Default: 2.

    Returns
    -------
    problems : list of str
        Problems of the jobs still failing after the last pass.
    """
    failed = run_pass(tasks)
    for attempt in range(1, max_retries + 1):
        if not failed:
            break
        print(
            "%d jobs did not complete; computing them again (retry %d/%d)"
            % (len(failed), attempt, max_retries)
        )
        failed = run_pass(list(failed))
    return sum(failed.values(), [])


def run_jobs_mpi(
    MPI,
    group_list,
    run_job,
    memory_per_job_gb=None,
    exact_ranks=True,
    max_retries=2,
):
    """Run the prepared job groups on MPI ranks, then recompute failed jobs.

    Parameters
    ----------
    MPI : module
        mpi4py.MPI.
    group_list : list of list of tuple
        Job groups; rank r runs task r of every group that has one.
    run_job : callable
        Function taking one task and returning a list of problem strings.
    memory_per_job_gb : float or None, optional
        Peak memory of one backend process in GiB; one rank per node issues a
        RuntimeWarning when the ranks of its node may not fit in the available
        memory, and the run goes on. None skips the check. Default: None.
    exact_ranks : bool, optional
        Require exactly as many ranks as tasks in the first group; False
        accepts more ranks, which then stay idle. Default: True.
    max_retries : int, optional
        Number of extra passes over the jobs that did not complete; each pass
        shares them among all ranks. Default: 2.

    Returns
    -------
    rank : int
        This rank, after every rank finished its jobs.
    problems : list of str
        Problems of the jobs still failing after the last pass, on every rank.

    Raises
    ------
    ValueError
        The rank count does not match the prepared group width.
    """
    comm = MPI.COMM_WORLD
    processes_num = comm.Get_size()
    width = len(group_list[0])
    if (processes_num != width) if exact_ranks else (processes_num < width):
        raise ValueError(
            "processes_num is %d, item num in group is %d. \n"
            "Pleasse check the process num!" % (processes_num, width)
        )
    rank = comm.Get_rank()
    # every rank runs one backend process at a time, so a node needs memory
    # for all of its ranks; one rank per node checks before any job starts
    node_comm = comm.Split_type(MPI.COMM_TYPE_SHARED)
    if node_comm.Get_rank() == 0:
        warn_if_memory_short(
            node_comm.Get_size(), memory_per_job_gb, "%s: " % MPI.Get_processor_name()
        )
    node_comm.Free()

    def run(task):
        problems = run_job(task)
        if problems:
            print("rank %d: %s" % (rank, problems[0]), flush=True)
        return problems

    failed = {}
    for ind_group in range(len(group_list)):
        # the last group holds the remainder and may be shorter than processes_num
        if rank >= len(group_list[ind_group]):
            continue
        print("ind_group:%d rank:%d" % (ind_group, rank))
        problems = run(group_list[ind_group][rank])
        if problems:
            failed[group_list[ind_group][rank]] = problems
    for attempt in range(max_retries + 1):
        # allgather waits for every rank and gives all of them the same jobs
        failed_all = {}
        for part in comm.allgather(failed):
            failed_all.update(part)
        if not failed_all or attempt == max_retries:
            break
        if rank == 0:
            print(
                "%d jobs did not complete; computing them again (retry %d/%d)"
                % (len(failed_all), attempt + 1, max_retries),
                flush=True,
            )
        failed = {}
        for ind, task in enumerate(failed_all):
            if ind % processes_num == rank:
                problems = run(task)
                if problems:
                    failed[task] = problems
    return rank, sum(failed_all.values(), [])


def finish_library(path_green, run_problems, check):
    """Check a finished library and raise RuntimeError if it is incomplete.

    Parameters
    ----------
    path_green : str
        Library root.
    run_problems : list of str
        Problems of the jobs still failing after the retries.
    check : callable
        Returns the list of library problems.

    Returns
    -------
    None
    """
    # a failed job's own message says more than the files it did not write
    failed_dirs = [
        problem.split(" failed: ")[0] for problem in run_problems if " failed: " in problem
    ]
    problems = run_problems + [
        problem
        for problem in check()
        if not any(problem.startswith(d + os.sep) for d in failed_dirs)
    ]
    raise_if_incomplete(list(dict.fromkeys(problems)), path_green)
    print("Library check passed: %s" % path_green)


def run_checked_job(job_dir, check_finished, run, check, stale=()):
    """Run one backend job unless check_finished finds it complete.

    Parameters
    ----------
    job_dir : str
        Directory holding the job's .finished marker and .failed log.
    check_finished : bool
        Skip the job when it has a .finished marker and check() finds no problem.
    run : callable
        Runs the executable and returns call_exe's success flag.
    check : callable
        Returns the list of problems with the job's output.
    stale : sequence of str, optional
        Glob patterns, relative to job_dir, of files derived from an earlier run;
        they are removed before the job runs so they cannot hide the new output.

    Returns
    -------
    problems : list of str
        Empty when the job completed.
    """
    try:
        if (
            check_finished
            and os.path.exists(os.path.join(job_dir, ".finished"))
            and not check()
        ):
            return []
        import glob

        for pattern in stale:
            for path in glob.glob(os.path.join(job_dir, pattern)):
                os.remove(path)
        if not run():
            return ["%s failed: %s" % (job_dir, read_failure_log(job_dir))]
        return check()
    except Exception as exc:
        return ["%s failed: %s: %s" % (job_dir, type(exc).__name__, exc)]


# every backend reads and writes paths through character*160 variables
FORTRAN_MAX_PATH = 160


def check_path_lengths(paths, limit=FORTRAN_MAX_PATH):
    """Raise ValueError when a backend would truncate the longest of the paths.

    Parameters
    ----------
    paths : iterable of str
        Paths a backend executable reads or writes.
    limit : int, optional
        Longest path the executable can hold. Default: 160.

    Returns
    -------
    None
    """
    longest = max(paths, key=len)
    if len(longest) > limit:
        raise ValueError(
            "%s has %d characters; the backend executable reads at most %d. "
            "Use a shorter path_green." % (longest, len(longest), limit)
        )


def check_ascii_table(path, n_rows, n_cols, skip_rows=1):
    """Return a problem when a text table is missing or incomplete, else None.

    The table must have n_rows rows after skip_rows header lines, the last one
    with n_cols values; a killed backend leaves fewer or shorter rows.
    """
    if not os.path.exists(path):
        return "%s is missing" % path
    with open(path, "rb") as fr:
        lines = fr.read().splitlines()[skip_rows:]
    while lines and not lines[-1].strip():
        lines.pop()
    if len(lines) != n_rows or len(lines[-1].split()) != n_cols:
        return "%s is incomplete: %d rows of %d values, expected %d rows of %d" % (
            path,
            len(lines),
            len(lines[-1].split()) if lines else 0,
            n_rows,
            n_cols,
        )
    return None


def check_file_size(path, expected, check_values=False):
    """Return a problem when a float32 file is missing or has the wrong size, else None."""
    if not os.path.exists(path):
        return "%s is missing" % path
    size = os.path.getsize(path)
    if size != expected:
        return "%s has %d bytes, expected %d" % (path, size, expected)
    if check_values and not np.all(np.isfinite(np.fromfile(path, dtype=np.float32))):
        return "%s contains NaN or infinite values" % path
    return None


def write_bin_atomic(array, path):
    """Write array.tofile(path) so an interruption never leaves a short file."""
    array.tofile(path + ".tmp")
    os.replace(path + ".tmp", path)


def raise_if_incomplete(problems, path_green, max_listed=20):
    """Raise RuntimeError listing library problems; do nothing when there are none."""
    if not problems:
        return
    lines = problems[:max_listed]
    if len(problems) > max_listed:
        lines.append("... and %d more" % (len(problems) - max_listed))
    raise RuntimeError(
        "The Green's function library at %s is incomplete (%d problems):\n  %s\n"
        "Fix the cause (reduce processes_num if jobs ran out of memory), then rerun "
        "the create_grnlib function with check_finished=True to compute only the "
        "unfinished jobs."
        % (path_green, len(problems), "\n  ".join(lines))
    )


def read_tpts_table(path_green, event_depth_km, receiver_depth_km, ind):
    """Read one pair of flat float32 TauP arrival-table entries.

    Parameters
    ----------
    path_green : str
        Absolute library root containing green_lib_info.json and backend subdirectories.
    event_depth_km : float
        Requested source depth in km, positive down.
    receiver_depth_km : float
        Requested receiver depth in km, positive down.
    ind : int
        Zero-based index into the stored distance sequence.

    Returns
    -------
    first_p, first_s : float

    Raises
    ------
    OSError
        Required inputs or outputs cannot be accessed.
    ValueError
        Parameters do not describe a supported grid or observable.

    Notes
    -----


    """
    fr_tp = open(
        os.path.join(
            path_green,
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "tp_table.bin",
        ),
        "rb",
    )
    tp = np.fromfile(file=fr_tp, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_tp.close()

    fr_ts = open(
        os.path.join(
            path_green,
            "%.2f" % event_depth_km,
            "%.2f" % receiver_depth_km,
            "ts_table.bin",
        ),
        "rb",
    )
    ts = np.fromfile(file=fr_ts, dtype=np.float32, count=1, offset=ind * 4)[0]
    fr_ts.close()
    return float(tp), float(ts)
