"""Shared paths, model preparation, and output checks for executable tutorials."""
import argparse
import importlib.metadata
import json
from pathlib import Path
import platform
import sys
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pygrnwang.ak135fc import s as AK135_ELASTIC_MODEL


MOMENT_NM = 1e15
MECHANISM = [30.0, 45.0, 90.0]


def parser_for(name):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path,
                        default=Path(__file__).resolve().parent / "output" / name,
                        help="Directory for the model, library, figures and summary.")
    parser.add_argument("--reuse", action="store_true",
                        help="Read the existing library without rerunning a solver.")
    return parser


def prepare(args, backend, extra=None):
    """Resolve paths before solvers change cwd, and write a six-column ND model."""
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    model = output / "ak135_tutorial.nd"
    # AK135 elastic velocities/density are bundled; constant Q is an explicit
    # tutorial choice, not the frequency-dependent AK135-F attenuation model.
    model_text = "\n".join(line + "  600.0  300.0" if len(line.split()) > 1 else line
                           for line in AK135_ELASTIC_MODEL.splitlines()) + "\n"
    if args.reuse:
        if not model.is_file() or model.read_text(encoding="utf-8") != model_text:
            raise ValueError("--reuse requires the model from a completed tutorial run")
    else:
        model.write_text(model_text, encoding="utf-8")
    library = output / "library"
    library.mkdir(exist_ok=True)
    report = {"backend": backend, "started_utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
              "platform": platform.platform(), "python": sys.version.split()[0],
              "dependencies": {name: importlib.metadata.version(name)
                               for name in ("numpy", "scipy", "pandas", "matplotlib", "obspy")},
              "source_depth_km": 10.0, "receiver_depth_km": 0.0,
              "moment_nm": MOMENT_NM, "strike_dip_rake_deg": MECHANISM,
              "azimuth_deg": 30.0, "reuse": args.reuse,
              "model": "Bundled AK135 elastic properties; illustrative constant Qp=600, Qs=300",
              "outputs": {}}
    report.update(extra or {})
    return output, str(library), str(model), report, time.perf_counter()


def save_waveforms(output, report, name, arrays, distances, dt, labels,
                   unit, time_label="Time since origin (s)", start_times=None):
    """Check shape/finiteness/nonzero output and save physical-unit arrays/plots."""
    values = np.asarray(arrays)
    if values.ndim != 3 or values.shape[:2] != (len(distances), len(labels)):
        raise AssertionError("Unexpected waveform shape: %s" % (values.shape,))
    if values.shape[2] != 256 or not np.isfinite(values).all() or not np.all(np.any(values != 0, axis=(1, 2))):
        raise AssertionError("Each distance must have 256 finite samples and a nonzero waveform")
    starts = np.zeros(len(distances)) if start_times is None else np.asarray(start_times)
    times = starts[:, None] + np.arange(values.shape[2])[None, :] * dt
    np.savez_compressed(output / (name + ".npz"), values=values, time_s=times,
                        distance_km=distances, components=labels, unit=unit)
    fig, axes = plt.subplots(len(labels), 1, figsize=(8, max(4, len(labels) * 1.45)),
                             sharex=True, constrained_layout=True, squeeze=False)
    for component, label in enumerate(labels):
        axis = axes[component, 0]
        for index, distance in enumerate(distances):
            axis.plot(times[index], values[index, component], lw=1.1,
                      label="%g km" % distance)
        axis.set_ylabel("%s (%s)" % (label, unit))
        axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
        axis.grid(alpha=0.2)
    axes[0, 0].legend(ncol=len(distances), fontsize=8)
    axes[0, 0].set_title("%s: %s, M0 = 10^15 N m" % (report["backend"], name))
    axes[-1, 0].set_xlabel(time_label)
    fig.savefig(output / (name + ".png"), dpi=140)
    plt.close(fig)
    report["outputs"][name] = {"shape": list(values.shape), "components": labels,
                              "unit": unit, "finite": True,
                              "peak_absolute": float(np.max(np.abs(values)))}


def finish(output, report, started):
    report["elapsed_seconds"] = round(time.perf_counter() - started, 3)
    files = [p for p in output.rglob("*") if p.is_file() and not p.name.startswith("summary")]
    report["file_count"] = len(files)
    report["output_bytes"] = sum(p.stat().st_size for p in files)
    summary_name = "summary-reuse.json" if report["reuse"] else "summary.json"
    (output / summary_name).write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2))
    print("Results: %s" % output)
