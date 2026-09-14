"""Compare five regional displacement examples without fitting shifts or scales.

Run qseis06.py and qseis2025.py with --regional, and run the three spherical
examples first. Input overrides accept an example output directory or disp.npz.
The saved comparison uses origin time, metres, and ENU components throughout.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from source_time_function import validate_qseis_stf


BACKENDS = ("qseis06", "qseis2025", "spgrn2012", "spgrn2020", "qssp2020")
DISTANCES = np.array([300.0, 600.0, 900.0])
COMPONENTS = ("E", "N", "U")
COLORS = ("#9467bd", "#ff7f0e", "#2ca02c", "#1f77b4", "#d62728")
STYLES = ("-", "--", "-.", "-", "--")


def load_result(name, supplied_path):
    """Read a completed example and reject incompatible or stale metadata."""
    path = supplied_path.expanduser().resolve()
    if path.is_dir():
        path = path / "disp.npz"
    candidates = [path.parent / filename
                  for filename in ("summary.json", "summary-reuse.json")]
    candidates = [candidate for candidate in candidates if candidate.is_file()]
    if not candidates:
        raise ValueError("%s: a completed summary.json or summary-reuse.json is required" % name)
    summary_path = max(candidates, key=lambda item: (item.stat().st_mtime_ns, item.name))
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("backend", "").lower() != name:
        raise ValueError("%s: summary backend does not match the input" % name)
    with np.load(path, allow_pickle=False) as archive:
        required = {"values", "time_s", "distance_km", "components", "unit"}
        if not required.issubset(archive.files):
            raise ValueError("%s: missing NPZ fields %s" % (name, sorted(required - set(archive.files))))
        values = archive["values"]
        times = archive["time_s"]
        if (not np.array_equal(archive["distance_km"], DISTANCES)
                or archive["components"].tolist() != list(COMPONENTS)
                or archive["unit"].shape != () or archive["unit"].item() != "m"):
            raise ValueError("%s: expected distances 300/600/900 km, ENU components, and metres" % name)
    if (values.ndim != 3 or values.shape[:2] != (3, 3) or values.shape[2] != 256
            or times.shape != (3, values.shape[2])
            or not np.issubdtype(values.dtype, np.floating)
            or not np.issubdtype(times.dtype, np.floating)
            or not np.isfinite(values).all() or not np.isfinite(times).all()
            or not np.all(np.diff(times, axis=1) > 0)
            or not np.all(np.any(values != 0, axis=(1, 2)))):
        raise ValueError("%s: expected finite nonzero (3, 3, 256) traces and increasing (3, 256) times" % name)
    if summary.get("outputs", {}).get("disp", {}).get("shape") != list(values.shape):
        raise ValueError("%s: summary displacement shape differs from the NPZ" % name)
    expected = {
        "source_depth_km": 10.0, "receiver_depth_km": 0.0,
        "moment_nm": 1e15, "strike_dip_rake_deg": [30.0, 45.0, 90.0],
        "azimuth_deg": 30.0, "distances_km": DISTANCES.tolist(),
        "sampling_interval_s": 4.0,
    }
    if name.startswith("qseis"):
        expected.update(regional=True, flat_earth_transform=True,
                        time_window_s=4092.0, wavelet_type=0,
                        wavelet_duration_samples=16, wavelet_duration_s=64.0,
                        earth_model_numeric_rows=24)
    else:
        expected.update(source_duration_s=64.0, spec_time_window_s=4092.0,
                        max_frequency_hz=0.0625)
    if name == "spgrn2020":
        expected.update(max_slowness_s_km=0.0, full_wavefield=True)
    if name == "qssp2020":
        expected["harmonic_bounds"] = [2000, 8000]
    mismatches = [key for key, value in expected.items() if summary.get(key) != value]
    if mismatches:
        raise ValueError("%s: incompatible or old summary settings: %s; rerun the example"
                         % (name, ", ".join(mismatches)))
    if name.startswith("qseis"):
        stf = summary.get("source_time_function", {})
        if stf.get("scheme_id") != "normalized_sin_squared_64s_damping_compensated_v1":
            raise ValueError("%s: rerun with the shared, damping-compensated source time function" % name)
    if not np.allclose(np.diff(times, axis=1), 4.0, rtol=0.0, atol=1e-10):
        raise ValueError("%s: NPZ time spacing differs from the declared 4 s" % name)
    if name == "spgrn2020":
        starts = summary.get("trace_start_times_s")
        if starts is None or not np.array_equal(times[:, 0], starts):
            raise ValueError("SPGRN2020: use actual native origin-time starts from the updated example")
    elif name == "spgrn2012":
        start_offset, velocity = summary.get("t0_s"), summary.get("v0_km_s")
        if start_offset is None or velocity is None or velocity <= 0:
            raise ValueError("SPGRN2012: missing positive reduction velocity and time offset")
        if not np.array_equal(times[:, 0], start_offset + DISTANCES / velocity):
            raise ValueError("SPGRN2012: NPZ time origin disagrees with reduction settings")
    elif not np.array_equal(times[:, 0], np.zeros(3)):
        raise ValueError("%s: regional examples must start at source origin time zero" % name)
    provenance = {
        "npz_path": str(path), "npz_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "summary_path": str(summary_path), "summary": summary,
        "native_time_ranges_s": times[:, [0, -1]].tolist(),
    }
    library_path = path.parent / "library" / "green_lib_info.json"
    if library_path.is_file():
        library_info = json.loads(library_path.read_text(encoding="utf-8"))
        required_settings = ({"max_slowness": 0.0} if name == "spgrn2020" else
                             {"min_harmonic": 2000, "max_harmonic": 8000}
                             if name == "qssp2020" else {})
        if name.startswith("qseis"):
            required_settings.update(flat_earth_transform=True, time_window=4092.0,
                                     sampling_interval=4.0, wavelet_duration=16,
                                     wavelet_type=0, earth_model_layer_num=24)
        if any(library_info.get(key) != value for key, value in required_settings.items()):
            raise ValueError("%s: library metadata contains old comparison settings" % name)
        if name.startswith("qseis"):
            validate_qseis_stf(library_path.parent)
        provenance.update(library_metadata_path=str(library_path), library_metadata=library_info)
    return {"values": values.astype(float), "times": times.astype(float),
            "provenance": provenance}


def metrics(values, reference, times):
    """Return unscaled L2 error, Pearson correlation, and absolute peak data."""
    denominator = float(np.linalg.norm(reference))
    centred = values.ravel() - float(np.mean(values))
    reference_centred = reference.ravel() - float(np.mean(reference))
    correlation_denominator = float(np.linalg.norm(centred) * np.linalg.norm(reference_centred))
    peak_index = int(np.argmax(np.abs(values)))
    return {
        "relative_l2": float(np.linalg.norm(values - reference) / denominator) if denominator else None,
        "correlation": float(np.clip(np.dot(centred, reference_centred) / correlation_denominator,
                                     -1.0, 1.0)) if correlation_denominator else None,
        "peak_absolute_m": float(np.max(np.abs(values))),
        "reference_peak_absolute_m": float(np.max(np.abs(reference))),
        "absolute_peak_time_s": float(times[peak_index % times.size]),
    }


def plot_comparison(output_path, selected, grids, interpolated, title):
    """Plot components in rows and distances in columns on physical axes."""
    figure, axes = plt.subplots(3, 3, figsize=(12, 8), sharex=True, constrained_layout=True)
    for distance_index, distance in enumerate(DISTANCES):
        for component_index, component in enumerate(COMPONENTS):
            axis = axes[component_index, distance_index]
            for name in selected:
                style_index = BACKENDS.index(name)
                axis.plot(grids[distance_index], interpolated[name][distance_index][component_index],
                          color=COLORS[style_index], ls=STYLES[style_index], lw=1.15,
                          label=name.upper())
            axis.set_xlim(0, 500)
            axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
            axis.grid(alpha=0.2)
            axis.set_ylabel("%s displacement (m)" % component)
            if component_index == 0:
                axis.set_title("%g km" % distance)
            if component_index == 2:
                axis.set_xlabel("Time since origin (s)")
    axes[0, 0].legend(fontsize=7, loc="best")
    figure.suptitle(title + "\nM0 = 10^15 N m; no fitted time shifts or amplitude factors")
    figure.savefig(output_path, dpi=150)
    plt.close(figure)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    default_root = Path(__file__).resolve().parent / "output"
    for name in BACKENDS:
        folder = name + "-regional" if name.startswith("qseis") else name
        parser.add_argument("--" + name, type=Path, default=default_root / folder,
                            help="%s output directory or disp.npz" % name.upper())
    parser.add_argument("--output-dir", type=Path, default=default_root / "backend-comparison")
    args = parser.parse_args()
    results = {name: load_result(name, getattr(args, name)) for name in BACKENDS}
    grids = []
    for index, distance in enumerate(DISTANCES):
        start = math.ceil(max(result["times"][index, 0] for result in results.values()))
        end = min(result["times"][index, -1] for result in results.values())
        if start >= 500 or end < 500:
            raise ValueError("%g km: all inputs must cover a common interval ending at 500 s" % distance)
        grids.append(np.arange(start, 501.0, 1.0))
    interpolated = {
        name: [np.array([np.interp(grid, result["times"][index], component)
                         for component in result["values"][index]])
               for index, grid in enumerate(grids)]
        for name, result in results.items()
    }
    report = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "reference": "SPGRN2020 full wavefield (max_slowness=0)",
        "method": "Linear interpolation to 1 s on all five backends' common valid origin-time interval at each distance through 500 s",
        "fitted_time_shifts": False, "fitted_amplitude_factors": False,
        "relative_l2_definition": "norm(values-reference)/norm(reference), over ENU and time for each distance",
        "correlation_definition": "Pearson correlation; null when either demeaned trace has zero norm",
        "components": list(COMPONENTS), "unit": "m",
        "comparison_time_ranges_s": [grid[[0, -1]].tolist() for grid in grids],
        "provenance": {name: result["provenance"] for name, result in results.items()},
        "limitations": [
            "QSEIS uses a transformed finite-depth model over a half-space; spherical examples use the complete Earth model.",
            "QSEIS uses a damping-compensated custom 64 s sin-squared rate matching the SPGRN2020/QSSP physical source to the recorded interpolation tolerance.",
            "SPGRN2012 retains its legacy real-frequency source; displacement integration windows still differ between backends.",
            "QSEIS Nyquist frequency is 0.125 Hz; spherical examples use a 0.0625 Hz spectral cutoff.",
            "Displacement integration baselines and native time windows differ between backends.",
            "Thresholds check these tutorial settings; they do not establish general backend equivalence.",
        ],
        "distances": [],
    }
    for index, distance in enumerate(DISTANCES):
        reference = interpolated["spgrn2020"][index]
        if not np.linalg.norm(reference):
            raise ValueError("SPGRN2020 reference is zero over the %g km comparison window" % distance)
        comparisons = {}
        for name in BACKENDS:
            values = interpolated[name][index]
            comparison = metrics(values, reference, grids[index])
            comparison["components"] = {
                component: metrics(values[j], reference[j], grids[index])
                for j, component in enumerate(COMPONENTS)
            }
            comparisons[name] = comparison
        report["distances"].append({"distance_km": float(distance), "backends": comparisons})
    first, second = results["qseis06"], results["qseis2025"]
    matching_axes = np.array_equal(first["times"], second["times"])
    matching_shape = first["values"].shape == second["values"].shape
    qseis_pass = bool(matching_axes and matching_shape and np.allclose(
        first["values"], second["values"], rtol=1e-5, atol=1e-20))
    qssp_errors = [item["backends"]["qssp2020"]["relative_l2"] for item in report["distances"]]
    checks = {
        "qseis06_vs_qseis2025": {
            "passed": qseis_pass, "same_native_time_axes": matching_axes,
            "rtol": 1e-5, "atol_m": 1e-20, "scope": "all saved native displacement samples",
            "max_absolute_difference_m": float(np.max(np.abs(first["values"] - second["values"])))
            if matching_shape else None,
        },
        "qssp2020_vs_spgrn2020": {
            "passed": all(error <= 0.05 for error in qssp_errors),
            "relative_l2_limit": 0.05, "relative_l2_by_distance": qssp_errors,
        },
    }
    report["checks"] = checks
    report["passed"] = all(check["passed"] for check in checks.values())
    output = args.output_dir.expanduser().resolve()
    output.mkdir(parents=True, exist_ok=True)
    groups = (
        ("all-backends.png", BACKENDS, "Five backend displacement comparison"),
        ("spherical-comparison.png", BACKENDS[2:], "Spherical backend displacement comparison"),
        ("qseis-comparison.png", BACKENDS[:2], "QSEIS06 and QSEIS2025 displacement comparison"),
    )
    for filename, selected, title in groups:
        plot_comparison(output / filename, selected, grids, interpolated, title)
    report["figures"] = [filename for filename, _, _ in groups]
    (output / "comparison.json").write_text(json.dumps(report, indent=2, allow_nan=False) + "\n",
                                             encoding="utf-8")
    print(json.dumps({"output_dir": str(output), "checks": checks, "passed": report["passed"]}, indent=2))
    if not report["passed"]:
        raise SystemExit("Backend comparison failed; see comparison.json and figures for diagnostics")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as error:
        print("Comparison input error: %s" % error, file=sys.stderr)
        raise SystemExit(1)
