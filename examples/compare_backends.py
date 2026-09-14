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

from common import REGIONAL_MAX_FREQUENCY_HZ, REGIONAL_STF
from source_time_function import validate_qseis_stf
from spectral_settings import qseis_spectral_settings, verify_spherical_spectrum
from spherical_source_time_function import validate_spgrn2012_stf


BACKENDS = ("qseis06", "qseis2025", "spgrn2012", "spgrn2020", "qssp2020")
DISTANCES = np.array([300.0, 600.0, 900.0])
COMPONENTS = ("E", "N", "U")
COLORS = ("#9467bd", "#ff7f0e", "#2ca02c", "#1f77b4", "#d62728")
STYLES = ("-", "--", "-.", "-", "--")


def load_result(name, supplied_path, source_radius_ratio=0.05):
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
        "sampling_interval_s": 4.0, "max_frequency_hz": REGIONAL_MAX_FREQUENCY_HZ,
        "physical_source_time_function": REGIONAL_STF,
    }
    if name.startswith("qseis"):
        expected.update(regional=True, flat_earth_transform=True,
                        time_window_s=4092.0, wavelet_type=0,
                        wavelet_duration_samples=16, wavelet_duration_s=64.0,
                        earth_model_numeric_rows=24, source_radius_ratio=source_radius_ratio)
    else:
        expected.update(source_duration_s=64.0, spec_time_window_s=4092.0)
    if name in ("spgrn2012", "spgrn2020"):
        expected.update(max_slowness_s_km=0.0, full_wavefield=True)
    if name == "spgrn2012":
        expected.update(native_source_duration_s=0.0, effective_source_duration_s=64.0,
                        time_window_s=4092.0, native_samples=1024, output_window_s=1020.0)
    elif name in ("spgrn2020", "qssp2020"):
        expected["time_window_s"] = 1020.0
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
        starts = summary.get("trace_start_times_s")
        if starts is None or not np.array_equal(times[:, 0], starts):
            raise ValueError("SPGRN2012: use the actual native origin-time starts")
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
    library = library_path.parent
    library_info = json.loads(library_path.read_text(encoding="utf-8"))
    required_settings = {"sampling_interval": 4.0, "anti_alias": 0.01}
    if name.startswith("qseis"):
        required_settings.update(flat_earth_transform=True, time_window=4092.0,
                                 sampling_num=1024, wavelet_duration=16, wavelet_type=0,
                                 earth_model_layer_num=24, time_reduction_velo=0,
                                 slowness_window=None, wavenumber_sampling_rate=12,
                                 free_surface=True if name == "qseis06" else 0)
    else:
        required_settings.update(max_frequency=REGIONAL_MAX_FREQUENCY_HZ,
                                 spec_time_window=4092.0, physical_dispersion=0,
                                 source_radius=0.0, gravity_fc=0.0, gravity_harmonic=0,
                                 cal_sph=1, cal_tor=1)
    if name == "spgrn2012":
        required_settings.update(source_duration=0.0, time_window=4092.0,
                                 samples_num=1024, max_slowness=0.0)
    elif name == "spgrn2020":
        required_settings.update(source_duration=64.0, time_window=1020.0, max_slowness=0.0)
    elif name == "qssp2020":
        required_settings.update(source_duration=64.0, time_window=1020.0,
                                 min_harmonic=2000, max_harmonic=8000, max_slowness=0.3)
    if any(library_info.get(key) != value for key, value in required_settings.items()):
        raise ValueError("%s: library metadata contains old comparison settings" % name)
    if name.startswith("qseis"):
        source = validate_qseis_stf(library)
        spectrum = qseis_spectral_settings(library)
        if name == "qseis2025":
            input_path = library / "10.00" / "0.00" / "0_0" / "grn.inp"
            with input_path.open(encoding="utf-8") as stream:
                records = []
                for line in stream:
                    if line.strip() and not line.lstrip().startswith("#"):
                        records.append(line.split())
                    if len(records) == 9:
                        break
            if len(records) != 9 or len(records[8]) != 2:
                raise ValueError("QSEIS2025: missing native wavenumber/source-radius record")
            eps, native_ratio = (float(value.replace("D", "e").replace("d", "e"))
                                 for value in records[8])
            if eps != 1e-6 or native_ratio != source_radius_ratio:
                raise ValueError("QSEIS2025: native source radius or wavenumber tolerance differs")
        elif source_radius_ratio != 0.05:
            raise ValueError("QSEIS06 has a fixed native source_radius_ratio of 0.05")
        if summary.get("source_time_function") != source:
            raise ValueError("%s: summary STF differs from the validated source evidence" % name)
    else:
        spectrum = verify_spherical_spectrum(library, name)
        if name == "spgrn2012":
            source = validate_spgrn2012_stf(library)
            if summary.get("source_time_function") != source:
                raise ValueError("SPGRN2012: summary STF differs from the validated source evidence")
    if summary.get("spectral_settings") != spectrum:
        raise ValueError("%s: summary spectral_settings differs from the verified native grid; "
                         "rerun the example to record current evidence" % name)
    provenance.update(library_metadata_path=str(library_path), library_metadata=library_info,
                      verified_spectral_settings=spectrum)
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
                if name == "qseis2025_point_source":
                    color, style, label = "#222222", ":", "QSEIS2025 point source"
                else:
                    style_index = BACKENDS.index(name)
                    color, style, label = COLORS[style_index], STYLES[style_index], name.upper()
                axis.plot(grids[distance_index], interpolated[name][distance_index][component_index],
                          color=color, ls=style, lw=1.15, label=label)
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
    parser.add_argument("--qseis2025-point-source", type=Path,
                        help="Optional QSEIS2025 source_radius_ratio=0 control output or disp.npz")
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
    control = None
    if args.qseis2025_point_source is not None:
        control = load_result("qseis2025", args.qseis2025_point_source, source_radius_ratio=0.0)
        for index, grid in enumerate(grids):
            if control["times"][index, 0] > grid[0] or control["times"][index, -1] < grid[-1]:
                raise ValueError("Point-source control must cover the five-backend comparison grid")
        interpolated["qseis2025_point_source"] = [
            np.array([np.interp(grid, control["times"][index], component)
                      for component in control["values"][index]])
            for index, grid in enumerate(grids)]
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
            "SPGRN2012 uses a forward-convolved full-period impulse response with the same physical 64 s moment-rate pulse.",
            "All backends use dt=4 s and the 0.125 Hz Nyquist band; the highest computed frequency is 0.124755859375 Hz and the Nyquist bin is zero.",
            "The standard QSEIS pair uses source_radius_ratio=0.05 spatial smoothing; the spherical examples use point sources.",
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
    if control is not None:
        report["source_radius_control"] = {
            "backend": "QSEIS2025", "source_radius_ratio": 0.0,
            "default_source_radius_ratio": 0.05, "threshold_applied": False,
            "method": "Same physical STF, frequency band, origin-time grids and unscaled metrics",
            "provenance": control["provenance"],
            "distances": [
                {"distance_km": float(distance),
                 "vs_spgrn2020": metrics(interpolated["qseis2025_point_source"][index],
                                        interpolated["spgrn2020"][index], grids[index]),
                 "vs_qseis2025_default": metrics(interpolated["qseis2025_point_source"][index],
                                                interpolated["qseis2025"][index], grids[index])}
                for index, distance in enumerate(DISTANCES)],
        }
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
    if control is not None:
        filename = "source-radius-comparison.png"
        plot_comparison(output / filename,
                        ("qseis2025", "qseis2025_point_source", "spgrn2020"),
                        grids, interpolated, "QSEIS2025 source-radius control against SPGRN2020")
        report["figures"].append(filename)
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
