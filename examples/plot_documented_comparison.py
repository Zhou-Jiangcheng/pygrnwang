"""Redraw the documented 2026-09-15 comparison from published NumPy arrays.

Requires Python >= 3.9, NumPy and Matplotlib. No pygrnwang installation is needed.
The input archives contain both raw and already filtered common-grid waveforms;
this script selects those arrays without running solvers or filtering again.

From a repository checkout::

    python examples/plot_documented_comparison.py --output-dir comparison-plots

For separately downloaded archives, add ``--data-dir PATH``. That directory must
contain ``point-source-waves.npz`` and ``baseline-waves.npz``. The QSEIS06 curve
uses the isolated point-source control, not its stock Gaussian source setting.

Eight PNG files are written: full and zoomed three-component displacement for
all distances, and full and zoomed six-component stress for each distance.
Existing output files are refused unless ``--overwrite`` is supplied.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter
import numpy as np

BACKENDS = ("qseis06", "qseis2025", "spgrn2012", "spgrn2020", "qssp2020")
STRESS_BACKENDS = ("qseis2025", "qssp2020")
DISPLACEMENT_COMPONENTS = ("East", "North", "Up")
STRESS_COMPONENTS = ("EE", "EN", "EU", "NN", "NU", "UU")
DISTANCES_KM = np.array([300.0, 600.0, 900.0])
TIME_S = np.arange(1601, dtype=float) * 0.25
ZOOM_WINDOWS_S = ((75.0, 115.0), (155.0, 220.0), (240.0, 310.0))
LABELS = {
    "qseis06": "QSEIS06 (deprecated)\npoint-source control, ratio 0",
    "qseis2025": "QSEIS2025\npoint source",
    "spgrn2012": "SPGRN2012 (deprecated)",
    "spgrn2020": "SPGRN2020",
    "qssp2020": "QSSP2020",
}
STYLES = {
    "qseis06": ("#d55e00", (0, (6, 3)), 1.6),
    "qseis2025": ("#d62728", "-", 1.0),
    "spgrn2012": ("#009e73", (0, (2, 2)), 1.5),
    "spgrn2020": ("#252525", (0, (8, 2, 2, 2)), 1.2),
    "qssp2020": ("#1759d1", "-", 0.8),
}


def read_array(archive, key: str, shape: tuple, filename: str) -> np.ndarray:
    """Load a real finite array and reject an unexpected archive structure."""
    if key not in archive:
        raise ValueError(f"{filename}: required array {key!r} is missing")
    values = archive[key]
    if values.shape != shape or values.dtype.kind not in "fiu":
        raise ValueError(f"{filename}: {key} must be real with shape {shape}; "
                         f"got {values.shape}, {values.dtype}")
    if not np.isfinite(values).all():
        raise ValueError(f"{filename}: {key} contains non-finite values")
    return np.asarray(values, dtype=float)


def check_coordinates(archive, filename: str) -> None:
    """Require the documented 300/600/900 km and 0:0.25:400 s grids."""
    times = read_array(archive, "time_s", (1601,), filename)
    distances = read_array(archive, "distances_km", (3,), filename)
    if not np.allclose(times, TIME_S, rtol=0, atol=1e-10):
        raise ValueError(f"{filename}: expected time_s = 0..400 s with dt=0.25 s")
    if not np.allclose(distances, DISTANCES_KM, rtol=0, atol=1e-8):
        raise ValueError(f"{filename}: expected distances_km = [300, 600, 900]")


def load_waveforms(data_dir: Path, processing: str) -> tuple:
    """Read displacement from the point control and stress from the baseline."""
    point_path = data_dir / "point-source-waves.npz"
    baseline_path = data_dir / "baseline-waves.npz"
    with np.load(point_path, allow_pickle=False) as archive:
        check_coordinates(archive, point_path.name)
        displacement = {
            backend: read_array(archive, f"{processing}_{backend}_disp",
                                (3, 3, 1601), point_path.name)
            for backend in BACKENDS
        }
    with np.load(baseline_path, allow_pickle=False) as archive:
        check_coordinates(archive, baseline_path.name)
        stress = {
            backend: read_array(archive, f"{processing}_{backend}_stress",
                                (3, 6, 1601), baseline_path.name)
            for backend in STRESS_BACKENDS
        }
        q25_baseline = read_array(archive, f"{processing}_qseis2025_disp",
                                  (3, 3, 1601), baseline_path.name)
    if not np.array_equal(displacement["qseis2025"], q25_baseline):
        raise ValueError("The archives disagree on QSEIS2025 displacement; "
                         "use both files from the same documented comparison")
    return displacement, stress


def style_axis(axis, ylabel: str, window: tuple) -> None:
    axis.set_xlim(*window)
    axis.set_xlabel("Time since earthquake origin (s)", fontsize=9)
    axis.set_ylabel(ylabel, fontsize=9)
    axis.axhline(0, color="0.82", linewidth=0.5, zorder=0)
    axis.grid(axis="x", color="0.9", linewidth=0.5)
    formatter = ScalarFormatter(useMathText=True)
    formatter.set_powerlimits((-2, 3))
    axis.yaxis.set_major_formatter(formatter)
    axis.tick_params(direction="in", top=True, right=True, labelsize=8)
    axis.margins(y=0.09)


def mode_label(processing: str) -> str:
    return "Raw archived waveforms" if processing == "raw" else "Archived 0.4 Hz low-pass"


def finish_plot(figure, axes, title: str, path: Path, overwrite: bool,
                legend_columns: int, footer: str) -> None:
    handles, labels = np.asarray(axes).flat[0].get_legend_handles_labels()
    figure.suptitle(title, fontsize=14, y=0.992)
    figure.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.961),
                  frameon=False, ncol=legend_columns, fontsize=9)
    figure.tight_layout(rect=(0, 0.05, 1, 0.915), h_pad=2, w_pad=2)
    figure.text(0.5, 0.015, footer, ha="center", va="bottom", fontsize=8, color="0.3")
    try:
        # Exclusive creation also prevents an accidental overwrite after preflight.
        with path.open("wb" if overwrite else "xb") as stream:
            figure.savefig(stream, format="png", dpi=160, facecolor="white")
    finally:
        plt.close(figure)


def plot_displacement(values: dict, processing: str, zoom: bool,
                      path: Path, overwrite: bool) -> None:
    figure, axes = plt.subplots(3, 3, figsize=(15, 9), sharex="col")
    for distance_index, distance in enumerate(DISTANCES_KM):
        window = ZOOM_WINDOWS_S[distance_index] if zoom else (0.0, 400.0)
        selected = (TIME_S >= window[0]) & (TIME_S <= window[1])
        for component_index, component in enumerate(DISPLACEMENT_COMPONENTS):
            axis = axes[component_index, distance_index]
            for backend in BACKENDS:
                color, line_style, width = STYLES[backend]
                axis.plot(TIME_S[selected],
                          values[backend][distance_index, component_index, selected] * 1e6,
                          color=color, linestyle=line_style, linewidth=width,
                          label=LABELS[backend])
            axis.set_title(f"{distance:g} km / {component}", fontsize=10)
            style_axis(axis, f"{component} displacement ($\\mu$m)", window)
    view = "zoom" if zoom else "0-400 s"
    finish_plot(figure, axes,
                f"Point-source displacement comparison | {mode_label(processing)} | {view}",
                path, overwrite, 5,
                "Physical displacement in micrometres; no trace normalization or fitted time shifts.\n"
                "QSEIS06 uses the isolated point-source control (stock source-radius ratio: 0.05).")


def plot_stress(values: dict, processing: str, distance_index: int, zoom: bool,
                path: Path, overwrite: bool) -> None:
    figure, axes = plt.subplots(3, 2, figsize=(12, 9), sharex=True)
    distance = DISTANCES_KM[distance_index]
    window = ZOOM_WINDOWS_S[distance_index] if zoom else (0.0, 400.0)
    selected = (TIME_S >= window[0]) & (TIME_S <= window[1])
    for component_index, component in enumerate(STRESS_COMPONENTS):
        axis = axes.flat[component_index]
        for backend in STRESS_BACKENDS:
            color, line_style, width = STYLES[backend]
            axis.plot(TIME_S[selected], values[backend][distance_index, component_index, selected],
                      color=color, linestyle=line_style, linewidth=width,
                      label=LABELS[backend])
        axis.set_title(f"{distance:g} km / {component}", fontsize=10)
        style_axis(axis, f"Stress {component} (Pa)", window)
    view = "zoom" if zoom else "0-400 s"
    finish_plot(figure, axes,
                f"Point-source stress comparison | {distance:g} km | {mode_label(processing)} | {view}",
                path, overwrite, 2,
                "ENU tensor components: EE, EN, EU, NN, NU, UU; tensile stress is positive (Pa).\n"
                "Archived amplitudes and origin times; no trace normalization or fitted time shifts.")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--data-dir", type=Path,
                        default=Path(__file__).resolve().parents[1] / "docs" / "_static" /
                        "comparisons" / "2026-09-15",
                        help="Directory containing the two documented NPZ archives")
    parser.add_argument("--output-dir", type=Path, required=True,
                        help="Explicit directory for the eight generated PNG files")
    parser.add_argument("--processing", choices=("raw", "lowpass_0p4Hz"),
                        default="lowpass_0p4Hz", help="Select archived arrays; no new filtering")
    parser.add_argument("--overwrite", action="store_true",
                        help="Allow replacing the eight named PNG files")
    args = parser.parse_args()
    try:
        data_dir = args.data_dir.expanduser().resolve()
        output_dir = args.output_dir.expanduser().resolve()
        displacement, stress = load_waveforms(data_dir, args.processing)
        jobs = [("displacement", None, zoom,
                 f"displacement-{'zoom' if zoom else 'full'}-{args.processing}.png")
                for zoom in (False, True)]
        jobs += [("stress", index, zoom,
                  f"stress-{distance:g}km-{'zoom' if zoom else 'full'}-{args.processing}.png")
                 for index, distance in enumerate(DISTANCES_KM) for zoom in (False, True)]
        paths = [output_dir / name for _, _, _, name in jobs]
        for path in paths:
            if path.is_symlink() or path.resolve().parent != output_dir:
                raise ValueError(f"Output must be an ordinary file within {output_dir}: {path}")
            if path.exists() and (not args.overwrite or not path.is_file()):
                raise ValueError(f"Refusing to overwrite {path}; use --overwrite for existing PNG files")
        output_dir.mkdir(parents=True, exist_ok=True)
        for (quantity, distance_index, zoom, _), path in zip(jobs, paths):
            if quantity == "displacement":
                plot_displacement(displacement, args.processing, zoom, path, args.overwrite)
            else:
                plot_stress(stress, args.processing, distance_index, zoom, path, args.overwrite)
            print(path)
    except (OSError, ValueError, KeyError, EOFError) as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
