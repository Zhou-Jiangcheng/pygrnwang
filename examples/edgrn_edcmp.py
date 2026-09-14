"""Build layered EDGRN tables, use EDCMP, and plot static displacement."""
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

from common import MECHANISM, MOMENT_NM, finish, parser_for, prepare
from pygrnwang.create_edgrn_bulk import pre_process_edgrn2, create_grnlib_edgrn2_sequential
from pygrnwang.create_edcmp_bulk import (
    pre_process_edcmp2, create_grnlib_edcmp2_sequential, convert_pd2bin_edcmp2_all)
from pygrnwang.read_edcmp import seek_edcmp2


def main():
    args = parser_for("edgrn_edcmp").parse_args()
    output, library, model, report, started = prepare(args, "EDGRN2 + EDCMP2")
    # Keep queries inside the table: finite-source corners can cross its edges.
    distances = [30.0, 60.0, 90.0]
    # EDGRN requires at least two source depths, even for one queried depth.
    grid = dict(processes_num=1, path_green=library,
                grn_source_depth_range=[10.0, 11.0], grn_source_delta_depth=1.0,
                grn_dist_range=[0.0, 120.0], grn_delta_dist=30.0,
                obs_depth_list=[0.0])
    if not args.reuse:
        pre_process_edgrn2(**grid, path_nd=model, earth_model_layer_num=24,
                          wavenumber_sampling_rate=12)
        create_grnlib_edgrn2_sequential(library)
        # Same geometry and root: EDCMP consumes the EDGRN tables from step one.
        pre_process_edcmp2(**grid, output_observables=(1, 0, 0, 0), layered=True)
        create_grnlib_edcmp2_sequential(library)
        # Serial EDCMP does not perform the bulk conversion automatically.
        convert_pd2bin_edcmp2_all(library, remove=False)
    # EDCMP normalizes the mechanism to M0=1. Scale its result explicitly.
    values = np.asarray([MOMENT_NM * seek_edcmp2(
        path_green=library, event_depth_km=10.0, receiver_depth_km=0.0,
        az_deg=30.0, dist_km=distance, focal_mechanism=MECHANISM,
        rotate=True, output_type="disp", times_mu=False, model_name=str(Path(library) / "noQ.nd"),
    ) for distance in distances])
    if values.shape != (3, 3) or not np.isfinite(values).all() or not np.any(values):
        raise AssertionError("Invalid static displacement output")
    labels = ["E", "N", "U"]
    np.savez_compressed(output / "disp.npz", values=values, distance_km=distances,
                        components=labels, unit="m")
    fig, axis = plt.subplots(figsize=(8, 4), constrained_layout=True)
    for index, label in enumerate(labels):
        axis.plot(distances, values[:, index], "o-", label=label)
    axis.set(xlabel="Epicentral distance (km)", ylabel="Static displacement (m)",
             title="EDGRN2 + EDCMP2, M0 = 10^15 N m")
    axis.ticklabel_format(axis="y", style="sci", scilimits=(-2, 2))
    axis.legend()
    axis.grid(alpha=0.2)
    fig.savefig(output / "disp.png", dpi=140)
    plt.close(fig)
    report["outputs"]["disp"] = {"shape": list(values.shape), "components": labels,
                                  "unit": "m", "finite": True,
                                  "peak_absolute": float(np.max(np.abs(values)))}
    report.update(distances_km=distances, earth_model_numeric_rows=24, source_grid_km=[10.0, 11.0], distance_grid_km=[0.0, 30.0, 60.0, 90.0, 120.0])
    finish(output, report, started)


if __name__ == "__main__":
    main()
