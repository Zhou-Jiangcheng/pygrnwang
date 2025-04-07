from .create_edgrn_bulk import (
    pre_process_edgrn2,
    create_grnlib_edgrn2_parallel_single_node,
)
from .create_edcmp_bulk import (
    pre_process_edcmp2,
    compute_static_stress_edcmp2_parallel_single_node,
)


def pre_process_edgrn2_edcmp2(
    processes_num,
    path_green,
    path_bin_edgrn,
    path_bin_edcmp,
    grn_source_depth_range,
    grn_source_delta_depth,
    grn_dist_range,
    grn_delta_dist,
    obs_depth_list,
    obs_x_range,
    obs_y_range,
    obs_delta_x,
    obs_delta_y,
    source_array,
    source_ref=None,
    obs_ref=None,
    wavenumber_sampling_rate=12,
    path_nd=None,
    earth_model_layer_num=None,
    layered=True,
    lam=30516224000,
    mu=33701888000,
):
    pre_process_edgrn2(
        processes_num,
        path_green,
        path_bin_edgrn,
        grn_source_depth_range,
        grn_source_delta_depth,
        grn_dist_range,
        grn_delta_dist,
        obs_depth_list,
        wavenumber_sampling_rate,
        path_nd,
        earth_model_layer_num,
    )
    pre_process_edcmp2(
        processes_num,
        path_green,
        path_bin_edcmp,
        obs_depth_list,
        obs_x_range,
        obs_y_range,
        obs_delta_x,
        obs_delta_y,
        source_array,
        source_ref,
        obs_ref,
        layered,
        lam,
        mu,
    )


def create_grnlib_edgrn2_edcmp2_parallel_single_node(path_green, check_finished=False):
    run_time = create_grnlib_edgrn2_parallel_single_node(path_green, check_finished)
    if run_time:
        compute_static_stress_edcmp2_parallel_single_node(path_green, check_finished)
