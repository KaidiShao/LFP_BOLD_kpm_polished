function params = build_blp_eigenfunction_reduction_params()
%BUILD_BLP_EIGENFUNCTION_REDUCTION_PARAMS Default batch/script params for pipeline 5.

params = struct();
params.autodl_root = 'E:\autodl_results';
params.processed_root = io_project.get_project_processed_root();

params.dataset_stems = {};
params.exclude_dataset_stems = {};
params.filename_pattern = '*_outputs_*.mat';
params.variable_name = 'EDMD_outputs';

params.condition_key_filter = {};
params.run_name_filter = {};
params.condition_output_tag = '';
params.condition_tag_mode = 'condition';
params.method_filter = {};
params.component_count_sweep = 4;
params.time_component_count_sweep = [];
params.spectrum_component_count_sweep = [];

params.continue_on_error = true;
params.close_figures = true;
params.write_manifest = true;

params.make_overview_plot = false;
params.make_state_space_plot = true;
params.make_consensus_state_space_plot = true;
params.make_spectrum_diagnostics = true;
params.make_top30_window_plots = true;
params.top30_n_windows = 30;
params.top30_force_recompute_norm_cache = false;

params.do_thresholded_density = true;
params.do_thresholded_events = false;
params.do_dimred_thresholded_density = false;
params.do_dimred_thresholded_events = false;
params.threshold_mode = 'quantile';
params.threshold_ratio = 0.7;
params.threshold_ratio_sweep = [];
params.save_thresholded_density_results = [];
params.save_dimred_thresholded_density_results = [];

params.feature_variant = 'abs';
params.feature_normalization = 'maxabs_per_mode';
params.save_payload = 'compact';
params.save_v7_3 = true;
end
