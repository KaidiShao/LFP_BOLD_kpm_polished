function params = build_bold_dimred_cross_modal_coupling_params()
%BUILD_BOLD_DIMRED_CROSS_MODAL_COUPLING_PARAMS Default parameters for pipeline 10.

base = build_bold_cross_modal_coupling_params();

params = struct();
params.processed_root = base.processed_root;
params.datapons_root = base.datapons_root;

params.dataset_stems = {};
params.exclude_dataset_stems = {};
params.observable_modes = {};
params.residual_forms = {};
params.run_name_filter = {};
params.run_name_contains = {};
params.feature_names = {};
params.method_tags = {};
params.path_kinds = {};
params.component_counts = [];
params.max_runs = [];
params.current_best_p7_only = base.current_best_p7_only;
params.current_p7_run_names = base.current_p7_run_names;
params.current_p7_autodl_roots = base.current_p7_autodl_roots;

params.density_sources = struct([]);
params.use_default_density_triplet = true;
params.density_source_kinds = base.density_source_kinds;
params.require_all_density_sources = false;
params.blp_density_threshold_ratio = base.blp_density_threshold_ratio;
params.blp_density_condition_suffix = base.blp_density_condition_suffix;
params.blp_dimred_method_tags = base.blp_dimred_method_tags;
params.blp_dimred_method_tag = base.blp_dimred_method_tag;
params.output_mode = 'both';

params.skip_existing_xcorr = true;
params.load_existing_xcorr = false;
params.continue_on_error = true;
params.headless = true;
params.close_figures = true;
params.return_heavy_outputs = false;

params.output_folder_name = io_project.get_pipeline_stage_name(10, ...
    'bold_dimred_density_cross_correlation');
params.figure_folder_name = io_project.get_pipeline_stage_name(10, ...
    'figures_bold_dimred_top_xcorr_activation_maps');
params.manifest_dir = fullfile(params.processed_root, ...
    'postprocessing_manifests', 'pipeline10_dimred_xcorr');

params.xcorr = base.xcorr;
params.xcorr.feature_names = {'efun_real'};
params.xcorr.component_value_modes = {'real'};
params.xcorr.save_tag = 'dimred_xcorr';
params.xcorr.save_full_source_results = false;
params.xcorr.save_source_results = false;
params.xcorr.max_reusable_mat_bytes = 25 * 1024^2;

params.activation = base.activation;
params.activation.top_n = params.xcorr.top_n;

params.roi_summary = base.roi_summary;
params.roi_summary.top_n = params.xcorr.top_n;
end
