function params = build_bold_reskoopnet_postprocessing_params()
%BUILD_BOLD_RESKOOPNET_POSTPROCESSING_PARAMS Default parameters for pipeline 7.

params = struct();
roi_partition = build_bold_cortical_subcortical_partition_defaults();

% Discovery / run selection
params.autodl_roots = {'E:\autodl_results_local\bold_wsl'};
params.processed_root = io_project.get_project_processed_root();
params.dataset_stems = {};
params.exclude_dataset_stems = {};
params.observable_modes = {};
params.residual_forms = {};
params.run_name_filter = {};
params.run_name_contains = {};
params.filename_pattern = '*_outputs_*.mat';
params.require_summary_file = true;
params.allow_summary_only_outputs = false;
params.include_smoke = false;
params.dedupe_by_condition = true;
params.run_selection_mode = 'best_val_metric';
params.training_state_priority = {'final', 'best'};

% Pipeline 7 outputs
params.output_folder_name = io_project.get_pipeline_stage_name(7, 'bold_postprocessing');
params.manifest_dir = fullfile(params.processed_root, ...
    'postprocessing_manifests', params.output_folder_name);

% EDMD chunk concatenation
params.concat = struct();
params.concat.filename_pattern = '*_outputs_*.mat';
params.concat.variable_name = 'EDMD_outputs';
params.concat.concat_fields = {'efuns'};
params.concat.scan_mode = 'uniform_except_last';
params.concat.verbose = true;

% Main efun postprocessing
params.post_opts = struct();
params.post_opts.abs_thresh = 0.01;
params.post_opts.sort_by = 'modulus';
params.post_opts.sort_dir = 'descend';
params.post_opts.max_basis = 80;
params.post_opts.max_plot_samples = 2000;

% Deconvolution
params.deconv = struct();
params.deconv.method = 'koopman_residual';
params.deconv.lambda_source = 'edmd';
params.deconv.lambdaType = 'discrete';
params.deconv.max_modes_all = 80;
params.deconv.max_modes_sel = 40;
params.deconv.max_plot_samples = 2000;
params.deconv.plot_normalize_scope = 'global';

% Timescale plots
params.timescale_max_modes_all = Inf;
params.timescale_max_modes_sel = params.post_opts.max_basis;
params.timescale_max_lag = [];
params.timescale_xlim_time = 240;
params.timescale_match_empirical_scale = true;
params.timescale = struct();
params.timescale.max_modes_all = params.timescale_max_modes_all;
params.timescale.max_modes_sel = params.timescale_max_modes_sel;
params.timescale.maxLag = params.timescale_max_lag;
params.timescale.xlim_time = params.timescale_xlim_time;
params.timescale.match_empirical_scale_to_theoretical = ...
    params.timescale_match_empirical_scale;

% Intrinsic activation maps
params.intrinsic_activation = struct();
params.intrinsic_activation.enabled = false;
params.intrinsic_activation.selection_mode = 'sorted';
params.intrinsic_activation.basis_indices = 1:5;
params.intrinsic_activation.value_mode = 'abs';
params.intrinsic_activation.slice_list = 1:20;
params.intrinsic_activation.tiles_per_row = 10;
params.intrinsic_activation.overlay_alpha = 0.86;
params.intrinsic_activation.feature_reduce = 'mean';
params.intrinsic_activation.skip_existing = true;
params.intrinsic_activation.save_png = true;
params.intrinsic_activation.save_fig = false;
params.intrinsic_activation.resolution = 220;
params.intrinsic_activation.update_bold_post = true;

% Intrinsic ROI bar summaries
params.intrinsic_roi_summary = struct();
params.intrinsic_roi_summary.enabled = false;
params.intrinsic_roi_summary.selection_mode = 'sorted';
params.intrinsic_roi_summary.basis_indices = 1:5;
params.intrinsic_roi_summary.feature_reduce = 'mean';
params.intrinsic_roi_summary.roi_reduce = 'mean';
params.intrinsic_roi_summary.roi_value_mode = 'mean_abs';
params.intrinsic_roi_summary.mode_normalization = 'range';
params.intrinsic_roi_summary.flip_roi_order = true;
params.intrinsic_roi_summary.layout_mode = 'row';
params.intrinsic_roi_summary.show_roi_labels = 'first';
params.intrinsic_roi_summary.highlight_regions = true;
params.intrinsic_roi_summary.highlight_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
params.intrinsic_roi_summary.highlight_contains_tokens = {'V1'};
params.intrinsic_roi_summary.highlight_row_color = [0.10 0.10 0.10];
params.intrinsic_roi_summary.highlight_row_alpha = 0.08;
params.intrinsic_roi_summary.highlight_label_prefix = '> ';
params.intrinsic_roi_summary.show_cortical_subcortical_separator = true;
params.intrinsic_roi_summary.separator_after_roi_name = roi_partition.separator_after_roi_name;
params.intrinsic_roi_summary.cortical_exact_names = roi_partition.cortical_exact_names;
params.intrinsic_roi_summary.cortical_contains_tokens = roi_partition.cortical_contains_tokens;
params.intrinsic_roi_summary.subcortical_exact_names = roi_partition.subcortical_exact_names;
params.intrinsic_roi_summary.subcortical_contains_tokens = roi_partition.subcortical_contains_tokens;
params.intrinsic_roi_summary.skip_existing = true;
params.intrinsic_roi_summary.save_png = true;
params.intrinsic_roi_summary.save_fig = false;
params.intrinsic_roi_summary.resolution = 220;
params.intrinsic_roi_summary.update_bold_post = true;

% Execution / plotting behavior
params.datapons_root = 'E:\DataPons';
params.default_dt = 2;
params.make_main_plot = true;
params.compute_deconv = true;
params.make_deconv_plot = true;
params.make_timescale_plot = true;
params.draw_session_borders = true;
params.save_png = true;
params.save_fig = false;
params.resolution = 180;
params.figure_export_method = 'print';
params.figure_renderer = 'opengl';
params.skip_existing = true;
params.continue_on_error = true;
params.headless = true;
params.close_figures = true;
params.verbose = true;
end
