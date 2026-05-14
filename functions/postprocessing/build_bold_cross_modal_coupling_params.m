function params = build_bold_cross_modal_coupling_params()
%BUILD_BOLD_CROSS_MODAL_COUPLING_PARAMS Default parameters for pipeline 8.

params = struct();
roi_partition = build_bold_cortical_subcortical_partition_defaults();

% Discovery / run selection
params.processed_root = io_project.get_project_processed_root();
params.datapons_root = 'E:\DataPons';
params.dataset_stems = {};
params.exclude_dataset_stems = {};
params.observable_modes = {};
params.residual_forms = {};
params.run_name_filter = {};
params.run_name_contains = {};
params.max_runs = [];

% Shared density input override
params.density_sources = struct([]);
params.use_default_density_triplet = true;
params.density_source_kinds = { ...
    'event_density', ...
    'raw_abs_density', ...
    'dimred_abs_density', ...
    'raw_complex_split_density', ...
    'dimred_complex_split_density'};
params.require_all_density_sources = false;
params.blp_density_threshold_ratio = 0.70;
params.blp_dimred_method_tag = 'umap_k08';
params.output_mode = 'both';

% Execution behavior
params.skip_existing_xcorr = true;
params.load_existing_xcorr = false;
params.continue_on_error = true;
params.headless = false;
params.close_figures = false;

% Output bookkeeping
params.output_folder_name = io_project.get_pipeline_stage_name(8, 'efun_density_cross_correlation');
params.figure_folder_name = io_project.get_pipeline_stage_name(8, 'figures_bold_top_xcorr_activation_maps');
params.manifest_dir = fullfile(params.processed_root, ...
    'postprocessing_manifests', 'pipeline8_xcorr');

% XCORR stage
params.xcorr = struct();
params.xcorr.max_lag_sec = 10;
params.xcorr.border_mask_sec = 10;
params.xcorr.min_valid_samples = 20;
params.xcorr.top_n = 5;
params.xcorr.feature_names = {'efun_abs', 'efun_real', 'deconv_abs', 'deconv_real'};
params.xcorr.save_results = true;
params.xcorr.save_mat = true;
params.xcorr.make_figures = true;
params.xcorr.save_tag = 'xcorr';
params.xcorr.export_combined = true;
params.xcorr.export_by_density = true;
params.xcorr.export_by_density_feature = true;

% Activation-map stage
params.activation = struct();
params.activation.enabled = true;
params.activation.top_n = params.xcorr.top_n;
params.activation.value_mode = 'abs';
params.activation.slice_list = 1:20;
params.activation.tiles_per_row = 10;
params.activation.overlay_alpha = 0.86;
params.activation.feature_reduce = 'mean';
params.activation.highlight_core_rois = true;
params.activation.core_roi_observable_modes = {'svd', 'global_svd100', 'gsvd100_ds', ...
    'global_slow_band_power_svd100', 'roi_mean', 'roi_mean_slow_band_power'};
params.activation.core_roi_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
params.activation.core_roi_contains_tokens = {'V1'};
params.activation.core_roi_outline_color = [1 1 1];
params.activation.core_roi_line_width = 1.6;
params.activation.annotate_regions = false;
params.activation.annotation_observable_modes = {};
params.activation.annotation_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
params.activation.annotation_contains_tokens = {'V1'};
params.activation.annotation_repeat_each_slice = false;
params.activation.annotation_min_voxels_per_slice = 5;
params.activation.annotation_text_color = [1 1 1];
params.activation.annotation_font_size = 8;
params.activation.annotation_font_weight = 'bold';
params.activation.annotation_background_color = [0 0 0];
params.activation.annotation_margin = 1;
params.activation.annotation_text_offset = [0 -2];
params.activation.annotation_min_text_spacing = 5;
params.activation.annotation_show_marker = false;
params.activation.outline_regions = false;
params.activation.outline_observable_modes = {};
params.activation.outline_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
params.activation.outline_contains_tokens = {'V1'};
params.activation.outline_use_annotation_names = true;
params.activation.outline_line_width = 1.6;
params.activation.outline_colors = [];
params.activation.outline_show_legend = true;
params.activation.roi_ts_file = '';
params.activation.background_limits = [];
params.activation.skip_existing = true;
params.activation.export_combined = true;
params.activation.export_by_density = true;
params.activation.export_by_density_feature = true;
params.activation.save_png = true;
params.activation.save_fig = false;
params.activation.resolution = 220;

% ROI bar-summary stage
params.roi_summary = struct();
params.roi_summary.enabled = true;
params.roi_summary.top_n = params.xcorr.top_n;
params.roi_summary.feature_reduce = 'mean';
params.roi_summary.roi_reduce = 'mean';
params.roi_summary.roi_value_mode = 'mean_abs';
params.roi_summary.mode_normalization = 'range';
params.roi_summary.flip_roi_order = true;
params.roi_summary.layout_mode = 'row';
params.roi_summary.show_roi_labels = 'first';
params.roi_summary.highlight_regions = true;
params.roi_summary.highlight_exact_names = {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'};
params.roi_summary.highlight_contains_tokens = {'V1'};
params.roi_summary.highlight_row_color = [0.10 0.10 0.10];
params.roi_summary.highlight_row_alpha = 0.08;
params.roi_summary.highlight_label_prefix = '> ';
params.roi_summary.show_cortical_subcortical_separator = true;
params.roi_summary.separator_after_roi_name = roi_partition.separator_after_roi_name;
params.roi_summary.cortical_exact_names = roi_partition.cortical_exact_names;
params.roi_summary.cortical_contains_tokens = roi_partition.cortical_contains_tokens;
params.roi_summary.subcortical_exact_names = roi_partition.subcortical_exact_names;
params.roi_summary.subcortical_contains_tokens = roi_partition.subcortical_contains_tokens;
params.roi_summary.skip_existing = true;
params.roi_summary.export_combined = true;
params.roi_summary.export_by_density = true;
params.roi_summary.export_by_density_feature = true;
params.roi_summary.save_png = true;
params.roi_summary.save_fig = false;
params.roi_summary.resolution = 220;
end
