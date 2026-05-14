function result = build_blp_eigenfunction_reduction_result( ...
        prep, reduction, method_name, cfg, source_info, concat_info)
%BUILD_BLP_EIGENFUNCTION_REDUCTION_RESULT Pack core reduction outputs into result schema.

result = struct();

result.meta = struct();
result.meta.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
result.meta.path_kind = cfg.path.kind;
result.meta.feature_family = 'eigenfunction';
result.meta.feature_variant = cfg.feature.variant;
result.meta.method = method_name;

result.cfg = strip_blp_preloaded_eigenfunction_source_payload(cfg);
result.source = source_info;
result.concat = concat_info;

result.input = struct();
result.input.dt = prep.dt;
result.input.dt_source = prep.dt_source;
result.input.time_axis = prep.time_axis;
result.input.mode_index = (1:numel(prep.evalues_discrete)).';
result.input.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;
result.input.selected_mode_mask_in_original = prep.selected_mode_mask_in_original;

result.data = struct();
result.data.evalues_discrete = prep.evalues_discrete;
result.data.evalues_bilinear = prep.evalues_bilinear;
result.data.efun_raw_time_by_mode = prep.efun_raw_time_by_mode;
result.data.efun_feature_time_by_mode = prep.efun_feature_time_by_mode;
result.data.kpm_modes_mode_by_dict = prep.kpm_modes_mode_by_dict;

result.feature = struct();
result.feature.family = 'eigenfunction';
result.feature.variant = cfg.feature.variant;
result.feature.normalization = cfg.feature.normalization;
result.feature.axis_order = 'time_by_mode';

result.core = reduction.core;
result.quality = reduction.quality;
result.aux = reduction.aux;
result.summary = build_blp_eigenfunction_reduction_summary(result.core, cfg.summary);

result.artifacts = struct();
result.artifacts.result_mat_file = '';
result.artifacts.save_payload = cfg.save.payload;
result.artifacts.thresholded_density_mat_file = '';
result.artifacts.thresholded_density_figure_file = '';
result.artifacts.thresholded_events_mat_file = '';
result.artifacts.thresholded_events_figure_file = '';
result.artifacts.dimred_thresholded_density_mat_file = '';
result.artifacts.dimred_thresholded_density_figure_file = '';
result.artifacts.dimred_thresholded_events_mat_file = '';
result.artifacts.dimred_thresholded_events_figure_file = '';

result.thresholded_density = struct();
result.thresholded_events = struct();
result.dimred_thresholded_density = struct();
result.dimred_thresholded_events = struct();
end
