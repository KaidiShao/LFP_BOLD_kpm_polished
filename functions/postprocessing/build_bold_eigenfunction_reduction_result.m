function result = build_bold_eigenfunction_reduction_result(prep, reduction, method_name, cfg, B)
%BUILD_BOLD_EIGENFUNCTION_REDUCTION_RESULT Pack P9 reduction output.

result = struct();

result.meta = struct();
result.meta.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
result.meta.pipeline = 9;
result.meta.path_kind = cfg.path.kind;
result.meta.feature_family = 'bold_eigenfunction';
result.meta.feature_name = prep.feature.name;
result.meta.feature_source = prep.feature.source_kind;
result.meta.feature_variant = prep.feature.variant;
result.meta.method = method_name;

result.cfg = cfg;
if isfield(result.cfg, 'source') && isfield(result.cfg.source, 'preloaded_BOLD_POST')
    result.cfg.source = rmfield(result.cfg.source, 'preloaded_BOLD_POST');
end

result.run_info = local_get_field(B, 'run_info', struct());
result.source = struct();
result.source.bold_post_file = local_get_field(B, 'source_file', cfg.source.bold_post_file);
result.source.observable_file = local_get_field(B, 'observable_file', prep.observable_file);

result.input = struct();
result.input.dt = prep.dt;
result.input.dt_source = prep.dt_source;
result.input.time_axis = prep.time_axis;
result.input.mode_index = (1:numel(prep.evalues_discrete)).';
result.input.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;
result.input.selected_mode_idx_in_feature = prep.selected_mode_idx_in_feature;
result.input.selected_mode_mask_in_original = prep.selected_mode_mask_in_original;

result.data = struct();
result.data.evalues_discrete = prep.evalues_discrete;
result.data.evalues_bilinear = prep.evalues_bilinear;
result.data.bold_feature_raw_time_by_mode = prep.bold_feature_raw_time_by_mode;
result.data.bold_feature_time_by_mode = prep.bold_feature_time_by_mode;
result.data.kpm_modes_mode_by_dict = prep.kpm_modes_mode_by_dict;

result.feature = prep.feature;
result.feature.normalization = cfg.feature.normalization;

result.core = reduction.core;
result.quality = reduction.quality;
result.aux = reduction.aux;
result.summary = build_blp_eigenfunction_reduction_summary(result.core, struct( ...
    'smooth', struct('enable', true, 'method', 'movmean', 'window', 10)));

result.artifacts = struct();
result.artifacts.result_mat_file = '';
result.artifacts.summary_png_file = '';
result.artifacts.summary_fig_file = '';
result.artifacts.save_payload = cfg.save.payload;
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
