function out = build_bold_efun_density_cross_correlation_output(ctx, stats)
%BUILD_BOLD_EFUN_DENSITY_CROSS_CORRELATION_OUTPUT Pack pipeline 8 xcorr output struct.

out = struct();
out.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
out.bold_post_file = local_get_field(ctx.B, 'source_file', '');
out.params = ctx.params;
out.dt = ctx.dt;
out.T = ctx.T;
out.session = ctx.session;
out.max_lag_bins = ctx.max_lag_bins;
out.border_pad_bins = ctx.border_pad_bins;
out.border_mask = ctx.base_mask;
out.lag_bins = ctx.lag_bins;
out.lag_sec = ctx.lag_sec;
out.positive_lag_definition = ctx.positive_lag_definition;
out.feature_specs = ctx.feature_specs;
out.density_sources = ctx.density_sources;
out.source_results = stats.source_results;
out.peak_table = stats.peak_table;
out.top_table = stats.top_table;
out.peak_table_by_density = stats.peak_table_by_density;
out.top_table_by_density = stats.top_table_by_density;
out.density_group_names = stats.density_group_names;
out.density_group_fields = stats.density_group_fields;
out.peak_table_by_density_feature = stats.peak_table_by_density_feature;
out.top_table_by_density_feature = stats.top_table_by_density_feature;
out.density_feature_group_names = stats.density_feature_group_names;
out.density_feature_group_fields = stats.density_feature_group_fields;
out.density_feature_density_names = stats.density_feature_density_names;
out.density_feature_feature_names = stats.density_feature_feature_names;
out.output_mode = local_output_mode_label(ctx.params);
end


function label = local_output_mode_label(params)
if params.export_combined && params.export_by_density
    label = 'both';
elseif params.export_combined
    label = 'combined';
elseif params.export_by_density
    label = 'separate';
else
    label = 'none';
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
