function params = apply_bold_efun_density_cross_correlation_defaults(params)
%APPLY_BOLD_EFUN_DENSITY_CROSS_CORRELATION_DEFAULTS Fill canonical pipeline 8 defaults.

if nargin < 1 || isempty(params)
    params = struct();
end

params = local_set_default(params, 'max_lag_sec', 10);
params = local_set_default(params, 'border_mask_sec', params.max_lag_sec);
params = local_set_default(params, 'min_valid_samples', 20);
params = local_set_default(params, 'top_n', 5);
params = local_set_default(params, 'feature_names', ...
    {'efun_abs', 'efun_real', 'deconv_abs', 'deconv_real'});
params = local_set_default(params, 'density_field_preference', ...
    {'smoothed_density_mean', 'density_mean', ...
    'density_time_by_mode', 'density_time_by_component', 'density'});
params = local_set_default(params, 'density_transform', 'zscore');
params = local_set_default(params, 'bold_transform', 'zscore');
params = local_set_default(params, 'save_results', true);
params = local_set_default(params, 'save_mat', true);
params = local_set_default(params, 'save_source_results', true);
params = local_set_default(params, 'save_dir', pwd);
params = local_set_default(params, 'save_tag', 'xcorr');
params = local_set_default(params, 'make_figures', true);
params = local_set_default(params, 'export_combined', true);
params = local_set_default(params, 'export_by_density', true);
params = local_set_default(params, 'export_by_density_feature', true);

if ~isfield(params, 'plot') || isempty(params.plot)
    params.plot = struct();
end
params.plot = local_set_default(params.plot, 'save_dir', params.save_dir);
params.plot = local_set_default(params.plot, 'save_tag', params.save_tag);
params.plot = local_set_default(params.plot, 'top_n', params.top_n);
params.plot = local_set_default(params.plot, 'resolution', 220);
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end
