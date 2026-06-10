function [component_plot_ctx, component_info, plot_ctx, result] = ...
        build_bold_dimred_component_plot_context(dimred_result_input, bold_post_input, params)
%BUILD_BOLD_DIMRED_COMPONENT_PLOT_CONTEXT Convert P9 component loadings to spatial modes.

if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

result = local_load_result(dimred_result_input);
if nargin < 2 || isempty(bold_post_input)
    bold_post_input = local_get_field(local_get_field(result, 'source', struct()), ...
        'bold_post_file', '');
end
if isempty(bold_post_input)
    error('bold_post_input could not be resolved.');
end

[plot_ctx, ~] = build_bold_activation_plot_context(bold_post_input, struct( ...
    'datapons_root', params.datapons_root, ...
    'roi_ts_file', params.roi_ts_file, ...
    'background_limits', params.background_limits, ...
    'feature_reduce', params.feature_reduce));

component_kpm = local_component_kpm_modes(result);
component_plot_ctx = plot_ctx;
component_plot_ctx.kpm_modes = component_kpm;
component_plot_ctx.evalues = ones(size(component_kpm, 1), 1);
component_plot_ctx.run_info = plot_ctx.run_info;
component_plot_ctx.run_info.dimred_feature_name = local_get_field( ...
    local_get_field(result, 'feature', struct()), 'name', '');
component_plot_ctx.run_info.dimred_method = local_get_field( ...
    local_get_field(result, 'meta', struct()), 'method', '');

component_info = struct();
component_info.component_kpm_modes = component_kpm;
component_info.n_components = size(component_kpm, 1);
component_info.dimred_result_file = local_get_field(local_get_field(result, 'artifacts', struct()), ...
    'result_mat_file', '');
component_info.feature_name = component_plot_ctx.run_info.dimred_feature_name;
component_info.method = component_plot_ctx.run_info.dimred_method;
component_info.path_kind = local_get_field(local_get_field(result, 'meta', struct()), ...
    'path_kind', '');
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'roi_ts_file', '');
params = local_set_default(params, 'background_limits', []);
params = local_set_default(params, 'feature_reduce', 'mean');
end


function result = local_load_result(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load_mat_file_with_short_path(file, 'result');
    if ~isfield(S, 'result')
        error('result variable missing in %s.', file);
    end
    result = S.result;
    if ~isfield(result, 'artifacts') || ~isstruct(result.artifacts)
        result.artifacts = struct();
    end
    result.artifacts.result_mat_file = file;
elseif isstruct(input)
    result = input;
else
    error('dimred_result_input must be a result path or struct.');
end
end


function component_kpm = local_component_kpm_modes(result)
W = local_get_field(local_get_field(result, 'core', struct()), ...
    'mode_weights_mode_by_comp', []);
K = local_get_field(local_get_field(result, 'data', struct()), ...
    'kpm_modes_mode_by_dict', []);
if isempty(W) || isempty(K)
    error('P9 result must contain core.mode_weights_mode_by_comp and data.kpm_modes_mode_by_dict.');
end
W = double(W);
K = double(K);
n_mode = min(size(W, 1), size(K, 1));
if n_mode < 1
    error('No overlapping modes were available to build component spatial loadings.');
end
component_kpm = W(1:n_mode, :).'*K(1:n_mode, :);
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
