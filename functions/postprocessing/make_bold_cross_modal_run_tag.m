function tag = make_bold_cross_modal_run_tag(run_info)
%MAKE_BOLD_CROSS_MODAL_RUN_TAG Short, stable folder tag for pipeline 8 BOLD runs.

if nargin < 1 || isempty(run_info)
    run_info = struct();
end

if ischar(run_info) || isstring(run_info)
    run_name = char(string(run_info));
    observable_mode = local_parse_observable_mode(run_name);
    residual_form = local_parse_residual_form(run_name);
elseif isstruct(run_info)
    run_name = local_get_field(run_info, 'run_name', '');
    observable_mode = local_get_field(run_info, 'observable_mode', '');
    residual_form = local_get_field(run_info, 'residual_form', '');
    if isempty(observable_mode) && ~isempty(run_name)
        observable_mode = local_parse_observable_mode(run_name);
    end
    if isempty(residual_form) && ~isempty(run_name)
        residual_form = local_parse_residual_form(run_name);
    end
else
    error('run_info must be a run-name string or struct.');
end

obs_tag = local_observable_tag(observable_mode);
res_tag = local_residual_tag(residual_form);
if isempty(obs_tag)
    obs_tag = local_filename_safe(run_name);
end
if isempty(res_tag)
    tag = obs_tag;
else
    tag = sprintf('%s_%s', res_tag, obs_tag);
end
end


function tag = local_observable_tag(mode)
switch lower(char(string(mode)))
    case 'hp_svd100'
        tag = 'hp100';
    case 'global_svd100'
        tag = 'gsvd100';
    case 'global_slow_band_power_svd100'
        tag = 'gsbp100';
    case 'roi_mean_slow_band_power'
        tag = 'roi_sbp';
    case 'slow_band_power_svd'
        tag = 'sbp_svd';
    case 'roi_mean'
        tag = 'roi';
    case 'elehp'
        tag = 'elehp';
    case 'svd'
        tag = 'svd';
    case 'hp'
        tag = 'hp';
    otherwise
        tag = local_filename_safe(mode);
end
end


function tag = local_residual_tag(form)
switch lower(char(string(form)))
    case 'projected_vlambda'
        tag = 'pv';
    case 'projected_kv'
        tag = 'pkv';
    otherwise
        tag = local_filename_safe(form);
end
end


function mode = local_parse_observable_mode(run_name)
known = {'global_slow_band_power_svd100', 'roi_mean_slow_band_power', ...
    'slow_band_power_svd', 'gsvd100_ds', 'global_svd100', 'HP_svd100', ...
    'slow_band_power', 'roi_mean', 'eleHP', 'HP', 'svd'};
mode = '';
for i = 1:numel(known)
    if contains(run_name, ['_' known{i}], 'IgnoreCase', true) || ...
            endsWith(run_name, known{i}, 'IgnoreCase', true)
        mode = known{i};
        return;
    end
end
end


function form = local_parse_residual_form(run_name)
forms = {'projected_vlambda', 'projected_kv'};
form = '';
for i = 1:numel(forms)
    if contains(run_name, forms{i}, 'IgnoreCase', true)
        form = forms{i};
        return;
    end
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = char(string(S.(name)));
else
    value = default_value;
end
end


function tag = local_filename_safe(value)
tag = lower(char(string(value)));
tag = regexprep(tag, '[^a-z0-9]+', '_');
tag = regexprep(tag, '^_+|_+$', '');
if isempty(tag)
    tag = 'run';
end
end
