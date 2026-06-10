function out = compute_bold_dimred_efun_density_cross_correlation(dimred_result_input, density_sources, params)
%COMPUTE_BOLD_DIMRED_EFUN_DENSITY_CROSS_CORRELATION P10 component-density xcorr.

if nargin < 1 || isempty(dimred_result_input)
    error('dimred_result_input is required.');
end
if nargin < 2
    density_sources = struct([]);
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = apply_bold_efun_density_cross_correlation_defaults(params);
params = local_apply_component_defaults(params);

result = local_load_dimred_result(dimred_result_input);
pseudo_B = local_build_component_bold_post(result, params);

ctx = prepare_bold_efun_density_cross_correlation_context( ...
    pseudo_B, density_sources, params);
stats = compute_bold_efun_density_cross_correlation_source_results(ctx);
out = build_bold_efun_density_cross_correlation_output(ctx, stats);
out = local_relabel_component_output(out, result, params);
out = publish_bold_efun_density_cross_correlation_output(out, params);
end


function params = local_apply_component_defaults(params)
if ~isfield(params, 'component_value_modes') || isempty(params.component_value_modes)
    params.component_value_modes = {'real'};
end
modes = cellstr(string(params.component_value_modes(:)).');
feature_names = cell(size(modes));
for i = 1:numel(modes)
    switch lower(modes{i})
        case 'real'
            feature_names{i} = 'efun_real';
        case 'abs'
            feature_names{i} = 'efun_abs';
        otherwise
            error('component_value_modes entries must be ''real'' or ''abs''.');
    end
end
params.feature_names = feature_names;
end


function result = local_load_dimred_result(input)
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


function B = local_build_component_bold_post(result, params)
components = local_component_matrix(result, params);
T = size(components, 1);
n_comp = size(components, 2);
bold_post_file = local_get_field(local_get_field(result, 'source', struct()), ...
    'bold_post_file', '');
bold_post = local_load_bold_post_if_available(bold_post_file);

B = struct();
dimred_result_file = local_get_field(local_get_field(result, 'artifacts', struct()), ...
    'result_mat_file', '');
if ~isempty(bold_post_file)
    B.source_file = bold_post_file;
else
    B.source_file = dimred_result_file;
end
B.run_info = local_get_field(bold_post, 'run_info', local_get_field(result, 'run_info', struct()));
B.dt = local_get_field(bold_post, 'dt', ...
    local_get_field(local_get_field(result, 'input', struct()), 'dt', []));
if isempty(B.dt)
    B.dt = 1;
end
B.time_vec = local_get_field(bold_post, 'time_vec', ...
    local_get_field(local_get_field(result, 'input', struct()), ...
    'time_axis', (0:T-1)' * B.dt));
B.session = local_get_field(bold_post, 'session', local_session_from_result(result, T));

E = struct();
E.efuns = components;
E.evalues = ones(n_comp, 1);
E.kpm_modes = eye(n_comp);
E.dt = B.dt;
E.dx = B.dt;
E.idx_final_in_original = (1:n_comp).';
B.EDMD_outputs = E;
end


function B = local_load_bold_post_if_available(file)
B = struct();
if isempty(file) || exist(file, 'file') ~= 2
    return;
end
try
    S = load_mat_file_with_short_path(file, 'BOLD_POST');
    if isfield(S, 'BOLD_POST') && isstruct(S.BOLD_POST)
        B = S.BOLD_POST;
        B.source_file = file;
    end
catch
    B = struct();
end
end


function C = local_component_matrix(result, params)
source_name = 'raw';
if isfield(params, 'component_source') && ~isempty(params.component_source)
    source_name = char(string(params.component_source));
end
switch lower(source_name)
    case {'smooth', 'smooth_if_available'}
        if isfield(result, 'summary') && ...
                isfield(result.summary, 'temporal_components_smooth_time_by_comp') && ...
                ~isempty(result.summary.temporal_components_smooth_time_by_comp)
            C = result.summary.temporal_components_smooth_time_by_comp;
        else
            C = result.core.temporal_components_time_by_comp;
        end
    case {'raw', 'component', 'components'}
        C = result.core.temporal_components_time_by_comp;
    otherwise
        error('Unsupported component_source: %s.', source_name);
end
C = double(C);
end


function session = local_session_from_result(result, T)
session = struct();
if isfield(result, 'run_info')
    session = local_get_field(result.run_info, 'session', struct());
end
if isempty(fieldnames(session)) && isfield(result, 'input') && isfield(result.input, 'session')
    session = result.input.session;
end
if ~isfield(session, 'session_start_idx') || isempty(session.session_start_idx)
    session.session_start_idx = 1;
    session.session_end_idx = T;
    session.session_lengths = T;
    session.session_ids = 1;
    session.border_idx = [];
end
end


function out = local_relabel_component_output(out, result, params)
out.meta = local_get_field(out, 'meta', struct());
out.meta.pipeline = 10;
out.meta.source_pipeline = 9;
out.bold_dimred_result_file = local_get_field(local_get_field(result, 'artifacts', struct()), ...
    'result_mat_file', '');
out.bold_dimred_feature = local_get_field(local_get_field(result, 'feature', struct()), ...
    'name', '');
out.bold_dimred_method = local_get_field(local_get_field(result, 'meta', struct()), ...
    'method', '');
out.bold_dimred_path_kind = local_get_field(local_get_field(result, 'meta', struct()), ...
    'path_kind', '');
out.bold_dimred_n_components = local_get_field(local_get_field(result, 'core', struct()), ...
    'n_components', NaN);

out.peak_table = local_relabel_table(out.peak_table, result);
out.top_table = local_relabel_table(out.top_table, result);
out.peak_table_by_density = local_relabel_table_struct(out.peak_table_by_density, result);
out.top_table_by_density = local_relabel_table_struct(out.top_table_by_density, result);
out.peak_table_by_density_feature = local_relabel_table_struct(out.peak_table_by_density_feature, result);
out.top_table_by_density_feature = local_relabel_table_struct(out.top_table_by_density_feature, result);

if isfield(out, 'density_feature_feature_names')
    out.density_feature_feature_names = local_relabel_feature_names( ...
        out.density_feature_feature_names, result);
end
if isfield(out, 'density_feature_group_names')
    out.density_feature_group_names = local_relabel_group_names( ...
        out.density_feature_group_names, result);
end
out.source_results = local_relabel_source_results(out.source_results, result);
out.params = params;
end


function T = local_relabel_table(T, result)
if isempty(T) || ~istable(T)
    return;
end
if ismember('bold_feature', T.Properties.VariableNames)
    labels = local_relabel_feature_names(cellstr(string(T.bold_feature)), result);
    T.bold_feature = labels(:);
end
if ~ismember('bold_component_index', T.Properties.VariableNames) && ...
        ismember('bold_mode_index', T.Properties.VariableNames)
    T.bold_component_index = T.bold_mode_index;
end
T.bold_dimred_feature = repmat({local_get_field(local_get_field(result, 'feature', struct()), ...
    'name', '')}, height(T), 1);
T.bold_dimred_method = repmat({local_get_field(local_get_field(result, 'meta', struct()), ...
    'method', '')}, height(T), 1);
end


function S = local_relabel_table_struct(S, result)
if ~isstruct(S)
    return;
end
names = fieldnames(S);
for i = 1:numel(names)
    S.(names{i}) = local_relabel_table(S.(names{i}), result);
end
end


function labels = local_relabel_feature_names(labels, result)
feature_name = local_get_field(local_get_field(result, 'feature', struct()), 'name', 'bold_dimred');
labels = cellstr(string(labels(:)));
for i = 1:numel(labels)
    switch lower(labels{i})
        case 'efun_real'
            labels{i} = sprintf('%s_component_real', feature_name);
        case 'efun_abs'
            labels{i} = sprintf('%s_component_abs', feature_name);
    end
end
end


function names = local_relabel_group_names(names, result)
names = cellstr(string(names(:)).');
for i = 1:numel(names)
    names{i} = strrep(names{i}, 'efun_real', ...
        sprintf('%s_component_real', local_get_field(local_get_field(result, 'feature', struct()), 'name', 'bold_dimred')));
    names{i} = strrep(names{i}, 'efun_abs', ...
        sprintf('%s_component_abs', local_get_field(local_get_field(result, 'feature', struct()), 'name', 'bold_dimred')));
end
end


function source_results = local_relabel_source_results(source_results, result)
if isempty(source_results)
    return;
end
if iscell(source_results)
    for i = 1:numel(source_results)
        source_results{i} = local_relabel_one_source_result(source_results{i}, result);
    end
else
    for i = 1:numel(source_results)
        source_results(i) = local_relabel_one_source_result(source_results(i), result);
    end
end
end


function item = local_relabel_one_source_result(item, result)
if ~isstruct(item) || ~isfield(item, 'best_feature') || ...
        ~isstruct(item.best_feature) || ~isfield(item.best_feature, 'name')
    return;
end
labels = local_relabel_feature_names({item.best_feature.name}, result);
item.best_feature.name = labels{1};
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
