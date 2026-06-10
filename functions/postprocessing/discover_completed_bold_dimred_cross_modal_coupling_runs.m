function candidates = discover_completed_bold_dimred_cross_modal_coupling_runs(params)
%DISCOVER_COMPLETED_BOLD_DIMRED_CROSS_MODAL_COUPLING_RUNS Find P9 results for P10.

if nargin < 1 || isempty(params)
    params = build_bold_dimred_cross_modal_coupling_params();
else
    params = local_apply_defaults(params);
end

candidates = repmat(local_empty_candidate(), 0, 1);
dataset_dirs = dir(params.processed_root);
dataset_dirs = dataset_dirs([dataset_dirs.isdir]);
dataset_dirs = dataset_dirs(~ismember({dataset_dirs.name}, {'.', '..'}));

for i_ds = 1:numel(dataset_dirs)
    dataset_stem = dataset_dirs(i_ds).name;
    if any(strcmpi(dataset_stem, params.exclude_dataset_stems))
        continue;
    end
    if ~isempty(params.dataset_stems) && ...
            ~any(strcmpi(dataset_stem, params.dataset_stems))
        continue;
    end

    p9_root = io_project.get_pipeline_stage_dir(params.processed_root, ...
        dataset_stem, 9, 'bold_eigenfunction_reduction');
    if exist(p9_root, 'dir') ~= 7
        continue;
    end

    L = dir(fullfile(p9_root, '*', '*', '*', 'mat', '*.mat'));
    for i_file = 1:numel(L)
        result_file = fullfile(L(i_file).folder, L(i_file).name);
        cand = local_candidate_from_result_file(result_file, params);
        if isempty(cand)
            continue;
        end
        if ~local_matches_filters(cand, params)
            continue;
        end
        candidates(end + 1, 1) = cand; %#ok<AGROW>
    end
end

if isempty(candidates)
    return;
end

if isfield(params, 'current_best_p7_only') && logical(params.current_best_p7_only)
    current_run_names = resolve_current_bold_p7_run_names(params);
    if isempty(current_run_names)
        candidates = candidates([]);
        return;
    end
    keep = ismember(lower(string({candidates.run_name})), lower(string(current_run_names)));
    candidates = candidates(keep);
    if isempty(candidates)
        return;
    end
end

[~, order] = sort(string({candidates.dataset_stem}) + "|" + ...
    string({candidates.run_tag}) + "|" + string({candidates.feature_name}) + "|" + ...
    string({candidates.method_tag}));
candidates = candidates(order);

if ~isempty(params.max_runs)
    candidates = candidates(1:min(numel(candidates), double(params.max_runs)));
end
end


function params = local_apply_defaults(params)
defaults = build_bold_dimred_cross_modal_coupling_params();
params = local_merge_defaults(defaults, params);
list_fields = {'dataset_stems', 'exclude_dataset_stems', 'observable_modes', ...
    'residual_forms', 'run_name_filter', 'run_name_contains', ...
    'feature_names', 'method_tags', 'path_kinds', ...
    'current_p7_run_names', 'current_p7_autodl_roots'};
for i = 1:numel(list_fields)
    name = list_fields{i};
    params.(name) = cellstr(string(params.(name)(:)).');
end
end


function out = local_merge_defaults(defaults, overrides)
out = defaults;
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    value = overrides.(name);
    if isstruct(value) && isfield(defaults, name) && isstruct(defaults.(name)) && ...
            isscalar(value) && isscalar(defaults.(name))
        out.(name) = local_merge_defaults(defaults.(name), value);
    elseif ~isempty(value)
        out.(name) = value;
    else
        out.(name) = defaults.(name);
    end
end
end


function cand = local_candidate_from_result_file(result_file, params)
cand = [];
try
    S = load_mat_file_with_short_path(result_file, 'result');
catch
    return;
end
if ~isfield(S, 'result') || ~isstruct(S.result)
    return;
end
R = S.result;

run_info = local_get_field(R, 'run_info', struct());
cfg = local_get_field(R, 'cfg', struct());
source = local_get_field(R, 'source', struct());

dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', ...
    local_get_field(local_get_field(cfg, 'dataset', struct()), 'name', ''))));
if isempty(dataset_stem)
    return;
end
dataset_id = char(string(local_get_field(run_info, 'dataset_id', ...
    io_project.get_dataset_id_from_stem(dataset_stem))));
run_name = char(string(local_get_field(run_info, 'run_name', '')));
observable_mode = char(string(local_get_field(run_info, 'observable_mode', ...
    local_parse_observable_mode(run_name))));
residual_form = char(string(local_get_field(run_info, 'residual_form', ...
    local_parse_residual_form(run_name))));
feature_name = char(string(local_get_field(local_get_field(R, 'feature', struct()), ...
    'name', local_get_field(local_get_field(R, 'meta', struct()), 'feature_name', ''))));
method_tag = local_method_tag_from_result(R, result_file);
path_kind = char(string(local_get_field(local_get_field(R, 'meta', struct()), ...
    'path_kind', '')));
n_components = local_get_field(local_get_field(R, 'core', struct()), 'n_components', NaN);
bold_post_file = char(string(local_get_field(source, 'bold_post_file', '')));
if isempty(bold_post_file) || exist(bold_post_file, 'file') ~= 2
    return;
end

run_tag = local_get_field(local_get_field(cfg, 'output', struct()), 'run_tag', '');
if isempty(run_tag)
    run_tag = make_bold_cross_modal_run_tag(struct( ...
        'dataset_stem', dataset_stem, 'run_name', run_name, ...
        'observable_mode', observable_mode, 'residual_form', residual_form));
end
run_tag = char(string(run_tag));
feature_tag = regexprep(lower(feature_name), '[^\w\-]+', '_');
p10_tag = sprintf('%s__%s__%s', run_tag, feature_tag, method_tag);

cand = local_empty_candidate();
cand.dataset_stem = dataset_stem;
cand.dataset_id = dataset_id;
cand.run_name = run_name;
cand.run_tag = run_tag;
cand.p10_tag = p10_tag;
cand.observable_mode = observable_mode;
cand.residual_form = residual_form;
cand.feature_name = feature_name;
cand.feature_tag = feature_tag;
cand.method_tag = method_tag;
cand.path_kind = path_kind;
cand.n_components = double(n_components);
cand.bold_post_file = bold_post_file;
cand.dimred_result_file = result_file;
cand.xcorr_dir = fullfile(io_project.get_pipeline_stage_dir(params.processed_root, ...
    dataset_stem, 10, 'bold_dimred_density_cross_correlation'), p10_tag);
cand.xcorr_file = fullfile(cand.xcorr_dir, [params.xcorr.save_tag, '.mat']);
cand.activation_root = fullfile(io_project.get_pipeline_stage_dir(params.processed_root, ...
    dataset_stem, 10, 'figures_bold_dimred_top_xcorr_activation_maps'), p10_tag);
end


function tf = local_matches_filters(cand, params)
tf = true;
if ~isempty(params.observable_modes) && ~any(strcmpi(cand.observable_mode, params.observable_modes))
    tf = false; return;
end
if ~isempty(params.residual_forms) && ~any(strcmpi(cand.residual_form, params.residual_forms))
    tf = false; return;
end
if ~isempty(params.run_name_filter) && ~any(strcmpi(cand.run_name, params.run_name_filter))
    tf = false; return;
end
for i = 1:numel(params.run_name_contains)
    if ~contains(cand.run_name, params.run_name_contains{i}, 'IgnoreCase', true)
        tf = false; return;
    end
end
if ~isempty(params.feature_names) && ~any(strcmpi(cand.feature_name, params.feature_names))
    tf = false; return;
end
if ~isempty(params.method_tags) && ~any(strcmpi(cand.method_tag, params.method_tags))
    tf = false; return;
end
if ~isempty(params.path_kinds) && ~any(strcmpi(cand.path_kind, params.path_kinds))
    tf = false; return;
end
if ~isempty(params.component_counts) && ...
        ~any(double(cand.n_components) == double(params.component_counts(:)))
    tf = false; return;
end
end


function method_tag = local_method_tag_from_result(R, result_file)
cfg = local_get_field(R, 'cfg', struct());
output = local_get_field(cfg, 'output', struct());
method_tag = char(string(local_get_field(output, 'method_tag', '')));
if ~isempty(method_tag)
    return;
end
artifacts = local_get_field(R, 'artifacts', struct());
result_path = char(string(local_get_field(artifacts, 'result_mat_file', result_file)));
[~, base] = fileparts(result_path);
tokens = regexp(base, '_(svd|logsvd|nmf|mds|umap)_k\d+$', 'match');
if ~isempty(tokens)
    method_tag = regexprep(tokens{1}, '^_', '');
else
    method_tag = char(string(local_get_field(local_get_field(R, 'meta', struct()), 'method', 'method')));
    method_tag = regexprep(lower(method_tag), '[^\w\-]+', '_');
end
end


function cand = local_empty_candidate()
cand = struct('dataset_stem', '', 'dataset_id', '', 'run_name', '', ...
    'run_tag', '', 'p10_tag', '', 'observable_mode', '', 'residual_form', '', ...
    'feature_name', '', 'feature_tag', '', 'method_tag', '', 'path_kind', '', ...
    'n_components', NaN, 'bold_post_file', '', 'dimred_result_file', '', ...
    'xcorr_dir', '', 'xcorr_file', '', 'activation_root', '');
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
    value = S.(name);
else
    value = default_value;
end
end
