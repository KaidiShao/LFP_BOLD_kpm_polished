function [plot_ctx, B] = build_bold_activation_plot_context(bold_post_input, params)
%BUILD_BOLD_ACTIVATION_PLOT_CONTEXT Build a reusable spatial plotting context from BOLD_POST.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

B = local_load_bold_post(bold_post_input);
run_info = local_normalize_run_info(B);
Obs = local_resolve_observable_struct(B, run_info);
roi_ts_file = local_resolve_roi_ts_file(B, run_info, params);

EDMD_outputs = B.EDMD_outputs;
if ~isfield(EDMD_outputs, 'kpm_modes') || isempty(EDMD_outputs.kpm_modes)
    error('BOLD_POST.EDMD_outputs.kpm_modes is required for activation maps.');
end
if ~isfield(EDMD_outputs, 'evalues') || isempty(EDMD_outputs.evalues)
    error('BOLD_POST.EDMD_outputs.evalues is required for activation maps.');
end

plot_ctx = struct();
plot_ctx.run_info = run_info;
plot_ctx.observable_file = local_get_field(B, 'observable_file', '');
plot_ctx.roi_ts_file = roi_ts_file;
plot_ctx.evalues = EDMD_outputs.evalues(:);
plot_ctx.kpm_modes = EDMD_outputs.kpm_modes;
plot_ctx.observable = Obs;
plot_ctx.target_region_names = local_target_region_names_from_observable(Obs);
plot_ctx.feature_reduce = char(string(local_get_field(params, 'feature_reduce', 'mean')));

plot_ctx.direct_voxel_mode = true;
plot_ctx.coeff_t = [];
if isfield(Obs, 'O') && isstruct(Obs.O) && isfield(Obs.O, 'model') && ...
        isstruct(Obs.O.model) && isfield(Obs.O.model, 'coeff') && ...
        size(Obs.O.model.coeff, 2) == size(plot_ctx.kpm_modes, 2)
    plot_ctx.direct_voxel_mode = false;
    plot_ctx.coeff_t = double(Obs.O.model.coeff).';
end

R = load(roi_ts_file, 'roiTs');
if ~isfield(R, 'roiTs')
    error('roiTs variable missing in %s.', roi_ts_file);
end
[plot_ctx.coords, plot_ctx.region_labels, plot_ctx.n_var_region, plot_ctx.ana, ...
    plot_ctx.voxel_variable_labels, plot_ctx.voxel_region_idx, ...
    plot_ctx.voxel_within_region_idx] = ...
    local_collect_roi_coords_and_anatomy(R.roiTs, plot_ctx.target_region_names);
[plot_ctx.source_variable_labels, plot_ctx.source_label_origin] = ...
    local_resolve_source_variable_labels(Obs, plot_ctx.direct_voxel_mode);
if isempty(plot_ctx.source_variable_labels) && ...
        size(plot_ctx.kpm_modes, 2) == numel(plot_ctx.voxel_variable_labels)
    plot_ctx.source_variable_labels = plot_ctx.voxel_variable_labels;
    plot_ctx.source_label_origin = 'roiTs voxel-order fallback';
end
if isempty(plot_ctx.source_variable_labels)
    error(['Could not resolve activation source-variable labels for %s. ' ...
        'Observable metadata is missing the labels needed to map mode weights into voxel space.'], ...
        run_info.run_name);
end
[plot_ctx.source_to_voxel_indices, plot_ctx.source_space_kind, plot_ctx.mapping_summary] = ...
    local_build_source_to_voxel_map(plot_ctx.source_variable_labels, ...
    plot_ctx.region_labels, plot_ctx.n_var_region, plot_ctx.voxel_variable_labels);

if isfield(params, 'background_limits') && ~isempty(params.background_limits)
    plot_ctx.bg_limits = params.background_limits;
else
    plot_ctx.bg_limits = local_robust_limits(double(plot_ctx.ana(:)));
end
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'roi_ts_file', '');
params = local_set_default(params, 'background_limits', []);
params = local_set_default(params, 'feature_reduce', 'mean');
end


function B = local_load_bold_post(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load_mat_file_with_short_path(file, 'BOLD_POST');
    if ~isfield(S, 'BOLD_POST')
        error('BOLD_POST variable missing in %s.', file);
    end
    B = S.BOLD_POST;
    B.source_file = file;
elseif isstruct(input)
    B = input;
    if ~isfield(B, 'source_file')
        B.source_file = '';
    end
else
    error('bold_post_input must be a path or struct.');
end

if ~isfield(B, 'EDMD_outputs') || ~isstruct(B.EDMD_outputs)
    error('BOLD post input does not contain EDMD_outputs.');
end
end


function run_info = local_normalize_run_info(B)
if ~isfield(B, 'run_info') || ~isstruct(B.run_info)
    error('BOLD_POST.run_info is required.');
end

run_info = B.run_info;
run_info.dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', '')));
run_info.run_name = char(string(local_get_field(run_info, 'run_name', '')));
if isempty(run_info.dataset_stem) || isempty(run_info.run_name)
    error('BOLD_POST.run_info.dataset_stem and run_name are required.');
end
if ~isfield(run_info, 'dataset_id') || isempty(run_info.dataset_id)
    run_info.dataset_id = io_project.get_dataset_id_from_stem(run_info.dataset_stem);
end
run_info.dataset_id = char(string(run_info.dataset_id));
if ~isfield(run_info, 'observable_mode') || isempty(run_info.observable_mode)
    run_info.observable_mode = local_parse_observable_mode(run_info.run_name);
end
if ~isfield(run_info, 'residual_form') || isempty(run_info.residual_form)
    run_info.residual_form = local_parse_residual_form(run_info.run_name);
end
end


function Obs = local_resolve_observable_struct(B, run_info)
observable_file = local_get_field(B, 'observable_file', '');
if isempty(observable_file) || exist(observable_file, 'file') ~= 2
    observable_file = local_infer_observable_file(run_info);
end
Obs = local_get_field(B, 'observable', struct());
if isstruct(Obs)
    Obs = local_fill_observable_struct_defaults(Obs);
else
    Obs = struct();
end
Obs = local_attach_observable_source_file(Obs, observable_file);

if ~isempty(observable_file) && exist(observable_file, 'file') == 2
    available = whos('-file', observable_file);
    available_names = {available.name};
    load_names = intersect({'O', 'cfg', 'cfg_saved', 'params', 'obs_info', ...
        'observable_info', 'observable_labels', 'dataset_id', 'file_stem'}, ...
        available_names, 'stable');
    file_obs = load(observable_file, load_names{:});
    file_obs = local_fill_observable_struct_defaults(file_obs);
    file_obs = local_attach_observable_source_file(file_obs, observable_file);
    Obs = local_merge_missing_struct_fields(Obs, file_obs);
end

if isfield(Obs, 'O') || isfield(Obs, 'cfg') || ...
        isfield(Obs, 'observable_info') || isfield(Obs, 'obs_info')
    return;
end

error('BOLD_POST does not contain a usable observable struct or observable_file.');
end


function Obs = local_fill_observable_struct_defaults(Obs)
if ~isfield(Obs, 'O')
    Obs.O = struct();
end
if ~isfield(Obs, 'cfg')
    if isfield(Obs, 'cfg_saved') && isstruct(Obs.cfg_saved)
        Obs.cfg = Obs.cfg_saved;
    else
        Obs.cfg = struct();
    end
end
if ~isfield(Obs, 'params')
    Obs.params = struct();
end
end


function Obs = local_attach_observable_source_file(Obs, observable_file)
if nargin < 2 || isempty(observable_file)
    return;
end
if ~isfield(Obs, 'source_file') || isempty(Obs.source_file)
    Obs.source_file = observable_file;
end
end


function observable_file = local_infer_observable_file(run_info)
observable_file = '';
if nargin < 1 || ~isstruct(run_info)
    return;
end

dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', '')));
dataset_id = char(string(local_get_field(run_info, 'dataset_id', '')));
observable_mode = char(string(local_get_field(run_info, 'observable_mode', '')));
if isempty(dataset_stem) || isempty(observable_mode)
    return;
end
if isempty(dataset_id)
    dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
end

obs_dir = io_project.get_pipeline_stage_dir( ...
    io_project.get_project_processed_root(), dataset_stem, 3, 'bold_observables');
candidate = fullfile(obs_dir, ...
    sprintf('%s_bold_observables_%s.mat', dataset_id, observable_mode));
if exist(candidate, 'file') == 2
    observable_file = candidate;
    return;
end

L = dir(fullfile(obs_dir, sprintf('*_bold_observables_%s.mat', observable_mode)));
if isscalar(L)
    observable_file = fullfile(L(1).folder, L(1).name);
end
end


function out = local_merge_missing_struct_fields(out, fallback)
if nargin < 1 || ~isstruct(out)
    out = struct();
end
if nargin < 2 || ~isstruct(fallback)
    return;
end

if ~isscalar(out) || ~isscalar(fallback)
    if local_is_empty_value(out) && ~local_is_empty_value(fallback)
        out = fallback;
    end
    return;
end

names = fieldnames(fallback);
for i_name = 1:numel(names)
    name = names{i_name};
    if ~isfield(out, name) || local_is_empty_value(out.(name))
        out.(name) = fallback.(name);
    elseif isstruct(out.(name)) && isstruct(fallback.(name)) && ...
            isscalar(out.(name)) && isscalar(fallback.(name))
        out.(name) = local_merge_missing_struct_fields(out.(name), fallback.(name));
    end
end
end


function tf = local_is_empty_value(value)
if isstruct(value)
    tf = isempty(value);
else
    tf = isempty(value);
end
end


function roi_ts_file = local_resolve_roi_ts_file(B, run_info, params)
roi_ts_file = char(string(local_get_field(params, 'roi_ts_file', '')));
if ~isempty(roi_ts_file)
    if exist(roi_ts_file, 'file') ~= 2
        error('roi_ts_file does not exist: %s', roi_ts_file);
    end
    return;
end

if isfield(B, 'activation_context') && isstruct(B.activation_context) && ...
        isfield(B.activation_context, 'roi_ts_file') && ...
        ~isempty(B.activation_context.roi_ts_file) && ...
        exist(B.activation_context.roi_ts_file, 'file') == 2
    roi_ts_file = B.activation_context.roi_ts_file;
    return;
end

roi_dir = fullfile(params.datapons_root, run_info.dataset_id, 'roits');
if exist(roi_dir, 'dir') ~= 7
    error('ROI directory not found: %s', roi_dir);
end

L = dir(fullfile(roi_dir, sprintf('%s_*_roits.mat', lower(run_info.dataset_stem))));
if isempty(L)
    L = dir(fullfile(roi_dir, '*_roits.mat'));
end
if isempty(L)
    error('No roiTs MAT files found in %s.', roi_dir);
end
[~, order] = sort({L.name});
L = L(order);
roi_ts_file = fullfile(L(1).folder, L(1).name);
end


function target_region_names = local_target_region_names_from_observable(Obs)
target_region_names = {};
if isfield(Obs, 'cfg') && isstruct(Obs.cfg) && isfield(Obs.cfg, 'bold') && ...
        isstruct(Obs.cfg.bold) && isfield(Obs.cfg.bold, 'selected_region_names')
    target_region_names = cellstr(string(Obs.cfg.bold.selected_region_names(:)).');
end
if isempty(target_region_names)
    target_region_names = local_region_names_from_observable_metadata(Obs);
end
end


function target_region_names = local_region_names_from_observable_metadata(Obs)
target_region_names = {};

obs_info = [];
if isfield(Obs, 'observable_info') && istable(Obs.observable_info)
    obs_info = Obs.observable_info;
elseif isfield(Obs, 'obs_info') && istable(Obs.obs_info)
    obs_info = Obs.obs_info;
end

if ~isempty(obs_info) && ismember('region_label', obs_info.Properties.VariableNames)
    labels = cellstr(string(obs_info.region_label(:)));
    labels = labels(~cellfun(@isempty, labels));
    if ~isempty(labels)
        target_region_names = unique(labels, 'stable');
        return;
    end
end

model = struct();
if isfield(Obs, 'O') && isstruct(Obs.O) && isfield(Obs.O, 'model') && isstruct(Obs.O.model)
    model = Obs.O.model;
elseif isfield(Obs, 'model') && isstruct(Obs.model)
    model = Obs.model;
end

if isfield(model, 'input_variable_labels') && ~isempty(model.input_variable_labels)
    target_region_names = local_infer_region_names_from_variable_labels(model.input_variable_labels);
end
end


function target_region_names = local_infer_region_names_from_variable_labels(variable_labels)
labels = cellstr(string(variable_labels(:)));
region_names = cell(size(labels));
keep = false(size(labels));

for i = 1:numel(labels)
    tok = regexp(labels{i}, '^(.*)_v\d+(?:_.+)?$', 'tokens', 'once');
    if isempty(tok)
        tok = regexp(labels{i}, '^(.*)_mean(?:_.+)?$', 'tokens', 'once');
    end
    if isempty(tok) || isempty(tok{1})
        continue;
    end
    region_names{i} = tok{1};
    keep(i) = true;
end

if any(keep)
    target_region_names = unique(region_names(keep), 'stable');
else
    target_region_names = {};
end
end


function [coords, region_labels, n_var_region, ana, voxel_variable_labels, ...
        voxel_region_idx, voxel_within_region_idx] = ...
        local_collect_roi_coords_and_anatomy(roiTs, target_region_names)
coords_cells = {};
region_labels = {};
n_var_region = [];
ana = [];
voxel_variable_labels = {};
voxel_region_idx = [];
voxel_within_region_idx = [];
target_region_names = cellstr(string(target_region_names(:)).');

for i_region = 1:size(roiTs, 2)
    R = roiTs{1, i_region};
    if ~isstruct(R) || ~isfield(R, 'coords') || isempty(R.coords)
        continue;
    end
    if ~isfield(R, 'dat') || size(R.dat, 2) ~= size(R.coords, 1)
        continue;
    end
    label = local_roi_label(R, i_region);
    if ~isempty(target_region_names) && ~any(strcmpi(label, target_region_names))
        continue;
    end
    coords_cells{end + 1, 1} = double(R.coords); %#ok<AGROW>
    n_var_region(end + 1, 1) = size(R.coords, 1); %#ok<AGROW>
    region_labels{end + 1, 1} = label; %#ok<AGROW>
    voxel_variable_labels = [voxel_variable_labels; arrayfun(@(k) ... %#ok<AGROW>
        sprintf('%s_v%04d', label, k), (1:size(R.coords, 1)).', 'UniformOutput', false)];
    voxel_region_idx = [voxel_region_idx; repmat(numel(region_labels), size(R.coords, 1), 1)]; %#ok<AGROW>
    voxel_within_region_idx = [voxel_within_region_idx; (1:size(R.coords, 1)).']; %#ok<AGROW>
    if isempty(ana) && isfield(R, 'ana') && ~isempty(R.ana)
        ana = R.ana;
    end
end

if isempty(coords_cells)
    error('No ROI coordinates were found in roiTs for target regions: %s', ...
        strjoin(target_region_names, ', '));
end
if isempty(ana)
    error('No anatomy field .ana was found in roiTs.');
end
coords = cat(1, coords_cells{:});
end


function [labels, origin] = local_resolve_source_variable_labels(Obs, direct_voxel_mode)
labels = {};
origin = '';

if ~direct_voxel_mode
    model = local_resolve_observable_model(Obs);
    if isfield(model, 'input_variable_labels') && ~isempty(model.input_variable_labels)
        labels = cellstr(string(model.input_variable_labels(:)));
        origin = 'model.input_variable_labels';
        return;
    end
end

obs_info = [];
if isfield(Obs, 'observable_info') && istable(Obs.observable_info)
    obs_info = Obs.observable_info;
elseif isfield(Obs, 'obs_info') && istable(Obs.obs_info)
    obs_info = Obs.obs_info;
elseif isfield(Obs, 'O') && isstruct(Obs.O) && ...
        isfield(Obs.O, 'observable_info') && istable(Obs.O.observable_info)
    obs_info = Obs.O.observable_info;
end

if ~isempty(obs_info) && ismember('observable_label', obs_info.Properties.VariableNames)
    labels = cellstr(string(obs_info.observable_label(:)));
    origin = 'observable_info.observable_label';
    return;
end

if isfield(Obs, 'observable_labels') && ~isempty(Obs.observable_labels)
    labels = cellstr(string(Obs.observable_labels(:)));
    origin = 'observable_labels';
    return;
end
if isfield(Obs, 'O') && isstruct(Obs.O) && ...
        isfield(Obs.O, 'observable_labels') && ~isempty(Obs.O.observable_labels)
    labels = cellstr(string(Obs.O.observable_labels(:)));
    origin = 'O.observable_labels';
end
end


function model = local_resolve_observable_model(Obs)
model = struct();
if isfield(Obs, 'O') && isstruct(Obs.O) && ...
        isfield(Obs.O, 'model') && isstruct(Obs.O.model)
    model = Obs.O.model;
elseif isfield(Obs, 'model') && isstruct(Obs.model)
    model = Obs.model;
end
end


function [source_to_voxel_indices, source_space_kind, summary] = ...
        local_build_source_to_voxel_map(source_variable_labels, region_labels, ...
        n_var_region, voxel_variable_labels)
source_variable_labels = cellstr(string(source_variable_labels(:)));
n_source = numel(source_variable_labels);
n_voxel = numel(voxel_variable_labels);

source_to_voxel_indices = cell(n_source, 1);
source_kinds = cell(n_source, 1);
unresolved = {};

voxel_label_to_index = containers.Map('KeyType', 'char', 'ValueType', 'double');
region_to_voxel_indices = containers.Map('KeyType', 'char', 'ValueType', 'any');

offset = 0;
for i_region = 1:numel(region_labels)
    region_name = char(string(region_labels{i_region}));
    idx = offset + (1:n_var_region(i_region));
    region_to_voxel_indices(local_normalize_label_key(region_name)) = idx;
    for j = 1:numel(idx)
        voxel_label_to_index(local_normalize_label_key(voxel_variable_labels{idx(j)})) = idx(j);
    end
    offset = offset + n_var_region(i_region);
end

for i = 1:n_source
    label = char(string(source_variable_labels{i}));
    [target_idx, kind] = local_match_source_label(label, voxel_label_to_index, region_to_voxel_indices);
    source_to_voxel_indices{i} = target_idx;
    source_kinds{i} = kind;
    if isempty(target_idx)
        unresolved{end + 1, 1} = label; %#ok<AGROW>
    end
end

if ~isempty(unresolved)
    preview = unresolved(1:min(6, numel(unresolved)));
    error(['Could not resolve %d activation source label(s) into voxel space. ' ...
        'First labels: %s'], numel(unresolved), strjoin(preview(:).', ', '));
end

if isempty(source_kinds)
    source_space_kind = 'empty';
elseif all(strcmp(source_kinds, 'voxel'))
    source_space_kind = 'voxel';
elseif all(strcmp(source_kinds, 'region_mean'))
    source_space_kind = 'region_mean';
elseif any(strcmp(source_kinds, 'voxel_feature'))
    source_space_kind = 'multi_feature_voxel';
else
    source_space_kind = 'mixed';
end

summary = struct();
summary.n_source_variables = n_source;
summary.n_voxels = n_voxel;
summary.n_region_mean_features = sum(strcmp(source_kinds, 'region_mean'));
summary.n_voxel_features = sum(strcmp(source_kinds, 'voxel_feature'));
summary.n_direct_voxel_features = sum(strcmp(source_kinds, 'voxel'));
summary.source_kinds = source_kinds;
end


function [target_idx, kind] = local_match_source_label(label, voxel_label_to_index, region_to_voxel_indices)
target_idx = [];
kind = '';
key = local_normalize_label_key(label);

if isKey(voxel_label_to_index, key)
    target_idx = voxel_label_to_index(key);
    kind = 'voxel';
    return;
end

tok = regexp(label, '^(.*_v\d+)(?:_.+)?$', 'tokens', 'once');
if ~isempty(tok)
    base_key = local_normalize_label_key(tok{1});
    if isKey(voxel_label_to_index, base_key)
        target_idx = voxel_label_to_index(base_key);
        kind = 'voxel_feature';
        return;
    end
end

tok = regexp(label, '^(.*)_mean(?:_.+)?$', 'tokens', 'once');
if ~isempty(tok)
    region_key = local_normalize_label_key(tok{1});
    if isKey(region_to_voxel_indices, region_key)
        target_idx = region_to_voxel_indices(region_key);
        kind = 'region_mean';
        return;
    end
end
end


function key = local_normalize_label_key(label)
key = lower(strtrim(char(string(label))));
end


function label = local_roi_label(R, i_region)
if ~isfield(R, 'name') || isempty(R.name)
    error('roiTs{1,%d} is missing required ROI label field .name.', i_region);
end
if ischar(R.name) || isstring(R.name)
    label = char(string(R.name));
    return;
end
if iscell(R.name) && isscalar(R.name) && (ischar(R.name{1}) || isstring(R.name{1}))
    label = char(string(R.name{1}));
    return;
end
error('roiTs{1,%d}.name must be a string-like scalar.', i_region);
end


function lim = local_robust_limits(x)
x = x(isfinite(x));
if isempty(x)
    lim = [0 1];
    return;
end
lim = [prctile(x, 1), prctile(x, 99.5)];
if lim(2) <= lim(1)
    lim = [min(x), max(x)];
end
if lim(2) <= lim(1)
    lim = [0 1];
end
end


function mode = local_parse_observable_mode(run_name)
patterns = {'global_slow_band_power_svd100', 'roi_mean_slow_band_power', ...
    'slow_band_power_svd', 'slow_band_power', 'gsvd100_ds', 'global_svd100', 'HP_svd100', ...
    'roi_mean', 'eleHP', 'HP', 'svd', 'abs', 'complex_split', ...
    'complex', 'identity'};
mode = '';
for i = 1:numel(patterns)
    pat = patterns{i};
    if endsWith(run_name, pat, 'IgnoreCase', true) || ...
            contains(run_name, ['_' pat], 'IgnoreCase', true)
        mode = pat;
        return;
    end
end
end


function form = local_parse_residual_form(run_name)
forms = {'projected_kv', 'projected_vlambda'};
form = '';
for i = 1:numel(forms)
    if contains(run_name, forms{i}, 'IgnoreCase', true)
        form = forms{i};
        return;
    end
end
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
