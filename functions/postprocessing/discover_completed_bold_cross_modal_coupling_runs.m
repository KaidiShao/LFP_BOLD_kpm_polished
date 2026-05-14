function candidates = discover_completed_bold_cross_modal_coupling_runs(params)
%DISCOVER_COMPLETED_BOLD_CROSS_MODAL_COUPLING_RUNS Discover pipeline 8 candidate runs.

if nargin < 1 || isempty(params)
    params = build_bold_cross_modal_coupling_params();
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

    post_root = io_project.get_pipeline_stage_dir( ...
        params.processed_root, dataset_stem, 7, 'bold_postprocessing');
    if exist(post_root, 'dir') ~= 7
        continue;
    end

    L = dir(fullfile(post_root, '*', 'mat', '*_bold_post.mat'));
    for i_file = 1:numel(L)
        bold_post_file = fullfile(L(i_file).folder, L(i_file).name);
        [~, base_name] = fileparts(bold_post_file);
        run_name = char(erase(string(base_name), "_bold_post"));
        if ~local_matches_exact_filter(run_name, params.run_name_filter)
            continue;
        end
        if ~local_matches_contains_filter(run_name, params.run_name_contains)
            continue;
        end

        observable_mode = local_parse_observable_mode(run_name);
        residual_form = local_parse_residual_form(run_name);
        if ~isempty(params.observable_modes) && ...
                ~any(strcmpi(observable_mode, params.observable_modes))
            continue;
        end
        if ~isempty(params.residual_forms) && ...
                ~any(strcmpi(residual_form, params.residual_forms))
            continue;
        end

        cand = local_candidate_from_bold_post_file(bold_post_file, params.processed_root);
        cand.observable_mode = observable_mode;
        cand.residual_form = residual_form;
        cand.last_write_time = L(i_file).datenum;
        candidates(end + 1, 1) = cand; %#ok<AGROW>
    end
end

if isempty(candidates)
    return;
end

[~, order] = sort(string({candidates.dataset_stem}) + "|" + ...
    string({candidates.observable_mode}) + "|" + ...
    string({candidates.residual_form}) + "|" + ...
    string({candidates.run_name}));
candidates = candidates(order);
end


function params = local_apply_defaults(params)
defaults = build_bold_cross_modal_coupling_params();
params = local_merge_defaults(defaults, params);
params.dataset_stems = cellstr(string(params.dataset_stems(:)).');
params.exclude_dataset_stems = cellstr(string(params.exclude_dataset_stems(:)).');
params.observable_modes = cellstr(string(params.observable_modes(:)).');
params.residual_forms = cellstr(string(params.residual_forms(:)).');
params.run_name_filter = cellstr(string(params.run_name_filter(:)).');
params.run_name_contains = cellstr(string(params.run_name_contains(:)).');
end


function out = local_merge_defaults(defaults, overrides)
out = defaults;
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    value = overrides.(name);
    if isstruct(value) && isfield(defaults, name) && isstruct(defaults.(name))
        out.(name) = local_merge_defaults(defaults.(name), value);
    elseif ~isempty(value)
        out.(name) = value;
    else
        out.(name) = defaults.(name);
    end
end
end


function tf = local_matches_exact_filter(value, filter_values)
if isempty(filter_values)
    tf = true;
else
    tf = any(strcmpi(string(value), string(filter_values)));
end
end


function tf = local_matches_contains_filter(value, patterns)
tf = true;
for i = 1:numel(patterns)
    if ~contains(value, patterns{i}, 'IgnoreCase', true)
        tf = false;
        return;
    end
end
end


function mode = local_parse_observable_mode(run_name)
known = {'global_slow_band_power_svd100', 'roi_mean_slow_band_power', ...
    'slow_band_power_svd', 'gsvd100_ds', 'global_svd100', 'HP_svd100', 'slow_band_power', ...
    'roi_mean', 'eleHP', 'HP', 'svd'};
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


function cand = local_candidate_from_bold_post_file(bold_post_file, processed_root)
if exist(bold_post_file, 'file') ~= 2
    error('bold_post_file does not exist: %s', bold_post_file);
end

[mat_dir, base_name, ~] = fileparts(bold_post_file);
run_root = fileparts(mat_dir);
post_root = fileparts(run_root);
dataset_dir = fileparts(post_root);
[~, dataset_stem] = fileparts(dataset_dir);
run_name = char(erase(string(base_name), "_bold_post"));
save_tag = 'xcorr';

cand = local_empty_candidate();
cand.dataset_stem = dataset_stem;
cand.dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
cand.run_name = run_name;
cand.observable_mode = local_parse_observable_mode(run_name);
cand.residual_form = local_parse_residual_form(run_name);
cand.run_tag = make_bold_cross_modal_run_tag(cand);
cand.bold_post_file = bold_post_file;
cand.xcorr_dir = fullfile( ...
    io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 8, ...
    'efun_density_cross_correlation'), ...
cand.run_tag);
cand.xcorr_file = fullfile(cand.xcorr_dir, [save_tag, '.mat']);
cand.activation_root = fullfile( ...
    io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 8, ...
    'figures_bold_top_xcorr_activation_maps'), ...
cand.run_tag);
end


function cand = local_empty_candidate()
cand = struct('dataset_stem', '', 'dataset_id', '', 'run_name', '', ...
    'observable_mode', '', 'residual_form', '', 'run_tag', '', 'bold_post_file', '', ...
    'xcorr_dir', '', 'xcorr_file', '', 'activation_root', '', ...
    'last_write_time', NaN);
end
