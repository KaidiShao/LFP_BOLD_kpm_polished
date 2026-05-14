function runs = discover_completed_bold_reskoopnet_runs(params)
%DISCOVER_COMPLETED_BOLD_RESKOOPNET_RUNS Discover completed BOLD MLP output runs.

if nargin < 1 || isempty(params)
    params = build_bold_reskoopnet_postprocessing_params();
else
    params = local_apply_defaults(params);
end

runs = repmat(local_empty_run(), 0, 1);
for i_root = 1:numel(params.autodl_roots)
    root_dir = params.autodl_roots{i_root};
    if exist(root_dir, 'dir') ~= 7
        continue;
    end

    dataset_dirs = dir(root_dir);
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

        output_parent = fullfile(root_dir, dataset_stem, 'mlp', 'outputs');
        if exist(output_parent, 'dir') ~= 7
            continue;
        end

        run_dirs = dir(output_parent);
        run_dirs = run_dirs([run_dirs.isdir]);
        run_dirs = run_dirs(~ismember({run_dirs.name}, {'.', '..'}));
        for i_run = 1:numel(run_dirs)
            run_name = run_dirs(i_run).name;
            if ~params.include_smoke && contains(lower(run_name), 'smoke')
                continue;
            end
            if ~local_matches_exact_filter(run_name, params.run_name_filter)
                continue;
            end
            if ~local_matches_contains_filter(run_name, params.run_name_contains)
                continue;
            end

            output_dir = fullfile(run_dirs(i_run).folder, run_name);
            chunk_files = local_collect_output_chunk_files(output_dir, params.filename_pattern);
            if isempty(chunk_files)
                continue;
            end

            summary_files = dir(fullfile(output_dir, '*_summary.mat'));
            if params.require_summary_file && isempty(summary_files)
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

            all_times = [chunk_files.datenum];
            if ~isempty(summary_files)
                all_times = [all_times, [summary_files.datenum]]; %#ok<AGROW>
            end

            run_info = local_empty_run();
            run_info.dataset_stem = dataset_stem;
            run_info.dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
            run_info.run_name = run_name;
            run_info.output_dir = output_dir;
            run_info.autodl_root = root_dir;
            run_info.observable_mode = observable_mode;
            run_info.residual_form = residual_form;
            run_info.n_output_files = numel(chunk_files);
            run_info.n_summary_files = numel(summary_files);
            run_info.last_write_time = max(all_times);
            runs(end + 1, 1) = run_info; %#ok<AGROW>
        end
    end
end

if isempty(runs)
    return;
end

if params.dedupe_by_condition
    runs = local_keep_latest_unique_conditions(runs);
end

[~, order] = sort(string({runs.dataset_stem}) + "|" + ...
    string({runs.observable_mode}) + "|" + string({runs.residual_form}) + "|" + ...
    string({runs.run_name}));
runs = runs(order);
end


function params = local_apply_defaults(params)
defaults = build_bold_reskoopnet_postprocessing_params();
names = fieldnames(defaults);
for i = 1:numel(names)
    name = names{i};
    if ~isfield(params, name) || isempty(params.(name))
        params.(name) = defaults.(name);
    end
end
params.autodl_roots = cellstr(string(params.autodl_roots(:)).');
params.dataset_stems = cellstr(string(params.dataset_stems(:)).');
params.exclude_dataset_stems = cellstr(string(params.exclude_dataset_stems(:)).');
params.observable_modes = cellstr(string(params.observable_modes(:)).');
params.residual_forms = cellstr(string(params.residual_forms(:)).');
params.run_name_filter = cellstr(string(params.run_name_filter(:)).');
params.run_name_contains = cellstr(string(params.run_name_contains(:)).');
end


function runs = local_keep_latest_unique_conditions(runs)
keys = strings(numel(runs), 1);
for i = 1:numel(runs)
    keys(i) = string(lower(strjoin({runs(i).dataset_stem, ...
        runs(i).observable_mode, runs(i).residual_form}, '|')));
end

keep = false(numel(runs), 1);
u = unique(keys, 'stable');
for i = 1:numel(u)
    idx = find(keys == u(i));
    [~, best_local] = max([runs(idx).last_write_time]);
    keep(idx(best_local)) = true;
end
runs = runs(keep);
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


function files = local_collect_output_chunk_files(output_dir, filename_pattern)
L = dir(fullfile(output_dir, filename_pattern));
files = repmat(struct('name', '', 'folder', '', 'fullpath', '', ...
    'chunk_id', NaN, 'datenum', NaN), numel(L), 1);
n_keep = 0;
for i = 1:numel(L)
    name = L(i).name;
    if contains(name, '_outputs_Psi_')
        continue;
    end
    tokens = regexp(name, '^(.*)_outputs_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tokens)
        continue;
    end
    n_keep = n_keep + 1;
    files(n_keep).name = name;
    files(n_keep).folder = L(i).folder;
    files(n_keep).fullpath = fullfile(L(i).folder, name);
    files(n_keep).chunk_id = str2double(tokens{2});
    files(n_keep).datenum = L(i).datenum;
end
files = files(1:n_keep);
if isempty(files)
    return;
end
[~, order] = sort([files.chunk_id]);
files = files(order);
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


function run_info = local_empty_run()
run_info = struct('dataset_stem', '', 'dataset_id', '', 'run_name', '', ...
    'output_dir', '', 'autodl_root', '', 'observable_mode', '', ...
    'residual_form', '', 'n_output_files', 0, 'n_summary_files', 0, ...
    'last_write_time', NaN);
end
