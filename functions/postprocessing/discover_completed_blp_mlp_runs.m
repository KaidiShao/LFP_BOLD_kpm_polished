function runs = discover_completed_blp_mlp_runs(params)
%DISCOVER_COMPLETED_BLP_MLP_RUNS Find completed per-condition MLP outputs.

runs = struct('dataset', {}, 'observable_mode', {}, 'residual_form', {}, ...
    'run_name', {}, 'output_dir', {}, 'summary_count', {}, 'chunk_count', {}, ...
    'last_write_time', {});

for i = 1:numel(params.dataset_stems)
    dataset = params.dataset_stems{i};
    if any(strcmpi(dataset, params.exclude_dataset_stems))
        continue;
    end

    output_parent = fullfile(params.autodl_root, dataset, 'mlp', 'outputs');
    if exist(output_parent, 'dir') ~= 7
        continue;
    end

    L = dir(output_parent);
    L = L([L.isdir]);
    L = L(~ismember({L.name}, {'.', '..'}));

    for k = 1:numel(L)
        run_name = L(k).name;
        if contains(lower(run_name), 'smoke')
            continue;
        end

        output_dir = fullfile(L(k).folder, run_name);
        summary_files = dir(fullfile(output_dir, '*_summary.mat'));
        chunk_files = local_collect_chunk_files(output_dir, params.filename_pattern);
        allow_summary_only = isfield(params, 'allow_summary_only_outputs') && ...
            logical(params.allow_summary_only_outputs);
        if isempty(summary_files)
            continue;
        end
        if isempty(chunk_files) && ...
                ~(allow_summary_only && local_has_summary_edmd_outputs(summary_files))
            continue;
        end

        [observable_mode, residual_form] = local_infer_condition_from_run_name(run_name);
        file_times = [summary_files.datenum, chunk_files.datenum];
        last_write = max(file_times);

        run = struct();
        run.dataset = dataset;
        run.observable_mode = observable_mode;
        run.residual_form = residual_form;
        run.run_name = run_name;
        run.output_dir = output_dir;
        run.summary_count = numel(summary_files);
        run.chunk_count = numel(chunk_files);
        run.last_write_time = last_write;
        runs(end+1) = run; %#ok<AGROW>
    end
end

if isempty(runs)
    return;
end

keys = strings(numel(runs), 1);
for i = 1:numel(runs)
    keys(i) = string(lower(strjoin({runs(i).dataset, ...
        runs(i).observable_mode, runs(i).residual_form}, '|')));
end

keep = false(numel(runs), 1);
u = unique(keys, 'stable');
for i = 1:numel(u)
    idx = find(keys == u(i)); %#ok<FNDSB>
    [~, best_local] = max([runs(idx).last_write_time]);
    keep(idx(best_local)) = true;
end
runs = runs(keep);

runs = filter_invalid_blp_mlp_runs(runs, params);
if isempty(runs)
    return;
end

runs = local_filter_discovered_runs(runs, params);
if isempty(runs)
    return;
end

[~, order] = sort(string({runs.dataset}) + "|" + string({runs.observable_mode}) + "|" + string({runs.residual_form}));
runs = runs(order);
end


function tf = local_has_summary_edmd_outputs(summary_files)
tf = false;
if isempty(summary_files)
    return;
end

for i_file = 1:numel(summary_files)
    file_name = fullfile(summary_files(i_file).folder, summary_files(i_file).name);
    try
        vars = who('-file', file_name);
    catch
        vars = {};
    end
    if any(strcmp(vars, 'EDMD_outputs'))
        tf = true;
        return;
    end
end
end


function runs = local_filter_discovered_runs(runs, params)
if isempty(runs)
    return;
end

if isfield(params, 'run_name_filter') && ~isempty(params.run_name_filter)
    keep = false(numel(runs), 1);
    for i = 1:numel(runs)
        keep(i) = local_matches_filter(runs(i).run_name, params.run_name_filter);
    end
    runs = runs(keep);
end

if isempty(runs)
    return;
end

if isfield(params, 'condition_key_filter') && ~isempty(params.condition_key_filter)
    keep = false(numel(runs), 1);
    for i = 1:numel(runs)
        keep(i) = local_matches_filter(local_condition_key(runs(i)), ...
            params.condition_key_filter);
    end
    runs = runs(keep);
end
end


function tf = local_matches_filter(value, filter_values)
tf = any(strcmpi(string(value), string(filter_values)));
end


function key = local_condition_key(run_info)
key = lower(strjoin({run_info.dataset, run_info.observable_mode, ...
    run_info.residual_form}, '|'));
end


function [observable_mode, residual_form] = local_infer_condition_from_run_name(run_name)
name_lower = lower(run_name);
if contains(name_lower, 'complex_split') || contains(name_lower, 'complexsplit')
    observable_mode = 'complex_split';
elseif contains(name_lower, '_abs') || endsWith(name_lower, 'abs')
    observable_mode = 'abs';
else
    observable_mode = 'unknown';
end

if contains(name_lower, 'projected_vlambda')
    residual_form = 'projected_vlambda';
elseif contains(name_lower, 'projected_kv')
    residual_form = 'projected_kv';
else
    residual_form = 'unknown';
end
end


function files = local_collect_chunk_files(output_dir, filename_pattern)
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
