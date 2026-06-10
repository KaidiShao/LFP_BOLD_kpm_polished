% Run peak-by-consensus-state analysis from the current pipeline5 reduction MAT files.
%
% This is the forward path after rerunning pipeline5_eigenfunction_reduction.
% It does not reuse old pipeline5_eigenfunction_peaks_by_state results.
%
% Optional variables before run(...):
%   dataset_stems      = {'e10gb1','e10fV1','e10gh1','f12m01'};
%   peak_mode          = 'max_abs'; % 'max' or 'max_abs'
%   variant_filter     = {'abs_projected_vlambda'};
%   method_tag_filter  = {'svd_k08','mds_k08','umap_k08'};
%   force_recompute    = false;
%   save_figures       = true;
%   peak_activity_suffix = 'rmsenv_adaptive';
%   target_variant_suffix = 'rmsenv_adaptive'; % source reduction is unsuffixed
%   exclude_dataset_variant_pairs = {'e10fV1|complex_split_projected_vlambda'};
%   lfp_activity_transform = 'rms_envelope';
%   lfp_activity_window_policy = 'eigenvalue_adaptive_rms_envelope';
%   continue_on_error  = true;
%   dry_run            = false;
%   max_tasks          = [];

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
set(groot, 'defaultFigureVisible', 'off');

if ~exist('dataset_stems', 'var') || isempty(dataset_stems)
    dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01', 'e10gw1'};
end
if ~exist('peak_mode', 'var') || isempty(peak_mode)
    peak_mode = 'max_abs';
end
if ~exist('variant_filter', 'var')
    variant_filter = {};
end
if ~exist('method_tag_filter', 'var')
    method_tag_filter = {};
end
if ~exist('force_recompute', 'var') || isempty(force_recompute)
    force_recompute = false;
end
if ~exist('save_figures', 'var') || isempty(save_figures)
    save_figures = true;
end
if ~exist('continue_on_error', 'var') || isempty(continue_on_error)
    continue_on_error = true;
end
if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = false;
end
if ~exist('max_tasks', 'var')
    max_tasks = [];
end
if ~exist('target_variant_suffix', 'var')
    target_variant_suffix = '';
end
if ~exist('exclude_dataset_variant_pairs', 'var')
    exclude_dataset_variant_pairs = {};
end

dataset_stems = cellstr(string(dataset_stems(:)).');
peak_mode = char(string(peak_mode));
variant_filter = cellstr(string(variant_filter(:)).');
method_tag_filter = cellstr(string(method_tag_filter(:)).');
exclude_dataset_variant_pairs = cellstr(string(exclude_dataset_variant_pairs(:)).');
force_recompute = logical(force_recompute);
save_figures = logical(save_figures);
continue_on_error = logical(continue_on_error);
dry_run = logical(dry_run);
if ~isempty(max_tasks)
    max_tasks = double(max_tasks);
end

processed_root = io_project.get_project_processed_root();
source_stage_name = io_project.get_pipeline_stage_name(5, 'eigenfunction_reduction');
switch lower(peak_mode)
    case 'max'
        target_stage_name = io_project.get_pipeline_stage_name(5, 'eigenfunction_peaks_by_state');
    case 'max_abs'
        target_stage_name = 'pipeline5_eigenfunction_peaks_by_state_maxabs';
    otherwise
        error('Unsupported peak_mode: %s', peak_mode);
end
if exist('peak_activity_suffix', 'var') && ~isempty(peak_activity_suffix)
    peak_activity_suffix = regexprep(lower(char(string(peak_activity_suffix))), '[^a-z0-9_]+', '_');
    peak_activity_suffix = regexprep(peak_activity_suffix, '^_+|_+$', '');
    target_stage_name = sprintf('pipeline5_eigenfunction_peaks_by_state_%s', ...
        peak_activity_suffix);
else
    peak_activity_suffix = '';
end
target_variant_suffix = char(string(target_variant_suffix));
target_variant_suffix = regexprep(target_variant_suffix, '[^a-zA-Z0-9_]+', '_');
target_variant_suffix = regexprep(target_variant_suffix, '^_+|_+$', '');

activity_opts = struct();
activity_opts.lfp_activity_transform = 'none';
activity_opts.lfp_activity_window_policy = 'none';
activity_opts.value_transform = 'abs';
activity_opts.envelope_enable = false;
activity_opts.envelope_policy = 'none';
activity_opts.envelope_alpha = 0.35;
activity_opts.envelope_min_window_sec = 0.03;
activity_opts.envelope_max_window_sec = 1.0;
activity_opts.envelope_fallback_window_sec = 0.10;
if exist('lfp_activity_transform', 'var') && ~isempty(lfp_activity_transform)
    activity_opts.lfp_activity_transform = char(string(lfp_activity_transform));
end
if exist('lfp_activity_window_policy', 'var') && ~isempty(lfp_activity_window_policy)
    activity_opts.lfp_activity_window_policy = char(string(lfp_activity_window_policy));
end
if exist('density_value_transform', 'var') && ~isempty(density_value_transform)
    activity_opts.value_transform = char(string(density_value_transform));
end
if exist('envelope_enable', 'var') && ~isempty(envelope_enable)
    activity_opts.envelope_enable = logical(envelope_enable);
end
if exist('envelope_policy', 'var') && ~isempty(envelope_policy)
    activity_opts.envelope_policy = char(string(envelope_policy));
end
if exist('envelope_alpha', 'var') && ~isempty(envelope_alpha)
    activity_opts.envelope_alpha = double(envelope_alpha);
end
if exist('envelope_min_window_sec', 'var') && ~isempty(envelope_min_window_sec)
    activity_opts.envelope_min_window_sec = double(envelope_min_window_sec);
end
if exist('envelope_max_window_sec', 'var') && ~isempty(envelope_max_window_sec)
    activity_opts.envelope_max_window_sec = double(envelope_max_window_sec);
end
if exist('envelope_fallback_window_sec', 'var') && ~isempty(envelope_fallback_window_sec)
    activity_opts.envelope_fallback_window_sec = double(envelope_fallback_window_sec);
end

fprintf('Processed root: %s\n', processed_root);
fprintf('Source stage : %s\n', source_stage_name);
fprintf('Target stage : %s\n', target_stage_name);
fprintf('Peak mode    : %s\n\n', peak_mode);

rows = cell(0, 1);
for i_ds = 1:numel(dataset_stems)
    dataset_stem = dataset_stems{i_ds};
    cfg = load_blp_pipeline_cfg_by_stem(dataset_stem);
    source_root = fullfile(processed_root, dataset_stem, source_stage_name);
    if exist(source_root, 'dir') ~= 7
        fprintf('[skip] %s has no %s\n', dataset_stem, source_stage_name);
        continue;
    end

    variants = local_subdirs(source_root);
    for i_var = 1:numel(variants)
        variant_name = variants(i_var).name;
        if ~local_matches_filter(variant_name, variant_filter)
            continue;
        end
        if local_is_excluded_dataset_variant( ...
                dataset_stem, variant_name, exclude_dataset_variant_pairs)
            continue;
        end
        target_variant_name = local_apply_variant_suffix( ...
            variant_name, target_variant_suffix);

        variant_dir = fullfile(variants(i_var).folder, variant_name);
        methods = local_subdirs(variant_dir);
        for i_method = 1:numel(methods)
            method_tag = methods(i_method).name;
            if ~local_matches_filter(method_tag, method_tag_filter)
                continue;
            end

            method_dir = fullfile(methods(i_method).folder, method_tag);
            result_file = local_latest_reduction_result_file(method_dir);
            if isempty(result_file)
                continue;
            end

            save_dir = fullfile(processed_root, dataset_stem, target_stage_name, ...
                target_variant_name, method_tag);
            save_tag = sprintf('%s_%s_%s_peaks', ...
                dataset_stem, target_variant_name, method_tag);
            stats_file = fullfile(save_dir, [save_tag, '_stats.csv']);

            rows{end + 1, 1} = struct( ... %#ok<AGROW>
                'dataset_stem', dataset_stem, ...
                'cfg', cfg, ...
                'variant_name', variant_name, ...
                'target_variant_name', target_variant_name, ...
                'method_tag', method_tag, ...
                'result_file', result_file, ...
                'save_dir', save_dir, ...
                'save_tag', save_tag, ...
                'stats_file', stats_file);
        end
    end
end

fprintf('Discovered %d current pipeline5 peak tasks.\n\n', numel(rows));
if isempty(rows)
    return;
end

if ~isempty(max_tasks) && isfinite(max_tasks) && max_tasks > 0
    rows = rows(1:min(numel(rows), round(max_tasks)));
    fprintf('Limiting to first %d task(s).\n\n', numel(rows));
end

if dry_run
    for i_task = 1:numel(rows)
        task = rows{i_task};
        fprintf('[dry-run %03d/%03d] %s / %s -> %s / %s\n', ...
            i_task, numel(rows), task.dataset_stem, task.variant_name, ...
            task.target_variant_name, task.method_tag);
        fprintf('  result: %s\n', task.result_file);
        fprintf('  stats : %s\n', task.stats_file);
    end
    return;
end

ok_count = 0;
skip_count = 0;
fail_count = 0;
consensus_cache = struct();

for i_task = 1:numel(rows)
    task = rows{i_task};
    fprintf('[%03d/%03d] %s / %s -> %s / %s\n', ...
        i_task, numel(rows), task.dataset_stem, task.variant_name, ...
        task.target_variant_name, task.method_tag);

    if ~force_recompute && exist(task.stats_file, 'file') == 2
        skip_count = skip_count + 1;
        fprintf('  skip existing: %s\n\n', task.stats_file);
        continue;
    end

    try
        [C, source_consensus_file, consensus_cache] = local_load_consensus_cached( ...
            task.dataset_stem, task.cfg, processed_root, consensus_cache);
        S = load_mat_file_with_short_path(task.result_file, 'result');
        if ~isfield(S, 'result')
            error('Missing variable result in %s.', task.result_file);
        end

        params = local_peak_params(task.cfg, peak_mode, save_figures, activity_opts);
        params.save_dir = task.save_dir;
        params.save_tag = task.save_tag;

        analyze_eigenfunction_component_peaks_by_consensus_state( ...
            S.result, C, params, task.result_file, source_consensus_file);

        ok_count = ok_count + 1;
        fprintf('  OK -> %s\n\n', task.stats_file);
    catch ME
        fail_count = fail_count + 1;
        fprintf(2, '  FAILED: %s\n', ME.message);
        if ~continue_on_error
            rethrow(ME);
        end
        fprintf('\n');
    end
end

fprintf('Done. OK: %d | skipped: %d | failed: %d\n', ok_count, skip_count, fail_count);


function dirs = local_subdirs(root_dir)
dirs = dir(root_dir);
dirs = dirs([dirs.isdir]);
dirs = dirs(~ismember({dirs.name}, {'.', '..'}));
[~, order] = sort({dirs.name});
dirs = dirs(order);
end


function tf = local_matches_filter(value, filter_values)
if isempty(filter_values)
    tf = true;
else
    tf = any(strcmpi(char(string(value)), filter_values));
end
end


function tf = local_is_excluded_dataset_variant(dataset_stem, variant_name, pairs)
tf = false;
if isempty(pairs)
    return;
end
key = lower(sprintf('%s|%s', char(string(dataset_stem)), ...
    char(string(variant_name))));
for i = 1:numel(pairs)
    if strcmpi(key, lower(char(string(pairs{i}))))
        tf = true;
        return;
    end
end
end


function out = local_apply_variant_suffix(variant_name, suffix)
variant_name = char(string(variant_name));
suffix = char(string(suffix));
if isempty(suffix)
    out = variant_name;
    return;
end
if endsWith(lower(variant_name), ['_' lower(suffix)])
    out = variant_name;
else
    out = sprintf('%s_%s', variant_name, suffix);
end
end


function result_file = local_latest_reduction_result_file(method_dir)
result_file = '';
mat_dir = fullfile(method_dir, 'mat');
if exist(mat_dir, 'dir') ~= 7
    return;
end
L = dir(fullfile(mat_dir, '*.mat'));
if isempty(L)
    return;
end
names = string({L.name});
keep = contains(names, "_efun__time__", 'IgnoreCase', true) | ...
    contains(names, "_efun__spectrum__", 'IgnoreCase', true);
L = L(keep);
if isempty(L)
    return;
end
[~, idx] = max([L.datenum]);
result_file = fullfile(L(idx).folder, L(idx).name);
end


function [C, source_consensus_file, cache] = local_load_consensus_cached( ...
        dataset_stem, cfg, processed_root, cache)
field_name = matlab.lang.makeValidName(dataset_stem);
if isfield(cache, field_name)
    cached = cache.(field_name);
    C = cached.C;
    source_consensus_file = cached.source_consensus_file;
    return;
end

loader_cfg = struct();
loader_cfg.file_stem = cfg.file_stem;
[C, source_consensus_file] = io_results.load_consensus_state_results( ...
    loader_cfg, processed_root, []);
cache.(field_name) = struct('C', C, 'source_consensus_file', source_consensus_file);
fprintf('  Consensus states: %s\n', source_consensus_file);
end


function params = local_peak_params(cfg, peak_mode, save_figures, activity_opts)
params = struct();
params.component_source = 'smooth_if_available';
params.state_start_idx = local_get_nested(cfg, {'viz', 'state_space_consensus', 'state_start_idx'}, 1);
params.peak_mode = peak_mode;
params.baseline_mode = 'matched';
params.baseline_code = local_get_nested(cfg, {'viz', 'state_space_consensus', 'baseline_code'}, 0);
params.baseline_label = local_get_nested(cfg, {'viz', 'state_space_consensus', 'baseline_label'}, 'baseline');
params.min_window_samples = 1;
params.baseline_random_seed = 1;
params.alpha = 0.05;
params.save_results = true;
params.write_csv = true;
params.save_figures = save_figures;
params.close_figures = true;
params.figure_visible = 'off';
params.lfp_activity_transform = activity_opts.lfp_activity_transform;
params.lfp_activity_window_policy = activity_opts.lfp_activity_window_policy;
params.value_transform = activity_opts.value_transform;
params.envelope_enable = activity_opts.envelope_enable;
params.envelope_policy = activity_opts.envelope_policy;
params.envelope_alpha = activity_opts.envelope_alpha;
params.envelope_min_window_sec = activity_opts.envelope_min_window_sec;
params.envelope_max_window_sec = activity_opts.envelope_max_window_sec;
params.envelope_fallback_window_sec = activity_opts.envelope_fallback_window_sec;
end


function value = local_get_nested(S, names, default_value)
value = S;
for i = 1:numel(names)
    name = names{i};
    if ~isstruct(value) || ~isfield(value, name) || isempty(value.(name))
        value = default_value;
        return;
    end
    value = value.(name);
end
end
