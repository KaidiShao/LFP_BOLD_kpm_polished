% Recompute pipeline5 peak-by-state outputs using peak_mode = 'max_abs'
% across every dataset that already has pipeline5_eigenfunction_peaks_by_state.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
set(groot, 'defaultFigureVisible', 'off');

processed_root = io_project.get_project_processed_root();
source_stage_name = 'pipeline5_eigenfunction_peaks_by_state';
target_stage_name = 'pipeline5_eigenfunction_peaks_by_state_maxabs';

fprintf('Processed root: %s\n', processed_root);
fprintf('Source stage : %s\n', source_stage_name);
fprintf('Target stage : %s\n\n', target_stage_name);

dataset_dirs = dir(processed_root);
dataset_dirs = dataset_dirs([dataset_dirs.isdir]);
dataset_dirs = dataset_dirs(~ismember({dataset_dirs.name}, {'.', '..'}));

task_rows = cell(0, 1);
for i = 1:numel(dataset_dirs)
    dataset_name = dataset_dirs(i).name;
    source_stage_dir = fullfile(processed_root, dataset_name, source_stage_name);
    if exist(source_stage_dir, 'dir') ~= 7
        continue;
    end

    variant_dirs = dir(source_stage_dir);
    variant_dirs = variant_dirs([variant_dirs.isdir]);
    variant_dirs = variant_dirs(~ismember({variant_dirs.name}, {'.', '..'}));
    for j = 1:numel(variant_dirs)
        variant_name = variant_dirs(j).name;
        variant_dir = fullfile(source_stage_dir, variant_name);

        method_dirs = dir(variant_dir);
        method_dirs = method_dirs([method_dirs.isdir]);
        method_dirs = method_dirs(~ismember({method_dirs.name}, {'.', '..'}));
        for k = 1:numel(method_dirs)
            method_name = method_dirs(k).name;
            method_dir = fullfile(variant_dir, method_name);
            mat_list = dir(fullfile(method_dir, '*_peaks.mat'));
            if isempty(mat_list)
                continue;
            end

            old_mat_file = fullfile(mat_list(1).folder, mat_list(1).name);
            new_save_dir = fullfile(processed_root, dataset_name, target_stage_name, variant_name, method_name);
            task_rows{end + 1, 1} = struct( ... %#ok<AGROW>
                'dataset', dataset_name, ...
                'variant', variant_name, ...
                'method', method_name, ...
                'old_mat_file', old_mat_file, ...
                'new_save_dir', new_save_dir);
        end
    end
end

fprintf('Discovered %d peak-analysis tasks.\n\n', numel(task_rows));
if isempty(task_rows)
    return;
end

ok_count = 0;
fail_count = 0;
for t = 1:numel(task_rows)
    task = task_rows{t};
    fprintf('[%03d/%03d] %s / %s / %s\n', t, numel(task_rows), ...
        task.dataset, task.variant, task.method);
    try
        local_rerun_one( ...
            task.old_mat_file, task.new_save_dir, processed_root, ...
            task.dataset, task.variant, task.method);
        ok_count = ok_count + 1;
        fprintf('  OK -> %s\n\n', task.new_save_dir);
    catch ME
        fail_count = fail_count + 1;
        fprintf(2, '  FAILED: %s\n', ME.message);
        for s = 1:numel(ME.stack)
            fprintf(2, '    at %s:%d\n', ME.stack(s).name, ME.stack(s).line);
        end
        fprintf('\n');
    end
end

fprintf('Done. Success: %d | Failed: %d\n', ok_count, fail_count);


function local_rerun_one(old_mat_file, new_save_dir, processed_root, dataset_name, variant_name, method_name)
S = load(old_mat_file, 'A');
if ~isfield(S, 'A')
    error('File %s does not contain variable A.', old_mat_file);
end

A_old = S.A;
if ~isfield(A_old, 'source') || ...
        ~isfield(A_old.source, 'result_file') || ...
        ~isfield(A_old.source, 'consensus_file')
    error('File %s is missing source.result_file or source.consensus_file.', old_mat_file);
end

result_file = local_resolve_source_result_file( ...
    char(string(A_old.source.result_file)), ...
    processed_root, dataset_name, variant_name, method_name);
result = local_load_result_struct(result_file);

params = A_old.params;
params.peak_mode = 'max_abs';
params.save_dir = new_save_dir;
params.save_results = true;
params.write_csv = true;
params.save_figures = false;
params.close_figures = true;
params.figure_visible = 'off';

if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

consensus_file = '';
try
    consensus_file = local_resolve_dataset_file( ...
        char(string(A_old.source.consensus_file)), ...
        processed_root, dataset_name);
catch ME
    fprintf('  Consensus file unavailable for direct rerun: %s\n', ME.message);
end

if ~isempty(consensus_file)
    try
        C = local_load_consensus_struct(consensus_file);
        analyze_eigenfunction_component_peaks_by_consensus_state( ...
            result, C, params, result_file, consensus_file);
        return;
    catch ME
        if ~contains(string(ME.message), "Consensus state vector length")
            rethrow(ME);
        end
        fprintf('  Falling back to stored window table due to consensus length mismatch.\n');
    end
else
    fprintf('  Falling back to stored window table because consensus file could not be resolved.\n');
end

local_rerun_from_saved_windows(A_old, result, params, result_file, consensus_file);
end


function resolved_path = local_resolve_source_result_file( ...
        original_path, processed_root, dataset_name, variant_name, method_name)
if exist(original_path, 'file') == 2
    resolved_path = original_path;
    return;
end

[~, base_name, ext] = fileparts(original_path);
file_name = [base_name ext];

candidate_pattern = fullfile( ...
    processed_root, dataset_name, 'pipeline5_eigenfunction_reduction', ...
    variant_name, sprintf('%s_*', method_name), 'mat', file_name);
candidates = dir(candidate_pattern);
if numel(candidates) == 1
    resolved_path = fullfile(candidates(1).folder, candidates(1).name);
    fprintf('  Remapped result file -> %s\n', resolved_path);
    return;
elseif numel(candidates) > 1
    resolved_path = local_select_best_candidate(candidates, file_name);
    fprintf('  Remapped result file (best match) -> %s\n', resolved_path);
    return;
end

resolved_path = local_resolve_dataset_file(original_path, processed_root, dataset_name);
if exist(resolved_path, 'file') == 2
    fprintf('  Resolved result file by dataset search -> %s\n', resolved_path);
    return;
end

error('Missing source result file: %s', original_path);
end


function resolved_path = local_resolve_dataset_file(original_path, processed_root, dataset_name)
if exist(original_path, 'file') == 2
    resolved_path = original_path;
    return;
end

[~, base_name, ext] = fileparts(original_path);
file_name = [base_name ext];
dataset_root = fullfile(processed_root, dataset_name);
candidates = dir(fullfile(dataset_root, '**', file_name));
if isempty(candidates)
    error('Missing source file: %s', original_path);
end

resolved_path = local_select_best_candidate(candidates, file_name);
fprintf('  Resolved missing file -> %s\n', resolved_path);
end


function resolved_path = local_select_best_candidate(candidates, file_name)
if isempty(candidates)
    error('No candidates found for %s.', file_name);
end

candidate_paths = strings(numel(candidates), 1);
scores = zeros(numel(candidates), 1);
for i = 1:numel(candidates)
    candidate_paths(i) = string(fullfile(candidates(i).folder, candidates(i).name));
    score = 0;
    if contains(candidate_paths(i), "pipeline5_eigenfunction_reduction", 'IgnoreCase', true)
        score = score + 100;
    end
    if contains(candidate_paths(i), "\mat\", 'IgnoreCase', true)
        score = score + 10;
    end
    score = score + candidates(i).datenum / 1e10;
    scores(i) = score;
end

[~, best_idx] = max(scores);
resolved_path = char(candidate_paths(best_idx));
end


function local_rerun_from_saved_windows(A_old, result, params, source_result_file, source_consensus_file)
if ~isfield(A_old, 'window_table') || ~istable(A_old.window_table)
    error('Stored analysis is missing A.window_table required for fallback rerun.');
end
if ~isfield(A_old, 'state_catalog') || isempty(A_old.state_catalog)
    error('Stored analysis is missing A.state_catalog required for fallback rerun.');
end

[components, component_source] = local_get_component_series(result, params);
window_table = A_old.window_table;
[event_peak_table, baseline_peak_table] = local_recompute_peaks_from_window_table( ...
    components, window_table, params.peak_mode);

state_catalog = A_old.state_catalog(:);
state_codes = double([state_catalog.code]).';
state_labels = string({state_catalog.label}).';
stats_table = local_build_stats_table( ...
    event_peak_table, baseline_peak_table, ...
    state_codes, state_labels, size(components, 2), params);

A = struct();
A.meta = struct();
A.meta.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
A.meta.analysis_name = 'eigenfunction_component_peaks_by_consensus_state';
A.meta.component_source = component_source;
A.meta.n_time_samples = size(components, 1);
A.meta.n_components = size(components, 2);
A.meta.state_start_idx = double(params.state_start_idx);
A.meta.state_end_idx = double(params.state_start_idx) + size(components, 1) - 1;
A.meta.peak_mode = params.peak_mode;
A.meta.baseline_mode = params.baseline_mode;
A.meta.baseline_code = double(params.baseline_code);
A.meta.baseline_label = params.baseline_label;
A.meta.window_source = 'stored_window_table';

A.source = struct();
A.source.result_file = source_result_file;
A.source.consensus_file = source_consensus_file;

A.params = params;
A.state_catalog = state_catalog;
A.state_codes = state_codes;
A.state_labels = state_labels;
A.window_table = window_table;
A.event_peak_table = event_peak_table;
A.baseline_peak_table = baseline_peak_table;
A.stats_table = stats_table;
A.save_paths = struct();

local_write_analysis_outputs(A, params);
end


function [event_peak_table, baseline_peak_table] = local_recompute_peaks_from_window_table( ...
        components, window_table, peak_mode)
T = size(components, 1);
K = size(components, 2);

event_window_id = zeros(0, 1);
event_state_code = zeros(0, 1);
event_state_label = strings(0, 1);
event_component_idx = zeros(0, 1);
event_peak_value = zeros(0, 1);
event_global_start = zeros(0, 1);
event_global_end = zeros(0, 1);
event_duration = zeros(0, 1);

baseline_window_id = zeros(0, 1);
baseline_matched_state_code = zeros(0, 1);
baseline_matched_state_label = strings(0, 1);
baseline_component_idx = zeros(0, 1);
baseline_peak_value = zeros(0, 1);
baseline_global_start = zeros(0, 1);
baseline_global_end = zeros(0, 1);
baseline_duration = zeros(0, 1);

for i = 1:height(window_table)
    win_id = double(window_table.window_id(i));
    code = double(window_table.state_code(i));
    label = string(window_table.state_label(i));

    local_start = double(window_table.local_start_idx(i));
    local_end = double(window_table.local_end_idx(i));
    local_validate_window_bounds(local_start, local_end, T, ...
        sprintf('event window_id=%d', win_id));
    local_idx = (local_start:local_end).';
    event_peak = local_window_peak(components(local_idx, :), peak_mode);

    baseline_start = double(window_table.baseline_local_start_idx(i));
    baseline_end = double(window_table.baseline_local_end_idx(i));
    if isfinite(baseline_start) && isfinite(baseline_end)
        local_validate_window_bounds(baseline_start, baseline_end, T, ...
            sprintf('baseline window_id=%d', win_id));
        baseline_idx = (baseline_start:baseline_end).';
        baseline_peak = local_window_peak(components(baseline_idx, :), peak_mode);
        baseline_n = numel(baseline_idx);
    else
        baseline_peak = nan(1, K);
        baseline_n = NaN;
    end

    event_n = numel(local_idx);
    event_g1 = double(window_table.global_start_idx(i));
    event_g2 = double(window_table.global_end_idx(i));
    baseline_g1 = double(window_table.baseline_global_start_idx(i));
    baseline_g2 = double(window_table.baseline_global_end_idx(i));

    for k = 1:K
        event_window_id(end + 1, 1) = win_id; %#ok<AGROW>
        event_state_code(end + 1, 1) = code; %#ok<AGROW>
        event_state_label(end + 1, 1) = label; %#ok<AGROW>
        event_component_idx(end + 1, 1) = k; %#ok<AGROW>
        event_peak_value(end + 1, 1) = event_peak(k); %#ok<AGROW>
        event_global_start(end + 1, 1) = event_g1; %#ok<AGROW>
        event_global_end(end + 1, 1) = event_g2; %#ok<AGROW>
        event_duration(end + 1, 1) = event_n; %#ok<AGROW>

        baseline_window_id(end + 1, 1) = win_id; %#ok<AGROW>
        baseline_matched_state_code(end + 1, 1) = code; %#ok<AGROW>
        baseline_matched_state_label(end + 1, 1) = label; %#ok<AGROW>
        baseline_component_idx(end + 1, 1) = k; %#ok<AGROW>
        baseline_peak_value(end + 1, 1) = baseline_peak(k); %#ok<AGROW>
        baseline_global_start(end + 1, 1) = baseline_g1; %#ok<AGROW>
        baseline_global_end(end + 1, 1) = baseline_g2; %#ok<AGROW>
        baseline_duration(end + 1, 1) = baseline_n; %#ok<AGROW>
    end
end

event_peak_table = table(event_window_id, event_state_code, event_state_label, ...
    event_component_idx, event_peak_value, event_global_start, event_global_end, ...
    event_duration, ...
    'VariableNames', {'window_id', 'state_code', 'state_label', ...
    'component_idx', 'peak_value', 'global_start_idx', 'global_end_idx', ...
    'duration_samples'});

baseline_peak_table = table(baseline_window_id, baseline_matched_state_code, ...
    baseline_matched_state_label, baseline_component_idx, baseline_peak_value, ...
    baseline_global_start, baseline_global_end, baseline_duration, ...
    'VariableNames', {'matched_window_id', 'matched_state_code', ...
    'matched_state_label', 'component_idx', 'peak_value', ...
    'global_start_idx', 'global_end_idx', 'duration_samples'});
end


function local_validate_window_bounds(start_idx, end_idx, T, label)
if ~isfinite(start_idx) || ~isfinite(end_idx)
    error('Non-finite indices for %s.', label);
end
if start_idx < 1 || end_idx > T || end_idx < start_idx
    error('Invalid index range [%d, %d] for %s with T=%d.', ...
        start_idx, end_idx, label, T);
end
end


function local_write_analysis_outputs(A, params)
if params.save_results || params.write_csv
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end
end

if params.save_results
    A.save_paths.main_mat = fullfile(params.save_dir, [params.save_tag, '.mat']);
    save(A.save_paths.main_mat, 'A', '-v7.3');
end

if params.write_csv
    A.save_paths.window_csv = fullfile(params.save_dir, [params.save_tag, '_win.csv']);
    A.save_paths.event_peak_csv = fullfile(params.save_dir, [params.save_tag, '_event.csv']);
    A.save_paths.baseline_peak_csv = fullfile(params.save_dir, [params.save_tag, '_base.csv']);
    A.save_paths.stats_csv = fullfile(params.save_dir, [params.save_tag, '_stats.csv']);

    writetable(A.window_table, A.save_paths.window_csv);
    writetable(A.event_peak_table, A.save_paths.event_peak_csv);
    writetable(A.baseline_peak_table, A.save_paths.baseline_peak_csv);
    writetable(A.stats_table, A.save_paths.stats_csv);
end
end


function [components, source_name] = local_get_component_series(result, params)
components = result.core.temporal_components_time_by_comp;
source_name = 'raw';

has_smooth = isfield(result, 'summary') && ...
    isfield(result.summary, 'temporal_components_smooth_time_by_comp') && ...
    ~isempty(result.summary.temporal_components_smooth_time_by_comp);

switch lower(params.component_source)
    case 'raw'
        source_name = 'raw';

    case 'smooth'
        if ~has_smooth
            error('No smoothed temporal components are available in result.summary.');
        end
        components = result.summary.temporal_components_smooth_time_by_comp;
        source_name = 'smooth';

    case 'smooth_if_available'
        if has_smooth
            components = result.summary.temporal_components_smooth_time_by_comp;
            source_name = 'smooth';
        end

    otherwise
        error('Unknown params.component_source = %s.', params.component_source);
end
end


function peak = local_window_peak(X, peak_mode)
peak = nan(1, size(X, 2));
for k = 1:size(X, 2)
    x = X(:, k);
    x = x(isfinite(x));
    if isempty(x)
        continue;
    end

    switch lower(peak_mode)
        case 'max'
            peak(k) = max(x);
        case 'max_abs'
            peak(k) = max(abs(x));
    end
end
end


function stats_table = local_build_stats_table(event_peak_table, baseline_peak_table, ...
        state_codes, state_labels, n_components, params)
rows = cell(0, 1);

for s = 1:numel(state_codes)
    code = state_codes(s);
    label = state_labels(s);

    for k = 1:n_components
        event_mask = event_peak_table.state_code == code & ...
            event_peak_table.component_idx == k;
        base_mask = baseline_peak_table.matched_state_code == code & ...
            baseline_peak_table.component_idx == k;

        x = event_peak_table.peak_value(event_mask);
        b = baseline_peak_table.peak_value(base_mask);
        pair_ok = isfinite(x) & isfinite(b);
        x_pair = x(pair_ok);
        b_pair = b(pair_ok);
        x = x(isfinite(x));
        b = b(isfinite(b));

        [p_zero, t_zero, df_zero] = local_one_sample_ttest(x, 0, 'right');
        [p_base, t_base, df_base] = local_paired_ttest(x_pair, b_pair);
        p_zero_signrank = local_signrank_right(x, 0);
        p_base_signrank = local_signrank_two_sided(x_pair, b_pair);

        row = table(code, label, k, numel(x), numel(b), numel(x_pair), ...
            local_nanmean(x), local_nanmedian(x), local_nanstd(x), local_sem(x), ...
            local_nanmean(b), local_nanmedian(b), local_nanstd(b), local_sem(b), ...
            local_nanmean(x_pair - b_pair), local_nanmedian(x_pair - b_pair), ...
            local_cohen_d_paired(x_pair, b_pair), ...
            p_zero, t_zero, df_zero, p_zero_signrank, ...
            p_base, t_base, df_base, p_base_signrank, ...
            'VariableNames', {'state_code', 'state_label', 'component_idx', ...
            'n_event_windows', 'n_baseline_windows', 'n_paired_windows', ...
            'event_mean_peak', 'event_median_peak', 'event_std_peak', 'event_sem_peak', ...
            'baseline_mean_peak', 'baseline_median_peak', 'baseline_std_peak', ...
            'baseline_sem_peak', 'mean_event_minus_baseline', ...
            'median_event_minus_baseline', 'cohen_d_paired_vs_baseline', ...
            'p_gt_zero_ttest_right', 't_gt_zero', 'df_gt_zero', ...
            'p_gt_zero_signrank_right', 'p_vs_baseline_paired_ttest_two_sided', ...
            't_vs_baseline', 'df_vs_baseline', ...
            'p_vs_baseline_signrank_two_sided'});
        rows{end + 1, 1} = row; %#ok<AGROW>
    end
end

if isempty(rows)
    stats_table = table();
else
    stats_table = vertcat(rows{:});
    stats_table.q_gt_zero_ttest_right = local_bh_fdr(stats_table.p_gt_zero_ttest_right);
    stats_table.q_vs_baseline_paired_ttest_two_sided = ...
        local_bh_fdr(stats_table.p_vs_baseline_paired_ttest_two_sided);
    stats_table.significant_gt_zero = stats_table.q_gt_zero_ttest_right < params.alpha;
    stats_table.significant_vs_baseline = ...
        stats_table.q_vs_baseline_paired_ttest_two_sided < params.alpha;
end
end


function [p, tstat, df] = local_one_sample_ttest(x, mu, tail)
x = x(isfinite(x));
n = numel(x);
df = n - 1;
tstat = NaN;
p = NaN;

if n < 2
    return;
end

sd = std(x);
if sd == 0
    tstat = Inf * sign(mean(x) - mu);
    p = local_degenerate_p(mean(x) - mu, tail);
    return;
end

tstat = (mean(x) - mu) / (sd / sqrt(n));
p = local_t_pvalue(tstat, df, tail);
end


function [p, tstat, df] = local_paired_ttest(x, y)
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
[p, tstat, df] = local_one_sample_ttest(x - y, 0, 'both');
end


function p = local_t_pvalue(tstat, df, tail)
if exist('tcdf', 'file') ~= 2 || ~isfinite(tstat) || df < 1
    p = NaN;
    return;
end

switch lower(tail)
    case 'right'
        p = 1 - tcdf(tstat, df);
    case 'left'
        p = tcdf(tstat, df);
    case 'both'
        p = 2 * (1 - tcdf(abs(tstat), df));
    otherwise
        error('Unknown tail = %s.', tail);
end

p = min(max(p, 0), 1);
end


function p = local_degenerate_p(delta, tail)
switch lower(tail)
    case 'right'
        p = double(delta <= 0);
    case 'left'
        p = double(delta >= 0);
    case 'both'
        p = double(delta == 0);
    otherwise
        error('Unknown tail = %s.', tail);
end
end


function p = local_signrank_right(x, mu)
x = x(isfinite(x));
if numel(x) < 2 || exist('signrank', 'file') ~= 2
    p = NaN;
    return;
end

try
    p = signrank(x, mu, 'tail', 'right');
catch
    p = NaN;
end
end


function p = local_signrank_two_sided(x, y)
ok = isfinite(x) & isfinite(y);
x = x(ok);
y = y(ok);
if numel(x) < 2 || exist('signrank', 'file') ~= 2
    p = NaN;
    return;
end

try
    p = signrank(x, y);
catch
    p = NaN;
end
end


function value = local_nanmean(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = mean(x);
end
end


function value = local_nanmedian(x)
x = x(isfinite(x));
if isempty(x)
    value = NaN;
else
    value = median(x);
end
end


function value = local_nanstd(x)
x = x(isfinite(x));
if numel(x) < 2
    value = NaN;
else
    value = std(x);
end
end


function value = local_sem(x)
x = x(isfinite(x));
if numel(x) < 2
    value = NaN;
else
    value = std(x) / sqrt(numel(x));
end
end


function d = local_cohen_d_paired(x, y)
ok = isfinite(x) & isfinite(y);
delta = x(ok) - y(ok);
if numel(delta) < 2 || std(delta) == 0
    d = NaN;
else
    d = mean(delta) / std(delta);
end
end


function q = local_bh_fdr(p)
q = nan(size(p));
valid = isfinite(p);
pv = p(valid);
if isempty(pv)
    return;
end

[ps, order] = sort(pv(:));
m = numel(ps);
qs = ps .* m ./ (1:m).';
qs = flipud(cummin(flipud(qs)));
qs = min(qs, 1);

qv = nan(size(pv(:)));
qv(order) = qs;
q(valid) = qv;
end


function result = local_load_result_struct(result_file)
S = load(result_file);
if isfield(S, 'result')
    result = S.result;
    return;
end

fields = fieldnames(S);
for i = 1:numel(fields)
    value = S.(fields{i});
    if isstruct(value) && isfield(value, 'core') && isfield(value.core, 'temporal_components_time_by_comp')
        result = value;
        return;
    end
end

error('Could not find a valid result struct in %s.', result_file);
end


function C = local_load_consensus_struct(consensus_file)
S = load(consensus_file);
if isfield(S, 'C')
    C = S.C;
    return;
end

fields = fieldnames(S);
for i = 1:numel(fields)
    value = S.(fields{i});
    if isstruct(value) && isfield(value, 'state_catalog') && ...
            isfield(value, 'state_windows') && isfield(value, 'state_code_by_time')
        C = value;
        return;
    end
end

error('Could not find a valid consensus struct in %s.', consensus_file);
end
