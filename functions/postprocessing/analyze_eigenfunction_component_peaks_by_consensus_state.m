function A = analyze_eigenfunction_component_peaks_by_consensus_state(result, C, params, source_result_file, source_consensus_file)
%ANALYZE_EIGENFUNCTION_COMPONENT_PEAKS_BY_CONSENSUS_STATE
% Quantify temporal-component peaks inside consensus-state windows.
%
% For each consensus-state window and each reduced temporal component, this
% analysis extracts a window peak. It also samples a duration-matched
% baseline window from periods where C.state_code_by_time == 0, then tests
% whether state-specific peak distributions are greater than zero and
% different from their matched baseline distributions.

if nargin < 1 || isempty(result)
    error('result must be provided.');
end

if nargin < 2 || isempty(C)
    error('Consensus-state result C must be provided.');
end

if nargin < 3
    params = struct();
end

if nargin < 4
    source_result_file = '';
end

if nargin < 5
    source_consensus_file = '';
end

params = local_apply_defaults(params, result);
local_validate_inputs(result, C, params);

if isempty(source_consensus_file) && isfield(C, 'save_file')
    source_consensus_file = C.save_file;
end

[components, component_source] = local_get_component_series(result, params);
[state_code_by_traj, state_time_idx] = local_align_state_codes( ...
    C.state_code_by_time(:), size(components, 1), params.state_start_idx);

state_catalog = C.state_catalog(:);
state_codes = double([state_catalog.code]).';
state_labels = string({state_catalog.label}).';

rng(params.baseline_random_seed, 'twister');
baseline_runs = local_find_baseline_runs_by_session( ...
    state_code_by_traj == params.baseline_code, C, params.state_start_idx);

[event_peak_table, baseline_peak_table, window_table] = local_collect_window_peaks( ...
    components, C.state_windows(:), state_catalog, state_code_by_traj, ...
    state_time_idx, baseline_runs, params);

stats_table = local_build_stats_table(event_peak_table, baseline_peak_table, ...
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

if params.save_results || params.write_csv || params.save_figures
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

if params.save_figures
    A.save_paths = local_save_figures(A, params);

    if params.save_results && isfield(A.save_paths, 'main_mat')
        save(A.save_paths.main_mat, 'A', '-v7.3');
    end
end
end


function params = local_apply_defaults(params, result)
if ~isfield(params, 'component_source') || isempty(params.component_source)
    params.component_source = 'smooth_if_available';
end

if ~isfield(params, 'state_start_idx') || isempty(params.state_start_idx)
    params.state_start_idx = 1;
end

if ~isfield(params, 'peak_mode') || isempty(params.peak_mode)
    params.peak_mode = 'max';
end

if ~isfield(params, 'baseline_mode') || isempty(params.baseline_mode)
    params.baseline_mode = 'matched';
end

if ~isfield(params, 'baseline_code') || isempty(params.baseline_code)
    params.baseline_code = 0;
end

if ~isfield(params, 'baseline_label') || isempty(params.baseline_label)
    params.baseline_label = 'baseline';
end

if ~isfield(params, 'min_window_samples') || isempty(params.min_window_samples)
    params.min_window_samples = 1;
end

if ~isfield(params, 'baseline_random_seed') || isempty(params.baseline_random_seed)
    params.baseline_random_seed = 1;
end

if ~isfield(params, 'alpha') || isempty(params.alpha)
    params.alpha = 0.05;
end

if ~isfield(params, 'save_results') || isempty(params.save_results)
    params.save_results = true;
end

if ~isfield(params, 'write_csv') || isempty(params.write_csv)
    params.write_csv = true;
end

if ~isfield(params, 'save_figures') || isempty(params.save_figures)
    params.save_figures = true;
end

if ~isfield(params, 'close_figures') || isempty(params.close_figures)
    params.close_figures = true;
end

if ~isfield(params, 'figure_visible') || isempty(params.figure_visible)
    if isfield(params, 'save_figures') && isequal(params.save_figures, true)
        params.figure_visible = 'off';
    else
        params.figure_visible = 'on';
    end
end

if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = local_default_save_dir(result);
end

if ~isfield(params, 'save_tag') || isempty(params.save_tag)
    params.save_tag = 'efun_peaks';
end

params.peak_mode = char(string(params.peak_mode));
params.baseline_mode = char(string(params.baseline_mode));
params.baseline_label = char(string(params.baseline_label));
end


function local_validate_inputs(result, C, params)
if ~isfield(result, 'core') || ~isfield(result.core, 'temporal_components_time_by_comp')
    error('result.core.temporal_components_time_by_comp is required.');
end

required_fields = {'state_catalog', 'state_windows', 'state_code_by_time'};
for i = 1:numel(required_fields)
    if ~isfield(C, required_fields{i})
        error('Consensus result C is missing field %s.', required_fields{i});
    end
end

if ~strcmpi(params.baseline_mode, 'matched')
    error('Only params.baseline_mode = ''matched'' is currently supported.');
end

if ~ismember(lower(params.peak_mode), {'max', 'max_abs'})
    error('params.peak_mode must be ''max'' or ''max_abs''.');
end

if params.state_start_idx < 1 || params.state_start_idx ~= round(params.state_start_idx)
    error('params.state_start_idx must be a positive integer.');
end
end


function save_dir = local_default_save_dir(result)
if isfield(result, 'cfg') && isfield(result.cfg, 'output') && ...
        isfield(result.cfg.output, 'root') && ~isempty(result.cfg.output.root)
    save_dir = fullfile(result.cfg.output.root, 'peaks');
elseif isfield(result, 'cfg') && isfield(result.cfg, 'save') && ...
        isfield(result.cfg.save, 'dir') && ~isempty(result.cfg.save.dir)
    save_dir = fullfile(fileparts(result.cfg.save.dir), 'peaks');
else
    save_dir = fullfile(pwd, 'peaks');
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


function [state_code_by_traj, state_time_idx] = local_align_state_codes(state_code_by_time, T, state_start_idx)
state_start_idx = double(state_start_idx);
state_time_idx = (state_start_idx:(state_start_idx + T - 1)).';

if state_time_idx(end) > numel(state_code_by_time)
    error(['Consensus state vector length (%d) is shorter than requested ', ...
        'trajectory state index range [%d, %d].'], ...
        numel(state_code_by_time), state_time_idx(1), state_time_idx(end));
end

state_code_by_traj = double(state_code_by_time(state_time_idx));
end


function [event_peak_table, baseline_peak_table, window_table] = local_collect_window_peaks( ...
    components, state_windows, state_catalog, state_code_by_traj, ...
    state_time_idx, baseline_runs, params)
T = size(components, 1);
K = size(components, 2);
state_start_idx = double(params.state_start_idx);
state_end_idx = state_start_idx + T - 1;

state_codes = double([state_catalog.code]);
state_labels = string({state_catalog.label});

win_id = zeros(0, 1);
win_state_code = zeros(0, 1);
win_state_label = strings(0, 1);
win_global_start = zeros(0, 1);
win_global_end = zeros(0, 1);
win_local_start = zeros(0, 1);
win_local_end = zeros(0, 1);
win_duration = zeros(0, 1);
win_source_index = zeros(0, 1);
win_baseline_global_start = zeros(0, 1);
win_baseline_global_end = zeros(0, 1);
win_baseline_local_start = zeros(0, 1);
win_baseline_local_end = zeros(0, 1);

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

out_id = 0;
for i = 1:numel(state_windows)
    if ~isfield(state_windows(i), 'win_global') || isempty(state_windows(i).win_global)
        continue;
    end

    g = double(state_windows(i).win_global(:)).';
    if numel(g) < 2
        continue;
    end

    g1 = max(g(1), state_start_idx);
    g2 = min(g(2), state_end_idx);
    if g2 < g1
        continue;
    end

    local_idx = (g1:g2).' - state_start_idx + 1;
    if numel(local_idx) < params.min_window_samples
        continue;
    end

    code = double(state_windows(i).state_code);
    label = local_label_for_code(code, state_codes, state_labels);
    event_peak = local_window_peak(components(local_idx, :), params.peak_mode);
    baseline_idx = local_sample_matched_baseline_window(baseline_runs, numel(local_idx));

    if isempty(baseline_idx)
        baseline_peak = nan(1, K);
        bg1 = NaN;
        bg2 = NaN;
        bl1 = NaN;
        bl2 = NaN;
    else
        baseline_peak = local_window_peak(components(baseline_idx, :), params.peak_mode);
        bl1 = baseline_idx(1);
        bl2 = baseline_idx(end);
        bg1 = state_time_idx(bl1);
        bg2 = state_time_idx(bl2);
    end

    out_id = out_id + 1;
    win_id(end + 1, 1) = out_id; %#ok<AGROW>
    win_state_code(end + 1, 1) = code; %#ok<AGROW>
    win_state_label(end + 1, 1) = label; %#ok<AGROW>
    win_global_start(end + 1, 1) = g1; %#ok<AGROW>
    win_global_end(end + 1, 1) = g2; %#ok<AGROW>
    win_local_start(end + 1, 1) = local_idx(1); %#ok<AGROW>
    win_local_end(end + 1, 1) = local_idx(end); %#ok<AGROW>
    win_duration(end + 1, 1) = numel(local_idx); %#ok<AGROW>
    win_source_index(end + 1, 1) = i; %#ok<AGROW>
    win_baseline_global_start(end + 1, 1) = bg1; %#ok<AGROW>
    win_baseline_global_end(end + 1, 1) = bg2; %#ok<AGROW>
    win_baseline_local_start(end + 1, 1) = bl1; %#ok<AGROW>
    win_baseline_local_end(end + 1, 1) = bl2; %#ok<AGROW>

    for k = 1:K
        event_window_id(end + 1, 1) = out_id; %#ok<AGROW>
        event_state_code(end + 1, 1) = code; %#ok<AGROW>
        event_state_label(end + 1, 1) = label; %#ok<AGROW>
        event_component_idx(end + 1, 1) = k; %#ok<AGROW>
        event_peak_value(end + 1, 1) = event_peak(k); %#ok<AGROW>
        event_global_start(end + 1, 1) = g1; %#ok<AGROW>
        event_global_end(end + 1, 1) = g2; %#ok<AGROW>
        event_duration(end + 1, 1) = numel(local_idx); %#ok<AGROW>

        baseline_window_id(end + 1, 1) = out_id; %#ok<AGROW>
        baseline_matched_state_code(end + 1, 1) = code; %#ok<AGROW>
        baseline_matched_state_label(end + 1, 1) = label; %#ok<AGROW>
        baseline_component_idx(end + 1, 1) = k; %#ok<AGROW>
        baseline_peak_value(end + 1, 1) = baseline_peak(k); %#ok<AGROW>
        baseline_global_start(end + 1, 1) = bg1; %#ok<AGROW>
        baseline_global_end(end + 1, 1) = bg2; %#ok<AGROW>
        baseline_duration(end + 1, 1) = numel(local_idx); %#ok<AGROW>
    end
end

window_table = table(win_id, win_source_index, win_state_code, win_state_label, ...
    win_global_start, win_global_end, win_local_start, win_local_end, win_duration, ...
    win_baseline_global_start, win_baseline_global_end, ...
    win_baseline_local_start, win_baseline_local_end, ...
    'VariableNames', {'window_id', 'source_state_window_idx', 'state_code', ...
    'state_label', 'global_start_idx', 'global_end_idx', 'local_start_idx', ...
    'local_end_idx', 'duration_samples', 'baseline_global_start_idx', ...
    'baseline_global_end_idx', 'baseline_local_start_idx', 'baseline_local_end_idx'});

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

if ~isempty(event_peak_table)
    invalid_state_samples = state_code_by_traj(event_peak_table.global_start_idx - state_start_idx + 1) == params.baseline_code;
    if any(invalid_state_samples)
        warning('Some event windows start on baseline samples after clipping to the trajectory range.');
    end
end
end


function label = local_label_for_code(code, state_codes, state_labels)
idx = find(state_codes == code, 1, 'first');
if isempty(idx)
    label = string(sprintf('state_%d', code));
else
    label = state_labels(idx);
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


function baseline_idx = local_sample_matched_baseline_window(baseline_runs, n_samples)
baseline_idx = [];
if isempty(baseline_runs)
    return;
end

run_lengths = baseline_runs(:, 2) - baseline_runs(:, 1) + 1;
eligible = find(run_lengths >= n_samples);
if isempty(eligible)
    return;
end

run_idx = eligible(randi(numel(eligible)));
start_min = baseline_runs(run_idx, 1);
start_max = baseline_runs(run_idx, 2) - n_samples + 1;
start_idx = randi([start_min, start_max]);
baseline_idx = (start_idx:(start_idx + n_samples - 1)).';
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


function runs = local_find_baseline_runs_by_session(mask, C, state_start_idx)
mask = logical(mask(:));
T = numel(mask);
state_start_idx = double(state_start_idx);
runs = zeros(0, 2);

if isfield(C, 'session_start_idx') && isfield(C, 'session_end_idx') && ...
        ~isempty(C.session_start_idx) && ~isempty(C.session_end_idx)
    session_start_idx = double(C.session_start_idx(:));
    session_end_idx = double(C.session_end_idx(:));

    for i = 1:numel(session_start_idx)
        local_start = max(1, session_start_idx(i) - state_start_idx + 1);
        local_end = min(T, session_end_idx(i) - state_start_idx + 1);
        if local_end < local_start
            continue;
        end

        session_runs = local_find_true_runs(mask(local_start:local_end));
        if ~isempty(session_runs)
            session_runs = session_runs + local_start - 1;
            runs = [runs; session_runs]; %#ok<AGROW>
        end
    end
else
    runs = local_find_true_runs(mask);
end
end


function runs = local_find_true_runs(mask)
mask = logical(mask(:));
if isempty(mask)
    runs = zeros(0, 2);
    return;
end

d = diff([false; mask; false]);
starts = find(d == 1);
ends = find(d == -1) - 1;
runs = [starts, ends];
end


function save_paths = local_save_figures(A, params)
save_paths = A.save_paths;
if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

fig1 = local_plot_peak_distributions(A, params);
save_paths.peak_distribution_png = fullfile(params.save_dir, ...
    [params.save_tag, '_dist.png']);
exportgraphics(fig1, save_paths.peak_distribution_png, 'Resolution', 220);
if params.close_figures
    close(fig1);
end

fig2 = local_plot_mean_peak_heatmap(A, params);
save_paths.mean_peak_heatmap_png = fullfile(params.save_dir, ...
    [params.save_tag, '_mean.png']);
exportgraphics(fig2, save_paths.mean_peak_heatmap_png, 'Resolution', 220);
if params.close_figures
    close(fig2);
end

fig3 = local_plot_baseline_effect_heatmap(A, params);
save_paths.baseline_effect_heatmap_png = fullfile(params.save_dir, ...
    [params.save_tag, '_effect.png']);
exportgraphics(fig3, save_paths.baseline_effect_heatmap_png, 'Resolution', 220);
if params.close_figures
    close(fig3);
end
end


function fig = local_plot_peak_distributions(A, params)
K = A.meta.n_components;
state_codes = A.state_codes(:).';
state_labels = A.state_labels(:).';

fig = figure('Color', 'w', 'Position', [100, 80, 1280, 820], ...
    'Visible', params.figure_visible, ...
    'Name', 'Temporal component peak distributions', 'NumberTitle', 'off');
tiled = tiledlayout(fig, ceil(K / 2), min(K, 2), ...
    'TileSpacing', 'compact', 'Padding', 'compact');

for k = 1:K
    ax = nexttile(tiled);
    [values, groups] = local_distribution_values_for_component(A, k, state_codes, state_labels, params);

    if isempty(values)
        text(ax, 0.5, 0.5, sprintf('Component %d: no data', k), ...
            'HorizontalAlignment', 'center');
        axis(ax, 'off');
        continue;
    end

    if exist('boxchart', 'file') == 2
        boxchart(ax, categorical(groups), values, 'BoxFaceAlpha', 0.55);
    else
        boxplot(ax, values, groups, 'LabelOrientation', 'inline');
    end
    yline(ax, 0, 'k--', 'LineWidth', 0.8);
    title(ax, sprintf('Component %d', k), 'Interpreter', 'none');
    ylabel(ax, sprintf('Peak (%s)', A.meta.peak_mode), 'Interpreter', 'none');
    grid(ax, 'on');
    set(ax, 'FontSize', 10);
end

title(tiled, sprintf('Temporal Component Peaks By Consensus State (%s components)', ...
    A.meta.component_source), 'Interpreter', 'none');
end


function [values, groups] = local_distribution_values_for_component(A, comp_idx, state_codes, state_labels, params)
values = zeros(0, 1);
groups = strings(0, 1);

base_mask = A.baseline_peak_table.component_idx == comp_idx;
base_values = A.baseline_peak_table.peak_value(base_mask);
base_values = base_values(isfinite(base_values));
values = [values; base_values];
groups = [groups; repmat(string(sprintf('%d %s', params.baseline_code, params.baseline_label)), numel(base_values), 1)];

for s = 1:numel(state_codes)
    code = state_codes(s);
    label = state_labels(s);
    mask = A.event_peak_table.component_idx == comp_idx & ...
        A.event_peak_table.state_code == code;
    state_values = A.event_peak_table.peak_value(mask);
    state_values = state_values(isfinite(state_values));
    values = [values; state_values]; %#ok<AGROW>
    groups = [groups; repmat(string(sprintf('%d %s', code, label)), numel(state_values), 1)]; %#ok<AGROW>
end
end


function fig = local_plot_mean_peak_heatmap(A, params)
[matrix, y_labels, x_labels] = local_stats_matrix(A.stats_table, ...
    'event_mean_peak', A.state_codes, A.state_labels, A.meta.n_components);

fig = figure('Color', 'w', 'Position', [120, 90, 980, 520], ...
    'Visible', params.figure_visible, ...
    'Name', 'Mean component peak heatmap', 'NumberTitle', 'off');
ax = axes(fig);
imagesc(ax, matrix);
colormap(ax, turbo(256));
colorbar(ax);
set(ax, 'XTick', 1:numel(x_labels), 'XTickLabel', x_labels, ...
    'YTick', 1:numel(y_labels), 'YTickLabel', y_labels, ...
    'TickLabelInterpreter', 'none');
xlabel(ax, 'Temporal component');
ylabel(ax, 'Consensus state');
title(ax, sprintf('Mean Window Peak (%s)', params.peak_mode), 'Interpreter', 'none');
local_annotate_heatmap(ax, matrix, []);
end


function fig = local_plot_baseline_effect_heatmap(A, params)
[effect_matrix, y_labels, x_labels] = local_stats_matrix(A.stats_table, ...
    'mean_event_minus_baseline', A.state_codes, A.state_labels, A.meta.n_components);
[q_matrix, ~, ~] = local_stats_matrix(A.stats_table, ...
    'q_vs_baseline_paired_ttest_two_sided', A.state_codes, A.state_labels, A.meta.n_components);

fig = figure('Color', 'w', 'Position', [120, 90, 1040, 540], ...
    'Visible', params.figure_visible, ...
    'Name', 'Peak effect vs matched baseline', 'NumberTitle', 'off');
ax = axes(fig);
imagesc(ax, effect_matrix);
colormap(ax, local_blue_white_red(256));
max_abs_effect = max(abs(effect_matrix(isfinite(effect_matrix))));
if ~isempty(max_abs_effect) && max_abs_effect > 0
    clim(ax, [-max_abs_effect, max_abs_effect]);
end
colorbar(ax);
set(ax, 'XTick', 1:numel(x_labels), 'XTickLabel', x_labels, ...
    'YTick', 1:numel(y_labels), 'YTickLabel', y_labels, ...
    'TickLabelInterpreter', 'none');
xlabel(ax, 'Temporal component');
ylabel(ax, 'Consensus state');
title(ax, 'Mean Peak Difference: Consensus Window - Matched Baseline', ...
    'Interpreter', 'none');
local_annotate_heatmap(ax, effect_matrix, q_matrix);
end


function [matrix, y_labels, x_labels] = local_stats_matrix(stats_table, field_name, state_codes, state_labels, K)
matrix = nan(numel(state_codes), K);
y_labels = strings(numel(state_codes), 1);
x_labels = strings(1, K);

for s = 1:numel(state_codes)
    y_labels(s) = sprintf('%d %s', state_codes(s), state_labels(s));
    for k = 1:K
        x_labels(k) = sprintf('Comp %d', k);
        mask = stats_table.state_code == state_codes(s) & stats_table.component_idx == k;
        if any(mask)
            matrix(s, k) = stats_table.(field_name)(find(mask, 1, 'first'));
        end
    end
end
end


function local_annotate_heatmap(ax, matrix, q_matrix)
for r = 1:size(matrix, 1)
    for c = 1:size(matrix, 2)
        value = matrix(r, c);
        if ~isfinite(value)
            label = 'n/a';
        elseif nargin >= 3 && ~isempty(q_matrix)
            label = sprintf('%.3g%s', value, local_stars(q_matrix(r, c)));
        else
            label = sprintf('%.3g', value);
        end
        text(ax, c, r, label, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 9, ...
            'Color', 'k');
    end
end
end


function stars = local_stars(q)
if ~isfinite(q)
    stars = '';
elseif q < 0.001
    stars = '***';
elseif q < 0.01
    stars = '**';
elseif q < 0.05
    stars = '*';
else
    stars = '';
end
end


function cmap = local_blue_white_red(n)
if nargin < 1
    n = 256;
end

x = linspace(0, 1, n).';
cmap = zeros(n, 3);
for i = 1:n
    if x(i) < 0.5
        t = x(i) / 0.5;
        cmap(i, :) = [t, t, 1];
    else
        t = (x(i) - 0.5) / 0.5;
        cmap(i, :) = [1, 1 - t, 1 - t];
    end
end
end
