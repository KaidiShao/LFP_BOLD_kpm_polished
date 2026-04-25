function out = compute_spike_residual_lagged_correlation(result_input, params)
%COMPUTE_SPIKE_RESIDUAL_LAGGED_CORRELATION Lagged spike-vs-residual correlation.
%
% Positive lag means spike leads residual:
%   corr_lag(k) = corr(spike(t), residual(t + k))
%
% Inputs
%   result_input : either a saved comparison MAT file path containing
%                  variable "result", or an in-memory result struct.
%   params       : optional struct
%
% Key params fields
%   .features       default {'abs_rms', 'abs_mean'}
%   .max_lag_bins   default 80   (2 s when spike dx = 0.025 s)
%   .top_n_rows     default 200
%   .save_results   default true
%   .save_dir       default same folder as result_input file
%   .save_tag       default derived from source result file name
%   .verbose        default true

if nargin < 1 || isempty(result_input)
    error('result_input must be provided.');
end

if nargin < 2
    params = struct();
end

[base_result, base_result_file] = local_load_base_result(result_input);
params = local_apply_defaults(params, base_result_file);

if params.verbose
    fprintf('Lagged-correlation source:\n  %s\n', base_result_file);
end

required_session_fields = { ...
    'spike_rate', ...
    'residual_abs_mean', ...
    'residual_abs_rms', ...
    'residual_real_mean', ...
    'residual_imag_mean'};
local_validate_base_result(base_result, required_session_fields);

session_spike_dx = arrayfun(@(s) double(s.spike_dx), base_result.sessions(:));
if any(abs(session_spike_dx - session_spike_dx(1)) > 1e-12)
    error('Spike dx is inconsistent across sessions in the saved result.');
end

lag_bins = (-params.max_lag_bins:params.max_lag_bins).';
lag_sec = lag_bins * session_spike_dx(1);
n_lags = numel(lag_bins);
n_channels = height(base_result.channel_table);
n_modes = height(base_result.mode_table);

lag_corr = struct();
peak_tables = cell(numel(params.features), 1);
site_tables = cell(numel(params.features), 1);

for i_feature = 1:numel(params.features)
    feature_name = char(string(params.features{i_feature}));
    if params.verbose
        fprintf('Computing lagged correlation for feature %s ...\n', feature_name);
    end

    [Xcat, Ycat] = local_build_pooled_with_gaps(base_result.sessions, feature_name, params.max_lag_bins);
    corr_cube = nan(n_channels, n_modes, n_lags);

    for i_lag = 1:n_lags
        lag_val = lag_bins(i_lag);
        [Xlag, Ylag] = local_align_with_lag(Xcat, Ycat, lag_val);
        corr_cube(:, :, i_lag) = local_corr_matrix(Xlag, Ylag);
    end

    lag_corr.(feature_name) = corr_cube;

    [peak_map, peak_table] = local_summarize_peak_map( ...
        corr_cube, lag_bins, lag_sec, feature_name, base_result.channel_table, base_result.mode_table);
    lag_corr.([feature_name, '_peak']) = peak_map;
    peak_tables{i_feature} = peak_table;
    site_tables{i_feature} = local_build_site_summary(peak_table);
end

peak_table = vertcat(peak_tables{:});
if ~isempty(peak_table)
    finite_peak = isfinite(peak_table.peak_abs_corr);
    peak_table = [ ...
        sortrows(peak_table(finite_peak, :), {'peak_abs_corr', 'feature'}, {'descend', 'ascend'}); ...
        peak_table(~finite_peak, :)];
end

site_summary = vertcat(site_tables{:});
top_peak_table = peak_table(1:min(height(peak_table), params.top_n_rows), :);

out = struct();
out.base_result_file = base_result_file;
out.params = params;
out.lag_bins = lag_bins;
out.lag_sec = lag_sec;
out.positive_lag_definition = 'spike leads residual';
out.features = cellstr(string(params.features(:)));
out.channel_table = base_result.channel_table;
out.mode_table = base_result.mode_table;
out.source = local_build_source_summary(base_result);
out.lag_corr = lag_corr;
out.peak_table = peak_table;
out.top_peak_table = top_peak_table;
out.site_summary = site_summary;

if params.save_results
    save_paths = local_save_out(out, params, base_result_file);
else
    save_paths = struct('main_mat', '', 'peak_csv', '', 'top_peak_csv', '', 'site_summary_csv', '');
end
out.save_paths = save_paths;
end


function [base_result, base_result_file] = local_load_base_result(result_input)
if ischar(result_input) || isstring(result_input)
    base_result_file = char(string(result_input));
    if exist(base_result_file, 'file') ~= 2
        error('Result file does not exist: %s', base_result_file);
    end

    S = load(base_result_file, 'result');
    if ~isfield(S, 'result')
        error('Result file %s does not contain variable "result".', base_result_file);
    end
    base_result = S.result;
elseif isstruct(result_input)
    base_result = result_input;
    base_result_file = '';
else
    error('result_input must be a path or a result struct.');
end
end


function params = local_apply_defaults(params, base_result_file)
if ~isfield(params, 'features') || isempty(params.features)
    params.features = {'abs_rms', 'abs_mean'};
end
params.features = cellstr(string(params.features(:)).');

if ~isfield(params, 'max_lag_bins') || isempty(params.max_lag_bins)
    params.max_lag_bins = 80;
end
if ~isscalar(params.max_lag_bins) || params.max_lag_bins < 1 || params.max_lag_bins ~= floor(params.max_lag_bins)
    error('params.max_lag_bins must be a positive integer.');
end

if ~isfield(params, 'top_n_rows') || isempty(params.top_n_rows)
    params.top_n_rows = 200;
end
if ~isfield(params, 'save_results')
    params.save_results = true;
end
if ~isfield(params, 'verbose')
    params.verbose = true;
end

if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    if ~isempty(base_result_file)
        params.save_dir = fileparts(base_result_file);
    else
        params.save_dir = pwd;
    end
end

if ~isfield(params, 'save_tag') || isempty(params.save_tag)
    if ~isempty(base_result_file)
        [~, stem] = fileparts(base_result_file);
        params.save_tag = [stem, '_lagged'];
    else
        params.save_tag = 'spike_residual_lagged';
    end
end
end


function local_validate_base_result(base_result, required_session_fields)
if ~isfield(base_result, 'sessions') || isempty(base_result.sessions)
    error('Base result does not contain any saved sessions.');
end
if ~isfield(base_result, 'channel_table') || isempty(base_result.channel_table)
    error('Base result does not contain channel_table.');
end
if ~isfield(base_result, 'mode_table') || isempty(base_result.mode_table)
    error('Base result does not contain mode_table.');
end

for i = 1:numel(required_session_fields)
    name = required_session_fields{i};
    if ~isfield(base_result.sessions(1), name)
        error('Base result sessions do not contain field %s.', name);
    end
end
end


function [Xcat, Ycat] = local_build_pooled_with_gaps(session_structs, feature_name, gap_len)
n_sessions = numel(session_structs);
n_channels = size(session_structs(1).spike_rate, 2);

residual_field = ['residual_', feature_name];
if ~isfield(session_structs(1), residual_field)
    error('Saved result sessions do not contain field %s.', residual_field);
end
n_modes = size(session_structs(1).(residual_field), 2);

session_lengths = arrayfun(@(s) size(s.spike_rate, 1), session_structs(:));
total_len = sum(session_lengths) + gap_len * max(0, n_sessions - 1);

Xcat = nan(total_len, n_channels);
Ycat = nan(total_len, n_modes);

cursor = 1;
for i = 1:n_sessions
    Xi = double(session_structs(i).spike_rate);
    Yi = double(session_structs(i).(residual_field));
    n_i = size(Xi, 1);

    if size(Yi, 1) ~= n_i
        error('Session %d has mismatched spike and residual lengths.', i);
    end

    idx = cursor:(cursor + n_i - 1);
    Xcat(idx, :) = Xi;
    Ycat(idx, :) = Yi;
    cursor = idx(end) + 1;

    if i < n_sessions
        cursor = cursor + gap_len;
    end
end
end


function [Xlag, Ylag] = local_align_with_lag(Xcat, Ycat, lag_val)
n = size(Xcat, 1);

if lag_val > 0
    Xlag = Xcat(1:(n - lag_val), :);
    Ylag = Ycat((1 + lag_val):n, :);
elseif lag_val < 0
    lag_abs = -lag_val;
    Xlag = Xcat((1 + lag_abs):n, :);
    Ylag = Ycat(1:(n - lag_abs), :);
else
    Xlag = Xcat;
    Ylag = Ycat;
end
end


function C = local_corr_matrix(X, Y)
X = double(X);
Y = double(Y);

if isempty(X) || isempty(Y)
    C = zeros(size(X, 2), size(Y, 2));
    return;
end

valid = all(isfinite(X), 2) & all(isfinite(Y), 2);
X = X(valid, :);
Y = Y(valid, :);

if size(X, 1) < 2
    C = nan(size(X, 2), size(Y, 2));
    return;
end

X = X - mean(X, 1);
Y = Y - mean(Y, 1);

sx = std(X, 0, 1);
sy = std(Y, 0, 1);
sx(sx == 0) = NaN;
sy(sy == 0) = NaN;

C = (X.' * Y) ./ ((size(X, 1) - 1) .* (sx(:) * sy(:).'));
end


function [peak_map, peak_table] = local_summarize_peak_map(corr_cube, lag_bins, lag_sec, feature_name, channel_table, mode_table)
[n_channels, n_modes, ~] = size(corr_cube);

peak_corr = nan(n_channels, n_modes);
peak_abs_corr = nan(n_channels, n_modes);
peak_lag_bins = nan(n_channels, n_modes);
peak_lag_sec = nan(n_channels, n_modes);
zero_lag_idx = find(lag_bins == 0, 1, 'first');
zero_corr = corr_cube(:, :, zero_lag_idx);

rows = {};
for i_ch = 1:n_channels
    for i_mode = 1:n_modes
        curve = squeeze(corr_cube(i_ch, i_mode, :));
        if all(~isfinite(curve))
            continue;
        end

        [~, idx_peak] = max(abs(curve));
        peak_corr(i_ch, i_mode) = curve(idx_peak);
        peak_abs_corr(i_ch, i_mode) = abs(curve(idx_peak));
        peak_lag_bins(i_ch, i_mode) = lag_bins(idx_peak);
        peak_lag_sec(i_ch, i_mode) = lag_sec(idx_peak);

        rows(end + 1, :) = { ... %#ok<AGROW>
            feature_name, ...
            channel_table.channel_index(i_ch), ...
            channel_table.channel_site(i_ch), ...
            channel_table.channel_label(i_ch), ...
            mode_table.mode_rank(i_mode), ...
            mode_table.mode_original_idx(i_mode), ...
            mode_table.evalue_real(i_mode), ...
            mode_table.evalue_imag(i_mode), ...
            mode_table.evalue_abs(i_mode), ...
            zero_corr(i_ch, i_mode), ...
            peak_corr(i_ch, i_mode), ...
            peak_abs_corr(i_ch, i_mode), ...
            peak_lag_bins(i_ch, i_mode), ...
            peak_lag_sec(i_ch, i_mode), ...
            peak_abs_corr(i_ch, i_mode) - abs(zero_corr(i_ch, i_mode))};
    end
end

peak_map = struct();
peak_map.peak_corr = peak_corr;
peak_map.peak_abs_corr = peak_abs_corr;
peak_map.peak_lag_bins = peak_lag_bins;
peak_map.peak_lag_sec = peak_lag_sec;
peak_map.zero_corr = zero_corr;

if isempty(rows)
    peak_table = table();
else
    peak_table = cell2table(rows, 'VariableNames', { ...
        'feature', ...
        'channel_index', 'channel_site', 'channel_label', ...
        'mode_rank', 'mode_original_idx', ...
        'evalue_real', 'evalue_imag', 'evalue_abs', ...
        'zero_corr', 'peak_corr', 'peak_abs_corr', ...
        'peak_lag_bins', 'peak_lag_sec', 'delta_abs_vs_zero'});
end
end


function T = local_build_site_summary(peak_table)
if isempty(peak_table)
    T = table();
    return;
end

sites = unique(string(peak_table.channel_site), 'stable');
features = unique(string(peak_table.feature), 'stable');
rows = {};

for i_feature = 1:numel(features)
    for i_site = 1:numel(sites)
        mask = string(peak_table.feature) == features(i_feature) & string(peak_table.channel_site) == sites(i_site);
        Ts = peak_table(mask, :);
        if isempty(Ts)
            continue;
        end

        finite_mask = isfinite(Ts.peak_abs_corr);
        Ts = Ts(finite_mask, :);
        if isempty(Ts)
            continue;
        end

        [max_val, idx_max] = max(Ts.peak_abs_corr);
        rows(end + 1, :) = { ... %#ok<AGROW>
            char(features(i_feature)), ...
            char(sites(i_site)), ...
            mean(Ts.peak_abs_corr), ...
            max_val, ...
            Ts.channel_label(idx_max), ...
            Ts.mode_rank(idx_max), ...
            Ts.peak_lag_bins(idx_max), ...
            Ts.peak_lag_sec(idx_max)};
    end
end

if isempty(rows)
    T = table();
else
    T = cell2table(rows, 'VariableNames', { ...
        'feature', 'channel_site', 'mean_peak_abs_corr', 'max_peak_abs_corr', ...
        'channel_label_at_max', 'mode_rank_at_max', 'peak_lag_bins_at_max', 'peak_lag_sec_at_max'});
end
end


function save_paths = local_save_out(out, params, base_result_file)
if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

tag = local_filename_safe(params.save_tag);

save_paths = struct();
save_paths.main_mat = fullfile(params.save_dir, [tag, '.mat']);
save_paths.peak_csv = fullfile(params.save_dir, [tag, '_peak.csv']);
save_paths.top_peak_csv = fullfile(params.save_dir, [tag, '_top_peak.csv']);
save_paths.site_summary_csv = fullfile(params.save_dir, [tag, '_site_summary.csv']);

save(save_paths.main_mat, 'out', '-v7.3');
if ~isempty(out.peak_table)
    writetable(out.peak_table, save_paths.peak_csv);
    writetable(out.top_peak_table, save_paths.top_peak_csv);
else
    writetable(table(), save_paths.peak_csv);
    writetable(table(), save_paths.top_peak_csv);
end
if ~isempty(out.site_summary)
    writetable(out.site_summary, save_paths.site_summary_csv);
else
    writetable(table(), save_paths.site_summary_csv);
end
end


function out = local_filename_safe(name_in)
out = regexprep(char(name_in), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled';
end
end


function source = local_build_source_summary(base_result)
source = struct();
source.n_sessions = numel(base_result.sessions);
source.n_channels = height(base_result.channel_table);
source.n_modes = height(base_result.mode_table);

if isfield(base_result, 'source')
    if isfield(base_result.source, 'outputs_dir')
        source.outputs_dir = base_result.source.outputs_dir;
    end
    if isfield(base_result.source, 'save_paths')
        source.base_save_paths = base_result.source.save_paths;
    end
end

if isfield(base_result, 'lfp_selected_channels')
    source.lfp_selected_channels = base_result.lfp_selected_channels;
end
if isfield(base_result, 'spike_selected_channels')
    source.spike_selected_channels = base_result.spike_selected_channels;
end
end
