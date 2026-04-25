function [figs, info] = plot_bold_pre_reskoopnet_qc(input_data, params)
%PLOT_BOLD_PRE_RESKOOPNET_QC Visualize BOLD observables before ResKoopNet.
%
%   [figs, info] = plot_bold_pre_reskoopnet_qc(input_data, params)
%
%   input_data can be:
%     1) a BOLD observable MAT file saved by script_batch_build_bold_observables.m
%     2) a struct with top-level observable fields such as obs, session_lengths
%     3) a standardized observable struct O with O.data and session metadata
%
%   The function makes practical QC figures only. It does not preprocess data
%   and does not create new snapshot pairs for training.

if nargin < 2
    params = struct();
end
params = fill_default_params(params);

[B, info] = load_bold_qc_input(input_data);
X = double(B.obs);
if isempty(X)
    error('No observable data found. Expected obs, O.data, or D.data.');
end

if size(X, 1) < 2
    error('Observable data must contain at least two time samples.');
end

[n_samples, n_vars] = size(X);
info.n_samples = n_samples;
info.n_variables = n_vars;
info.output_files = strings(0, 1);

[session_start_idx, session_end_idx, session_lengths] = resolve_session_bounds(B, n_samples);
session_ids = resolve_session_ids(B, numel(session_lengths));
dx = resolve_dx(B, session_lengths);
t = build_time_axis_from_bounds(session_start_idx, session_end_idx, dx, n_samples);

snapshot_info = summarize_snapshot_pairs(B, session_start_idx, session_end_idx);
info.snapshot = snapshot_info;
info.session_ids = session_ids;
info.session_lengths = session_lengths;
info.session_start_idx = session_start_idx;
info.session_end_idx = session_end_idx;
info.dx = dx;

labels = resolve_variable_labels(B, n_vars);
plot_idx = select_variable_indices(n_vars, params.max_variables, params.variable_idx);
heatmap_idx = select_variable_indices(n_vars, params.max_heatmap_variables, params.heatmap_variable_idx);
segment_idx = choose_segment_indices(t, session_start_idx, session_end_idx, dx, params);
thin_idx = choose_thinned_indices(n_samples, params.max_heatmap_samples);

figs = struct();
figs.summary = plot_summary_figure( ...
    X, t, session_ids, session_lengths, session_start_idx, session_end_idx, ...
    dx, labels, plot_idx, heatmap_idx, segment_idx, thin_idx, snapshot_info, ...
    params, info);

if params.save
    info.output_files(end+1, 1) = save_qc_figure(figs.summary, params, 'summary');
end

if params.make_session_gallery
    [figs.session_gallery, gallery_files] = plot_session_gallery( ...
        X, t, session_ids, session_start_idx, session_end_idx, labels, plot_idx, params, info);
    info.output_files = [info.output_files; gallery_files(:)];
end

if params.close_after_save && params.save
    close_all_figs(figs);
end
end

function params = fill_default_params(params)
if ~isfield(params, 'save') || isempty(params.save)
    params.save = false;
end
if ~isfield(params, 'output_dir') || isempty(params.output_dir)
    params.output_dir = fullfile(pwd, 'bold_pre_reskoopnet_qc');
end
if ~isfield(params, 'figure_prefix') || isempty(params.figure_prefix)
    params.figure_prefix = 'bold_pre_reskoopnet_qc';
end
if ~isfield(params, 'visible') || isempty(params.visible)
    params.visible = 'on';
end
if ~isfield(params, 'close_after_save') || isempty(params.close_after_save)
    params.close_after_save = false;
end
if ~isfield(params, 'max_variables') || isempty(params.max_variables)
    params.max_variables = 12;
end
if ~isfield(params, 'max_heatmap_variables') || isempty(params.max_heatmap_variables)
    params.max_heatmap_variables = Inf;
end
if ~isfield(params, 'max_heatmap_samples') || isempty(params.max_heatmap_samples)
    params.max_heatmap_samples = 4000;
end
if ~isfield(params, 'variable_idx')
    params.variable_idx = [];
end
if ~isfield(params, 'heatmap_variable_idx')
    params.heatmap_variable_idx = [];
end
if ~isfield(params, 'segment_session') || isempty(params.segment_session)
    params.segment_session = 'longest';
end
if ~isfield(params, 'segment_start_sec') || isempty(params.segment_start_sec)
    params.segment_start_sec = 0;
end
if ~isfield(params, 'segment_duration_sec') || isempty(params.segment_duration_sec)
    params.segment_duration_sec = 300;
end
if ~isfield(params, 'trace_scale') || isempty(params.trace_scale)
    params.trace_scale = 1;
end
if ~isfield(params, 'trace_clip') || isempty(params.trace_clip)
    params.trace_clip = 4;
end
if ~isfield(params, 'heatmap_clim') || isempty(params.heatmap_clim)
    params.heatmap_clim = [-3, 3];
end
if ~isfield(params, 'make_session_gallery') || isempty(params.make_session_gallery)
    params.make_session_gallery = false;
end
if ~isfield(params, 'session_gallery_max_sessions') || isempty(params.session_gallery_max_sessions)
    params.session_gallery_max_sessions = 6;
end
if ~isfield(params, 'save_formats') || isempty(params.save_formats)
    params.save_formats = {'png'};
end
end

function [B, info] = load_bold_qc_input(input_data)
info = struct();
info.source_file = '';

if ischar(input_data) || isstring(input_data)
    input_file = char(input_data);
    if exist(input_file, 'file') ~= 2
        error('Input file not found: %s', input_file);
    end
    info.source_file = input_file;
    names = {'obs', 'O', 'D', 'obs_info', 'params', 'cfg', ...
        'dx', 'dt', 'fs', 'sampling_period', 'sample_period', 'sampling_frequency', ...
        'session_ids', 'session_lengths', 'session_dx', ...
        'session_start_idx', 'session_end_idx', 'border_idx', ...
        'snapshot_lag_saved', 'snapshot_valid_idx_x', 'snapshot_valid_idx_y', ...
        'snapshot_session_idx', 'snapshot_session_id'};
    B = load_existing_variables(input_file, names);
else
    B = input_data;
end

if isfield(B, 'obs')
    B.obs = B.obs;
elseif isfield(B, 'O') && isstruct(B.O) && isfield(B.O, 'data')
    B.obs = B.O.data;
elseif isfield(B, 'D') && isstruct(B.D) && isfield(B.D, 'data')
    B.obs = B.D.data;
elseif isfield(B, 'data')
    B.obs = B.data;
else
    B.obs = [];
end

if isfield(B, 'O') && isstruct(B.O)
    B = merge_missing_fields(B, B.O);
end
if isfield(B, 'D') && isstruct(B.D)
    B = merge_missing_fields(B, B.D);
end
end

function S = load_existing_variables(input_file, names)
file_info = whos('-file', input_file);
available = {file_info.name};
keep = names(ismember(names, available));
if isempty(keep)
    S = struct();
else
    S = load(input_file, keep{:});
end
end

function A = merge_missing_fields(A, B)
fields = fieldnames(B);
for i = 1:numel(fields)
    key = fields{i};
    if ~isfield(A, key) || isempty(A.(key))
        A.(key) = B.(key);
    end
end
end

function [start_idx, end_idx, lengths] = resolve_session_bounds(B, n_samples)
if isfield(B, 'session_start_idx') && isfield(B, 'session_end_idx') && ...
        ~isempty(B.session_start_idx) && ~isempty(B.session_end_idx)
    start_idx = double(B.session_start_idx(:));
    end_idx = double(B.session_end_idx(:));
    lengths = end_idx - start_idx + 1;
elseif isfield(B, 'session_lengths') && ~isempty(B.session_lengths)
    lengths = double(B.session_lengths(:));
    end_idx = cumsum(lengths);
    start_idx = [1; end_idx(1:end-1) + 1];
else
    lengths = n_samples;
    start_idx = 1;
    end_idx = n_samples;
end

valid = lengths > 0 & start_idx >= 1 & end_idx <= n_samples & start_idx <= end_idx;
start_idx = start_idx(valid);
end_idx = end_idx(valid);
lengths = lengths(valid);

if isempty(lengths)
    error('No valid session bounds were found.');
end
end

function session_ids = resolve_session_ids(B, n_sessions)
if isfield(B, 'session_ids') && ~isempty(B.session_ids)
    session_ids = double(B.session_ids(:));
else
    session_ids = (1:n_sessions).';
end
if numel(session_ids) ~= n_sessions
    session_ids = (1:n_sessions).';
end
end

function dx = resolve_dx(B, session_lengths)
dx = [];
candidate_fields = {'dx', 'dt', 'sampling_period', 'sample_period'};
for i = 1:numel(candidate_fields)
    key = candidate_fields{i};
    if isfield(B, key) && ~isempty(B.(key))
        value = double(B.(key));
        value = value(isfinite(value) & value > 0);
        if ~isempty(value)
            dx = value(1);
            break;
        end
    end
end
if isempty(dx) && isfield(B, 'session_dx') && ~isempty(B.session_dx)
    value = double(B.session_dx(:));
    value = value(isfinite(value) & value > 0);
    if ~isempty(value)
        dx = median(value);
    end
end
if isempty(dx) && isfield(B, 'fs') && ~isempty(B.fs)
    fs = double(B.fs);
    fs = fs(isfinite(fs) & fs > 0);
    if ~isempty(fs)
        dx = 1 / fs(1);
    end
end
if isempty(dx)
    dx = 1;
end
if numel(session_lengths) == 1 %#ok<ISMAT>
    dx = double(dx);
end
end

function t = build_time_axis_from_bounds(start_idx, end_idx, dx, n_samples)
t = ((1:n_samples).' - 1) * dx;
for i = 1:numel(start_idx)
    idx = start_idx(i):end_idx(i);
    t(idx) = (idx - 1) * dx;
end
end

function labels = resolve_variable_labels(B, n_vars)
labels = strings(n_vars, 1);
if isfield(B, 'variable_labels') && numel(B.variable_labels) >= n_vars
    labels = string(B.variable_labels(:));
    labels = labels(1:n_vars);
elseif isfield(B, 'obs_info') && istable(B.obs_info)
    if any(strcmp(B.obs_info.Properties.VariableNames, 'label'))
        labels = string(B.obs_info.label);
    elseif any(strcmp(B.obs_info.Properties.VariableNames, 'variable_label'))
        labels = string(B.obs_info.variable_label);
    end
elseif isfield(B, 'variable_info') && istable(B.variable_info)
    if any(strcmp(B.variable_info.Properties.VariableNames, 'variable_label'))
        labels = string(B.variable_info.variable_label);
    elseif any(strcmp(B.variable_info.Properties.VariableNames, 'region_label'))
        labels = string(B.variable_info.region_label);
    end
end
if numel(labels) < n_vars || all(strlength(labels) == 0)
    labels = compose("var%04d", (1:n_vars).');
else
    labels = labels(1:n_vars);
    empty_mask = strlength(labels) == 0;
    labels(empty_mask) = compose("var%04d", find(empty_mask));
end
end

function idx = select_variable_indices(n_vars, max_vars, requested_idx)
if ~isempty(requested_idx)
    idx = unique(double(requested_idx(:)).', 'stable');
    idx = idx(idx >= 1 & idx <= n_vars);
else
    idx = 1:min(max_vars, n_vars);
end
if isempty(idx)
    idx = 1:min(max_vars, n_vars);
end
end

function idx = choose_segment_indices(t, start_idx, end_idx, dx, params)
session_lengths = end_idx - start_idx + 1;
if isnumeric(params.segment_session)
    i_session = max(1, min(numel(session_lengths), round(params.segment_session)));
elseif strcmpi(params.segment_session, 'first')
    i_session = 1;
elseif strcmpi(params.segment_session, 'last')
    i_session = numel(session_lengths);
else
    [~, i_session] = max(session_lengths);
end

s1 = start_idx(i_session);
s2 = end_idx(i_session);
session_t0 = t(s1);
plot_t1 = session_t0 + max(0, params.segment_start_sec);
plot_t2 = plot_t1 + max(dx, params.segment_duration_sec);
idx = find((1:numel(t)).' >= s1 & (1:numel(t)).' <= s2 & t >= plot_t1 & t <= plot_t2);
if isempty(idx)
    idx = s1:min(s2, s1 + min(session_lengths(i_session), 300) - 1);
end
end

function idx = choose_thinned_indices(n_samples, max_samples)
if n_samples <= max_samples
    idx = 1:n_samples;
else
    idx = unique(round(linspace(1, n_samples, max_samples)));
end
idx = idx(:);
end

function snapshot_info = summarize_snapshot_pairs(B, session_start_idx, session_end_idx)
if isfield(B, 'snapshot_valid_idx_x') && isfield(B, 'snapshot_valid_idx_y') && ...
        ~isempty(B.snapshot_valid_idx_x) && ~isempty(B.snapshot_valid_idx_y)
    ix = double(B.snapshot_valid_idx_x(:));
    iy = double(B.snapshot_valid_idx_y(:));
else
    ix = [];
    iy = [];
    for i = 1:numel(session_start_idx)
        if session_end_idx(i) > session_start_idx(i)
            x_i = (session_start_idx(i):(session_end_idx(i)-1)).';
            ix = [ix; x_i]; %#ok<AGROW>
            iy = [iy; x_i + 1]; %#ok<AGROW>
        end
    end
end

pair_session = zeros(numel(ix), 1);
cross_border = false(numel(ix), 1);
for k = 1:numel(ix)
    sx = find(ix(k) >= session_start_idx & ix(k) <= session_end_idx, 1, 'first');
    sy = find(iy(k) >= session_start_idx & iy(k) <= session_end_idx, 1, 'first');
    if isempty(sx) || isempty(sy) || sx ~= sy
        cross_border(k) = true;
    else
        pair_session(k) = sx;
    end
end

snapshot_info = struct();
snapshot_info.valid_idx_x = ix;
snapshot_info.valid_idx_y = iy;
snapshot_info.n_pairs = numel(ix);
snapshot_info.n_cross_border_pairs = nnz(cross_border);
snapshot_info.lag_values = iy - ix;
snapshot_info.pair_session = pair_session;
snapshot_info.pairs_per_session = zeros(numel(session_start_idx), 1);
valid_pair_session = pair_session(pair_session > 0);
if ~isempty(valid_pair_session)
    snapshot_info.pairs_per_session = accumarray(valid_pair_session, 1, [numel(session_start_idx), 1], @sum, 0);
end
end

function hfig = plot_summary_figure(X, t, session_ids, session_lengths, session_start_idx, session_end_idx, ...
    dx, labels, plot_idx, heatmap_idx, segment_idx, thin_idx, snapshot_info, params, info)

hfig = figure('Color', 'w', 'Position', [80, 80, 1500, 950], 'Visible', params.visible);
tl = tiledlayout(hfig, 3, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

title(tl, compose_title(info), 'Interpreter', 'none', 'FontWeight', 'bold');

ax = nexttile(tl, 1);
bar(ax, session_ids, session_lengths, 0.85, 'FaceColor', [0.2 0.45 0.7]);
box(ax, 'off');
xlabel(ax, 'Session ID');
ylabel(ax, 'Samples');
title(ax, sprintf('Session lengths (TR = %.4g s)', dx));

ax = nexttile(tl, 2, [1, 2]);
plot_offset_traces(ax, X(segment_idx, plot_idx), t(segment_idx), labels(plot_idx), params);
title(ax, sprintf('Representative observable traces (%d of %d)', numel(plot_idx), size(X, 2)));

ax = nexttile(tl, 4, [1, 2]);
Z = robust_zscore(X(thin_idx, heatmap_idx));
imagesc(ax, t(thin_idx), 1:numel(heatmap_idx), Z.');
set(ax, 'YDir', 'normal');
clim(ax, params.heatmap_clim);
colormap(ax, parula(256));
colorbar(ax);
hold(ax, 'on');
plot_session_borders(ax, t, session_end_idx, [0.8 0.8 0.8]);
xlabel(ax, 'Global time (s)');
ylabel(ax, 'Observable');
title(ax, sprintf('Displayed observable heatmap (%d shown of %d total, z-scored, thinned)', ...
    numel(heatmap_idx), size(X, 2)));

ax = nexttile(tl, 6);
std_vals = std(X, 0, 1, 'omitnan');
nan_frac = mean(~isfinite(X), 1);
scatter(ax, std_vals, nan_frac, 12, 'filled', 'MarkerFaceAlpha', 0.5);
box(ax, 'off');
xlabel(ax, 'Variable std');
ylabel(ax, 'NaN/Inf fraction');
title(ax, 'Variable quality summary');

ax = nexttile(tl, 7);
if isempty(snapshot_info.lag_values)
    text(ax, 0.5, 0.5, 'No snapshot pairs found', 'HorizontalAlignment', 'center');
    axis(ax, 'off');
else
    histogram(ax, snapshot_info.lag_values, 'BinMethod', 'integers', 'FaceColor', [0.25 0.55 0.35]);
    box(ax, 'off');
    xlabel(ax, 'idx_y - idx_x');
    ylabel(ax, 'Pair count');
    title(ax, sprintf('Snapshot lags; cross-session pairs = %d', snapshot_info.n_cross_border_pairs));
end

ax = nexttile(tl, 8);
bar(ax, session_ids, snapshot_info.pairs_per_session, 0.85, 'FaceColor', [0.5 0.35 0.7]);
box(ax, 'off');
xlabel(ax, 'Session ID');
ylabel(ax, 'Pairs');
title(ax, 'Snapshot pairs per session');

ax = nexttile(tl, 9);
plot_state_space_preview(ax, X, segment_idx);
title(ax, 'First observables state-space preview');
end

function title_text = compose_title(info)
if isfield(info, 'source_file') && ~isempty(info.source_file)
    [~, name, ext] = fileparts(info.source_file);
    title_text = sprintf('BOLD pre-ResKoopNet QC: %s%s', name, ext);
else
    title_text = 'BOLD pre-ResKoopNet QC';
end
if isfield(info, 'n_samples') && isfield(info, 'n_variables')
    title_text = sprintf('%s | samples=%d, observables=%d', ...
        title_text, info.n_samples, info.n_variables);
end
end

function plot_offset_traces(ax, Xseg, tseg, labels, params)
Xseg = robust_zscore(Xseg);
Xseg = max(min(Xseg, params.trace_clip), -params.trace_clip);
n_plot = size(Xseg, 2);
offsets = (n_plot-1:-1:0) * params.trace_scale * 2.5;
hold(ax, 'on');
for i = 1:n_plot
    plot(ax, tseg, params.trace_scale * Xseg(:, i) + offsets(i), 'k', 'LineWidth', 0.7);
end
[ytick_values, order] = sort(offsets, 'ascend');
set(ax, 'YTick', ytick_values, 'YTickLabel', labels(order));
box(ax, 'off');
xlabel(ax, 'Time (s)');
xlim(ax, [tseg(1), tseg(end)]);
end

function Z = robust_zscore(X)
mu = mean(X, 1, 'omitnan');
sigma = std(X, 0, 1, 'omitnan');
sigma(~isfinite(sigma) | sigma == 0) = 1;
Z = (X - mu) ./ sigma;
Z(~isfinite(Z)) = 0;
end

function plot_session_borders(ax, t, session_end_idx, color_value)
for i = 1:numel(session_end_idx)-1
    idx = session_end_idx(i);
    if idx >= 1 && idx <= numel(t)
        xline(ax, t(idx), '--', 'Color', color_value, 'LineWidth', 0.7);
    end
end
end

function plot_state_space_preview(ax, X, segment_idx)
n_dim = min(3, size(X, 2));
Z = robust_zscore(X(segment_idx, 1:n_dim));
if n_dim >= 3
    plot3(ax, Z(:, 1), Z(:, 2), Z(:, 3), 'Color', [0.1 0.1 0.1], 'LineWidth', 0.7);
    xlabel(ax, 'obs 1');
    ylabel(ax, 'obs 2');
    zlabel(ax, 'obs 3');
    grid(ax, 'on');
else
    plot(ax, Z(:, 1), Z(:, min(2, n_dim)), 'k', 'LineWidth', 0.7);
    xlabel(ax, 'obs 1');
    ylabel(ax, sprintf('obs %d', min(2, n_dim)));
end
box(ax, 'off');
end

function [figs, output_files] = plot_session_gallery(X, t, session_ids, session_start_idx, session_end_idx, labels, plot_idx, params, info)
n_sessions = min(numel(session_ids), params.session_gallery_max_sessions);
figs = gobjects(n_sessions, 1);
output_files = strings(0, 1);
for i = 1:n_sessions
    idx = session_start_idx(i):session_end_idx(i);
    if numel(idx) > params.max_heatmap_samples
        idx = idx(choose_thinned_indices(numel(idx), params.max_heatmap_samples));
    end
    figs(i) = figure('Color', 'w', 'Position', [100, 100, 1200, 520], 'Visible', params.visible);
    tl = tiledlayout(figs(i), 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('%s | session %g', compose_title(info), session_ids(i)), ...
        'Interpreter', 'none', 'FontWeight', 'bold');

    ax = nexttile(tl, 1);
    plot_offset_traces(ax, X(idx, plot_idx), t(idx), labels(plot_idx), params);
    title(ax, 'Session traces');

    ax = nexttile(tl, 2);
    imagesc(ax, t(idx), 1:numel(plot_idx), robust_zscore(X(idx, plot_idx)).');
    set(ax, 'YDir', 'normal');
    clim(ax, params.heatmap_clim);
    colormap(ax, parula(256));
    colorbar(ax);
    xlabel(ax, 'Time (s)');
    ylabel(ax, 'Observable');
    title(ax, 'Session heatmap');

    if params.save
        output_files(end+1, 1) = save_qc_figure(figs(i), params, sprintf('session_%g', session_ids(i))); %#ok<AGROW>
    end
end
end

function output_file = save_qc_figure(hfig, params, tag)
if exist(params.output_dir, 'dir') ~= 7
    mkdir(params.output_dir);
end
safe_tag = regexprep(char(tag), '[^\w.-]+', '_');
base_file = string(fullfile(params.output_dir, sprintf('%s_%s', params.figure_prefix, safe_tag)));
output_file = base_file + ".png";
for i = 1:numel(params.save_formats)
    fmt = lower(char(params.save_formats{i}));
    switch fmt
        case 'png'
            exportgraphics(hfig, base_file + ".png", 'Resolution', 180);
        case 'pdf'
            exportgraphics(hfig, base_file + ".pdf", 'ContentType', 'vector');
        case 'fig'
            savefig(hfig, base_file + ".fig");
        otherwise
            warning('Unsupported save format skipped: %s', fmt);
    end
end
end

function close_all_figs(figs)
fields = fieldnames(figs);
for i = 1:numel(fields)
    h = figs.(fields{i});
    if all(isgraphics(h))
        close(h);
    end
end
end
