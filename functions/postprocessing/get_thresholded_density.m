function [D, fig] = get_thresholded_density(phi, dt, params)
%GET_THRESHOLDED_DENSITY Compute above-threshold eigenfunction density.
%
%   D = get_thresholded_density(phi, dt, params)
%
% Inputs
%   phi    [T x K] eigenfunction time series.
%   dt     Sampling interval in seconds. Required when window/step sizes are
%          specified in seconds.
%   params Struct with optional fields:
%          window_sec, step_sec, window_samples, step_samples
%          threshold_ratio, threshold_mode
%          value_transform: 'auto' | 'none' | 'abs' | 'real'
%          smooth_density, smooth_window_sec, smooth_window_bins
%          save_results, save_dir, save_stem, save_tag
%          make_figure, save_figure, figure_dir
%
% Output
%   D.density_time_by_mode is [Nwin x K], the fraction of samples in each
%   window where the transformed eigenfunction is above its threshold.

if nargin < 2
    dt = [];
end

if nargin < 3 || isempty(params)
    params = struct();
end

params = local_apply_defaults(params);
local_validate_phi(phi);

if ~isempty(dt)
    if ~isnumeric(dt) || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
        error('dt must be empty or a positive finite scalar.');
    end
    dt = double(dt);
end

[X, value_transform_used] = local_transform_phi(phi, params.value_transform);
[T, K] = size(X);
[win_len, step_len, window_sec_used, step_sec_used] = ...
    local_resolve_window_lengths(params, dt, T);

[session_info, session_source] = local_resolve_sessions(params, T);
if params.require_session_metadata && strcmp(session_source, 'single_trace_fallback')
    error(['Thresholded density requires session metadata, but none ', ...
        'was found in params or params.observable_file.']);
end

threshold_ratio = local_resolve_threshold_ratio(params.threshold_ratio, K);
[start_idx, end_idx, center_idx, window_session_idx, window_session_id] = ...
    local_session_window_indices(T, win_len, step_len, session_info);
Nwin = numel(start_idx);

density = zeros(Nwin, K, params.output_class);
threshold_by_mode = nan(1, K);
valid_count_by_window = zeros(Nwin, K, params.output_class);

mode_str = lower(char(params.threshold_mode));
for k = 1:K
    v = double(X(:, k));
    vv = v(isfinite(v));

    threshold_by_mode(k) = local_compute_threshold(v, vv, ...
        threshold_ratio(k), mode_str);

    if isnan(threshold_by_mode(k))
        density(:, k) = cast(NaN, params.output_class);
        valid_count_by_window(:, k) = cast(0, params.output_class);
        continue;
    end

    valid = isfinite(v);
    high = (v > threshold_by_mode(k)) & valid;

    cs_high = [0; cumsum(double(high))];
    cs_valid = [0; cumsum(double(valid))];
    high_counts = cs_high(end_idx + 1) - cs_high(start_idx);
    valid_counts = cs_valid(end_idx + 1) - cs_valid(start_idx);

    switch lower(params.density_denominator)
        case 'window_samples'
            denom = repmat(double(win_len), Nwin, 1);
        case 'valid_samples'
            denom = valid_counts;
        otherwise
            error('Unknown params.density_denominator = %s.', ...
                params.density_denominator);
    end

    dens_k = high_counts ./ max(denom, eps);
    dens_k(valid_counts == 0) = NaN;

    density(:, k) = cast(dens_k, params.output_class);
    valid_count_by_window(:, k) = cast(valid_counts, params.output_class);
end

if params.smooth_density
    smooth_bins = local_resolve_smooth_bins(params, dt, step_len, step_sec_used);
    if smooth_bins > 1
        density = movmean(density, smooth_bins, 1, 'Endpoints', 'shrink');
        density = cast(density, params.output_class);
    end
else
    smooth_bins = 0;
end

[t_start, t_end, t_center] = local_window_times( ...
    start_idx, end_idx, center_idx, dt, params.time_axis);

D = struct();
D.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
D.density_time_by_mode = density;
D.threshold_by_mode = threshold_by_mode(:).';
D.valid_count_time_by_mode = valid_count_by_window;
D.window_start_idx = start_idx(:);
D.window_end_idx = end_idx(:);
D.window_center_idx = center_idx(:);
D.window_session_idx = window_session_idx(:);
D.window_session_id = window_session_id(:);
D.t_start = t_start(:);
D.t_end = t_end(:);
D.t_centers = t_center(:);
D.mode_index = local_mode_index(params, K);
D.selected_mode_idx_in_original = local_get_field(params, ...
    'selected_mode_idx_in_original', []);

D.input = struct();
D.input.n_samples = T;
D.input.n_modes = K;
D.input.dt = dt;
D.input.value_transform_requested = params.value_transform;
D.input.value_transform_used = value_transform_used;
D.input.axis_order = 'time_by_mode';

D.params = params;
D.params.window_samples = win_len;
D.params.step_samples = step_len;
D.params.window_sec_resolved = window_sec_used;
D.params.step_sec_resolved = step_sec_used;
D.params.threshold_ratio_by_mode = threshold_ratio(:).';
D.params.smooth_window_bins_resolved = smooth_bins;

D.session = session_info;
D.session.source = session_source;

D.summary = struct();
D.summary.n_windows = Nwin;
D.summary.density_size = size(density);
density_double = double(density);
nonzero_density = double(density_double > 0);
nonzero_density(~isfinite(density_double)) = NaN;
D.summary.mean_density_by_mode = local_nanmean(density_double, 1);
D.summary.max_density_by_mode = local_nanmax(density_double, 1);
D.summary.nonzero_fraction_by_mode = local_nanmean(nonzero_density, 1);
D.summary.threshold_min = local_nanmin_vector(threshold_by_mode);
D.summary.threshold_max = local_nanmax_vector(threshold_by_mode);

D.artifacts = struct();
D.artifacts.mat_file = '';
D.artifacts.figure_file = '';

[D, mat_path, fig_path] = local_prepare_artifact_paths(D, params);
D.artifacts.mat_file = mat_path;
D.artifacts.figure_file = fig_path;

fig = [];
if params.make_figure || params.save_figure
    fig = local_plot_density(D, params);
    if params.save_figure
        if exist(params.figure_dir, 'dir') ~= 7
            mkdir(params.figure_dir);
        end
        exportgraphics(fig, fig_path, 'Resolution', params.figure_resolution);
    end
end

if params.save_results
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end
    density = D.density_time_by_mode; %#ok<NASGU>
    thresholds = D.threshold_by_mode; %#ok<NASGU>
    t_centers = D.t_centers; %#ok<NASGU>
    params_saved = D.params; %#ok<NASGU>
    if params.save_v7_3
        save(mat_path, 'D', 'density', 'thresholds', 't_centers', ...
            'params_saved', '-v7.3');
    else
        save(mat_path, 'D', 'density', 'thresholds', 't_centers', ...
            'params_saved');
    end
end
end


function params = local_apply_defaults(params)
if ~isfield(params, 'window_sec'), params.window_sec = 2; end
if ~isfield(params, 'step_sec'), params.step_sec = params.window_sec; end
if ~isfield(params, 'window_samples'), params.window_samples = []; end
if ~isfield(params, 'step_samples'), params.step_samples = []; end
if ~isfield(params, 'threshold_ratio'), params.threshold_ratio = 0.7; end
if ~isfield(params, 'threshold_mode'), params.threshold_mode = 'meanplusstd'; end
if ~isfield(params, 'value_transform'), params.value_transform = 'auto'; end
if ~isfield(params, 'density_denominator')
    params.density_denominator = 'window_samples';
end
if ~isfield(params, 'output_class'), params.output_class = 'single'; end
if ~isfield(params, 'smooth_density'), params.smooth_density = false; end
if ~isfield(params, 'smooth_window_sec'), params.smooth_window_sec = 0; end
if ~isfield(params, 'smooth_window_bins'), params.smooth_window_bins = []; end
if ~isfield(params, 'time_axis'), params.time_axis = []; end
if ~isfield(params, 'mode_index'), params.mode_index = []; end
if ~isfield(params, 'selected_mode_idx_in_original')
    params.selected_mode_idx_in_original = [];
end
if ~isfield(params, 'observable_file'), params.observable_file = ''; end
if ~isfield(params, 'require_session_metadata')
    params.require_session_metadata = false;
end
if ~isfield(params, 'session_start_idx'), params.session_start_idx = []; end
if ~isfield(params, 'session_end_idx'), params.session_end_idx = []; end
if ~isfield(params, 'session_lengths'), params.session_lengths = []; end
if ~isfield(params, 'session_ids'), params.session_ids = []; end
if ~isfield(params, 'session_dx'), params.session_dx = []; end

if ~isfield(params, 'save_results'), params.save_results = false; end
if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = pwd;
end
if ~isfield(params, 'save_stem') || isempty(params.save_stem)
    params.save_stem = 'thresholded_density';
end
if ~isfield(params, 'save_tag'), params.save_tag = ''; end
if ~isfield(params, 'save_v7_3'), params.save_v7_3 = true; end

if ~isfield(params, 'make_figure'), params.make_figure = false; end
if ~isfield(params, 'save_figure'), params.save_figure = false; end
if ~isfield(params, 'figure_dir') || isempty(params.figure_dir)
    params.figure_dir = params.save_dir;
end
if ~isfield(params, 'figure_visible') || isempty(params.figure_visible)
    if params.save_figure
        params.figure_visible = 'off';
    else
        params.figure_visible = 'on';
    end
end
if ~isfield(params, 'figure_position') || isempty(params.figure_position)
    params.figure_position = [120, 120, 1260, 560];
end
if ~isfield(params, 'figure_resolution'), params.figure_resolution = 200; end
if ~isfield(params, 'max_plot_modes') || isempty(params.max_plot_modes)
    params.max_plot_modes = Inf;
end
if ~isfield(params, 'title') || isempty(params.title)
    params.title = 'Thresholded Eigenfunction Density';
end
if ~isfield(params, 'colormap') || isempty(params.colormap)
    params.colormap = 'turbo';
end
if ~isfield(params, 'background_color') || isempty(params.background_color)
    params.background_color = [0, 0, 0];
end
if ~isfield(params, 'axes_color') || isempty(params.axes_color)
    params.axes_color = params.background_color;
end
if ~isfield(params, 'text_color') || isempty(params.text_color)
    params.text_color = [1, 1, 1];
end
if ~isfield(params, 'grid_color') || isempty(params.grid_color)
    params.grid_color = [0.75, 0.75, 0.75];
end
end


function local_validate_phi(phi)
if ~isnumeric(phi) || ndims(phi) ~= 2 || isempty(phi)
    error('phi must be a nonempty numeric [T x K] matrix.');
end
end


function [X, value_transform_used] = local_transform_phi(phi, value_transform)
value_transform = lower(char(value_transform));

switch value_transform
    case 'auto'
        if ~isreal(phi)
            X = abs(phi);
            value_transform_used = 'abs';
        else
            X = phi;
            value_transform_used = 'none';
        end
    case {'none', 'asis', 'as_is'}
        if ~isreal(phi)
            error(['value_transform="none" requires real-valued phi. ', ...
                'Use "abs", "real", or "auto" for complex eigenfunctions.']);
        end
        X = phi;
        value_transform_used = 'none';
    case 'abs'
        X = abs(phi);
        value_transform_used = 'abs';
    case 'real'
        X = real(phi);
        value_transform_used = 'real';
    otherwise
        error('Unknown value_transform = %s.', value_transform);
end
end


function [win_len, step_len, window_sec_used, step_sec_used] = ...
    local_resolve_window_lengths(params, dt, T)
if ~isempty(params.window_samples)
    win_len = round(double(params.window_samples));
    if win_len <= 0
        error('params.window_samples must be positive.');
    end
else
    if isempty(dt)
        error(['dt is required when params.window_samples is empty. ', ...
            'Set params.window_samples or provide dt.']);
    end
    win_len = max(1, round(double(params.window_sec) / dt));
end

if ~isempty(params.step_samples)
    step_len = round(double(params.step_samples));
    if step_len <= 0
        error('params.step_samples must be positive.');
    end
else
    if isempty(dt)
        error(['dt is required when params.step_samples is empty. ', ...
            'Set params.step_samples or provide dt.']);
    end
    step_len = max(1, round(double(params.step_sec) / dt));
end

win_len = min(win_len, T);
step_len = max(1, step_len);

if isempty(dt)
    window_sec_used = [];
    step_sec_used = [];
else
    window_sec_used = win_len * dt;
    step_sec_used = step_len * dt;
end
end


function ratio = local_resolve_threshold_ratio(threshold_ratio, K)
if ~isnumeric(threshold_ratio) || isempty(threshold_ratio)
    error('params.threshold_ratio must be a numeric scalar or 1xK vector.');
end

if isscalar(threshold_ratio)
    ratio = repmat(double(threshold_ratio), 1, K);
else
    ratio = double(threshold_ratio(:)).';
    if numel(ratio) ~= K
        error('params.threshold_ratio must be scalar or have one value per mode.');
    end
end
end


function thr = local_compute_threshold(v, vv, ratio, mode_str)
if isempty(vv)
    thr = NaN;
    return;
end

switch mode_str
    case 'quantile'
        if ratio <= 0 || ratio >= 1
            error('For quantile mode, threshold_ratio must be in (0,1).');
        end
        thr = quantile(vv, ratio);

    case 'maxfrac'
        if ratio < 0
            error('For maxfrac mode, threshold_ratio must be >= 0.');
        end
        thr = ratio * max(vv);

    case 'meanplusstd'
        thr = mean(vv) + ratio * std(vv, 0);

    otherwise
        error(['Unknown threshold_mode = %s. Use quantile, maxfrac, ', ...
            'or meanplusstd.'], mode_str);
end
end


function y = local_nanmean(A, dim)
valid = isfinite(A);
A(~valid) = 0;
counts = sum(valid, dim);
y = sum(A, dim) ./ max(counts, 1);
y(counts == 0) = NaN;
end


function y = local_nanmax(A, dim)
valid = isfinite(A);
A(~valid) = -Inf;
y = max(A, [], dim);
y(~any(valid, dim)) = NaN;
end


function y = local_nanmin_vector(v)
v = v(isfinite(v));
if isempty(v)
    y = NaN;
else
    y = min(v);
end
end


function y = local_nanmax_vector(v)
v = v(isfinite(v));
if isempty(v)
    y = NaN;
else
    y = max(v);
end
end


function [start_idx, end_idx, center_idx] = local_window_indices(T, win_len, step_len)
last_start = T - win_len + 1;
start_idx = (1:step_len:last_start).';
if isempty(start_idx)
    start_idx = 1;
end
end_idx = start_idx + win_len - 1;
center_idx = start_idx + (win_len - 1) / 2;
end


function [start_idx, end_idx, center_idx, session_idx, session_id] = ...
    local_session_window_indices(T, win_len, step_len, session_info)
start_idx = zeros(0, 1);
end_idx = zeros(0, 1);
center_idx = zeros(0, 1);
session_idx = zeros(0, 1);
session_id = zeros(0, 1);

if isempty(session_info.start_idx) || isempty(session_info.end_idx)
    [start_idx, end_idx, center_idx] = local_window_indices(T, win_len, step_len);
    session_idx = ones(numel(start_idx), 1);
    session_id = ones(numel(start_idx), 1);
    return;
end

for s = 1:numel(session_info.start_idx)
    idx1 = session_info.start_idx(s);
    idx2 = session_info.end_idx(s);
    session_len = idx2 - idx1 + 1;
    this_win_len = min(win_len, session_len);
    if this_win_len <= 0
        continue;
    end

    local_last_start = session_len - this_win_len + 1;
    local_start = (1:step_len:local_last_start).';
    if isempty(local_start)
        local_start = 1;
    end

    this_start = idx1 + local_start - 1;
    this_end = this_start + this_win_len - 1;
    this_center = this_start + (this_win_len - 1) / 2;
    n_win = numel(this_start);

    start_idx = [start_idx; this_start(:)]; %#ok<AGROW>
    end_idx = [end_idx; this_end(:)]; %#ok<AGROW>
    center_idx = [center_idx; this_center(:)]; %#ok<AGROW>
    session_idx = [session_idx; repmat(s, n_win, 1)]; %#ok<AGROW>
    session_id = [session_id; repmat(session_info.session_ids(s), n_win, 1)]; %#ok<AGROW>
end

if isempty(start_idx)
    [start_idx, end_idx, center_idx] = local_window_indices(T, win_len, step_len);
    session_idx = ones(numel(start_idx), 1);
    session_id = ones(numel(start_idx), 1);
end
end


function [session_info, source] = local_resolve_sessions(params, T)
session_info = struct();
session_info.start_idx = params.session_start_idx(:);
session_info.end_idx = params.session_end_idx(:);
session_info.session_ids = params.session_ids(:);
session_info.session_dx = params.session_dx(:);
source = 'params';

if (isempty(session_info.start_idx) || isempty(session_info.end_idx)) && ...
        isfield(params, 'session_lengths') && ~isempty(params.session_lengths)
    lengths = round(double(params.session_lengths(:)));
    lengths = lengths(isfinite(lengths) & lengths > 0);
    if ~isempty(lengths)
        session_info.end_idx = cumsum(lengths);
        session_info.start_idx = [1; session_info.end_idx(1:end-1) + 1];
        source = 'params.session_lengths';
    end
end

if isempty(session_info.start_idx) || isempty(session_info.end_idx)
    [session_info, loaded] = local_load_sessions_from_observable(params.observable_file);
    if loaded
        source = 'observable_file';
    end
end

if isempty(session_info.start_idx) || isempty(session_info.end_idx)
    session_info.start_idx = 1;
    session_info.end_idx = T;
    session_info.session_ids = 1;
    session_info.session_dx = [];
    source = 'single_trace_fallback';
end

if isempty(session_info.session_ids)
    session_info.session_ids = (1:numel(session_info.start_idx)).';
end

if numel(session_info.start_idx) ~= numel(session_info.end_idx)
    error('session_start_idx and session_end_idx must have the same length.');
end

session_info.start_idx = round(double(session_info.start_idx(:)));
session_info.end_idx = round(double(session_info.end_idx(:)));
session_info.session_ids = double(session_info.session_ids(:));

if numel(session_info.session_ids) ~= numel(session_info.start_idx)
    error('session_ids must have one value per session.');
end

bad = session_info.start_idx < 1 | session_info.end_idx > T | ...
    session_info.end_idx < session_info.start_idx;
if any(bad)
    error('Session boundaries are invalid for a trace with T=%d samples.', T);
end
end


function [session_info, loaded] = local_load_sessions_from_observable(observable_file)
loaded = false;
session_info = struct('start_idx', [], 'end_idx', [], ...
    'session_ids', [], 'session_dx', []);

if isempty(observable_file) || exist(observable_file, 'file') ~= 2
    return;
end

try
    vars = who('-file', observable_file);
catch
    return;
end

fields = intersect({'session_start_idx', 'session_end_idx', ...
    'session_lengths', 'session_ids', 'session_dx'}, vars, 'stable');
if isempty(fields)
    return;
end

S = load(observable_file, fields{:});
if isfield(S, 'session_start_idx') && isfield(S, 'session_end_idx')
    session_info.start_idx = S.session_start_idx(:);
    session_info.end_idx = S.session_end_idx(:);
elseif isfield(S, 'session_lengths')
    lengths = double(S.session_lengths(:));
    session_info.end_idx = cumsum(lengths);
    session_info.start_idx = [1; session_info.end_idx(1:end-1) + 1];
else
    return;
end

if isfield(S, 'session_ids')
    session_info.session_ids = S.session_ids(:);
end
if isfield(S, 'session_dx')
    session_info.session_dx = S.session_dx(:);
end
loaded = true;
end


function smooth_bins = local_resolve_smooth_bins(params, dt, step_len, step_sec_used)
if ~isempty(params.smooth_window_bins)
    smooth_bins = max(1, round(double(params.smooth_window_bins)));
    return;
end

if params.smooth_window_sec <= 0
    smooth_bins = 1;
    return;
end

if ~isempty(step_sec_used)
    smooth_bins = max(1, round(double(params.smooth_window_sec) / step_sec_used));
elseif ~isempty(dt)
    smooth_bins = max(1, round(double(params.smooth_window_sec) / (step_len * dt)));
else
    error('Cannot resolve smooth_window_sec without dt; set smooth_window_bins.');
end
end


function [t_start, t_end, t_center] = local_window_times( ...
    start_idx, end_idx, center_idx, dt, time_axis)
if ~isempty(time_axis)
    time_axis = time_axis(:);
    if numel(time_axis) < max(end_idx)
        error('params.time_axis must have at least T samples.');
    end
    t_start = time_axis(start_idx);
    t_end = time_axis(end_idx);
    t_center = interp1((1:numel(time_axis)).', time_axis, center_idx, ...
        'linear', 'extrap');
elseif ~isempty(dt)
    t_start = (start_idx - 1) * dt;
    t_end = (end_idx - 1) * dt;
    t_center = (center_idx - 1) * dt;
else
    t_start = start_idx;
    t_end = end_idx;
    t_center = center_idx;
end
end


function mode_index = local_mode_index(params, K)
if isfield(params, 'mode_index') && ~isempty(params.mode_index)
    mode_index = params.mode_index(:);
    if numel(mode_index) ~= K
        error('params.mode_index must have one entry per mode.');
    end
else
    mode_index = (1:K).';
end
end


function value = local_get_field(S, field_name, default_value)
if isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function [D, mat_path, fig_path] = local_prepare_artifact_paths(D, params)
stem = char(params.save_stem);
tag = char(params.save_tag);

ratio_tag = local_ratio_tag(params.threshold_ratio);
pieces = {stem, char(params.threshold_mode), ratio_tag};
if ~isempty(tag)
    pieces{end+1} = tag;
end

base = strjoin(pieces, '__');
base = regexprep(base, '[^\w\-]+', '_');
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = sprintf('%s__%s', base, timestamp);

mat_path = '';
fig_path = '';
if params.save_results
    mat_path = fullfile(params.save_dir, [filename, '.mat']);
end
if params.save_figure
    fig_path = fullfile(params.figure_dir, [filename, '.png']);
end

D.artifacts.file_stem = filename;
end


function tag = local_ratio_tag(ratio)
if isscalar(ratio)
    tag = sprintf('ratio_%03d', round(double(ratio) * 100));
else
    tag = sprintf('ratio_vector%d', numel(ratio));
end
end


function fig = local_plot_density(D, params)
density = D.density_time_by_mode;
K = size(density, 2);
n_modes = min(K, params.max_plot_modes);
mode_idx = D.mode_index(1:n_modes);

fig = figure( ...
    'Color', params.background_color, ...
    'Position', params.figure_position, ...
    'Visible', params.figure_visible, ...
    'Name', params.title, ...
    'NumberTitle', 'off');

ax = axes(fig);
imagesc(ax, D.t_centers, 1:n_modes, double(density(:, 1:n_modes)).');
set(ax, 'YDir', 'reverse');
set(ax, ...
    'Color', params.axes_color, ...
    'XColor', params.text_color, ...
    'YColor', params.text_color, ...
    'GridColor', params.grid_color, ...
    'LineWidth', 0.8, ...
    'FontSize', 11);
grid(ax, 'on');
box(ax, 'off');

if n_modes <= 40
    yticks(ax, 1:n_modes);
    yticklabels(ax, compose('%d', mode_idx));
end

xlabel(ax, 'Time (s)', 'Color', params.text_color);
ylabel(ax, 'Mode', 'Color', params.text_color);
title(ax, local_plot_title(D, params), ...
    'Color', params.text_color, 'Interpreter', 'none');

cb = colorbar(ax);
cb.Color = params.text_color;
cb.Label.String = 'Density (fraction above threshold)';
cb.Label.Color = params.text_color;

colormap(ax, local_colormap(params.colormap));
caxis(ax, [0, 1]);
end


function title_str = local_plot_title(D, params)
ratio_tag = local_ratio_tag(params.threshold_ratio);
if isempty(D.params.window_sec_resolved)
    win_txt = sprintf('%d samples', D.params.window_samples);
    step_txt = sprintf('%d samples', D.params.step_samples);
else
    win_txt = sprintf('%.3g s', D.params.window_sec_resolved);
    step_txt = sprintf('%.3g s', D.params.step_sec_resolved);
end

title_str = sprintf('%s | win=%s, step=%s, mode=%s, %s', ...
    params.title, win_txt, step_txt, char(params.threshold_mode), ratio_tag);
end


function cmap = local_colormap(name)
name = lower(char(name));
switch name
    case 'turbo'
        if exist('turbo', 'file') == 2 || exist('turbo', 'builtin') == 5
            cmap = turbo(256);
        else
            cmap = parula(256);
        end
    case 'parula'
        cmap = parula(256);
    case 'hot'
        cmap = hot(256);
    otherwise
        cmap = parula(256);
end
end
