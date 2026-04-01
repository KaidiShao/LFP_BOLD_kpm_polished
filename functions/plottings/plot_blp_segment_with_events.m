function hfig = plot_blp_segment_with_events(cfg, output_root, time_range_sec, event_input, band_colors)
% Plot raw BLP traces with detected event windows overlaid in color.

%% =========================
%  Input defaults
%  =========================
if nargin < 2 || isempty(output_root)
    output_root = 'D:\DataPons_processed\';
end

if nargin < 3 || isempty(time_range_sec)
    error('time_range_sec must be provided as [t_start, t_end].');
end

if nargin < 4
    event_input = [];
end

if numel(time_range_sec) ~= 2 || time_range_sec(2) <= time_range_sec(1)
    error('time_range_sec must be [t_start, t_end] with t_end > t_start.');
end

t1 = double(time_range_sec(1));
t2 = double(time_range_sec(2));

%% =========================
%  Plot settings
%  =========================
trace_scale = 0.18;
trace_clip = 4;
within_gap = 1.4;
between_gap = 2.2;
trace_linewidth = 0.4;
event_linewidth = 1.1;
border_color = [0.6 0.6 0.6];

if isfield(cfg, 'plot')
    if isfield(cfg.plot, 'trace_scale')
        trace_scale = cfg.plot.trace_scale;
    end
    if isfield(cfg.plot, 'trace_clip')
        trace_clip = cfg.plot.trace_clip;
    end
    if isfield(cfg.plot, 'within_gap')
        within_gap = cfg.plot.within_gap;
    end
    if isfield(cfg.plot, 'between_gap')
        between_gap = cfg.plot.between_gap;
    end
    if isfield(cfg.plot, 'trace_linewidth')
        trace_linewidth = cfg.plot.trace_linewidth;
    end
    if isfield(cfg.plot, 'event_linewidth')
        event_linewidth = cfg.plot.event_linewidth;
    end
end

%% =========================
%  Load raw data
%  =========================
D = load_blp_dataset(cfg);

t_raw = build_global_time_axis(D.session_lengths, D.session_dx);
t_raw = double(t_raw(:));

raw_idx = find(t_raw >= t1 & t_raw <= t2);
if isempty(raw_idx)
    error('The requested time range does not overlap the raw data.');
end

x_seg = double(D.data(raw_idx, :));
t_seg = double(t_raw(raw_idx));
idx_plot_1 = raw_idx(1);
idx_plot_2 = raw_idx(end);

%% =========================
%  Load event results
%  =========================
R = load_event_results(cfg, output_root, event_input);
n_bands = size(R.DetectResults, 1);
n_channels = size(R.DetectResults, 2);

if n_channels ~= size(D.data, 2)
    error('Event result channel count (%d) does not match loaded data (%d).', ...
        n_channels, size(D.data, 2));
end

if nargin < 5 || isempty(band_colors)
    band_colors = default_band_colors(n_bands);
end

if size(band_colors, 1) < n_bands || size(band_colors, 2) ~= 3
    error('band_colors must be an N-by-3 array with at least one row per band.');
end

%% =========================
%  Prepare raw-trace display
%  =========================
channel_sites = cfg.channels.sites(D.selected_channels);

x_disp = zeros(size(x_seg));

for c = 1:n_channels
    xc = x_seg(:, c);
    xc = xc - median(xc, 'omitnan');

    s = mad(xc, 1);
    if ~isfinite(s) || s == 0
        s = std(xc, 0, 'omitnan');
    end
    if ~isfinite(s) || s == 0
        s = 1;
    end

    xc = xc / s;
    xc = max(min(xc, trace_clip), -trace_clip);

    x_disp(:, c) = trace_scale * xc;
end

offsets = zeros(1, n_channels);
offset_val = 0;
offsets(n_channels) = offset_val;

for c = n_channels-1:-1:1
    if strcmp(channel_sites{c}, channel_sites{c+1})
        offset_val = offset_val + within_gap;
    else
        offset_val = offset_val + between_gap;
    end
    offsets(c) = offset_val;
end

offsets = double(offsets);

channel_labels = cell(1, n_channels);
for c = 1:n_channels
    channel_labels{c} = sprintf('Ch. %d', D.selected_channels(c));
end

[ytick_pos, tick_order] = sort(offsets, 'ascend');
ytick_labels = channel_labels(tick_order);

%% =========================
%  Build figure
%  =========================
hfig = figure('Color', 'w', 'Position', [100, 100, 1100, 460]);
ax1 = axes('Parent', hfig);
disableDefaultInteractivity(ax1);
hold(ax1, 'on');

%% =========================
%  Plot black traces first
%  =========================
for c = 1:n_channels
    plot(ax1, t_seg, x_disp(:, c) + offsets(c), 'k', 'LineWidth', trace_linewidth);
end

%% =========================
%  Overlay event windows in band colors
%  =========================
for c = 1:n_channels
    for b = 1:n_bands
        win = R.DetectResults{b, c}.event_win;

        if isempty(win)
            continue;
        end

        overlap_mask = win(:, 2) >= idx_plot_1 & win(:, 1) <= idx_plot_2;
        win = win(overlap_mask, :);

        for i = 1:size(win, 1)
            g1 = max(win(i, 1), idx_plot_1);
            g2 = min(win(i, 2), idx_plot_2);

            if g2 < g1
                continue;
            end

            local1 = g1 - idx_plot_1 + 1;
            local2 = g2 - idx_plot_1 + 1;

            plot(ax1, ...
                t_seg(local1:local2), ...
                x_disp(local1:local2, c) + offsets(c), ...
                'Color', band_colors(b, :), ...
                'LineWidth', event_linewidth);
        end
    end
end

%% =========================
%  Session borders
%  =========================
border_t_raw = double(t_raw(D.border_idx));
border_t_raw = border_t_raw(border_t_raw >= t1 & border_t_raw <= t2);

for i = 1:numel(border_t_raw)
    xline(ax1, border_t_raw(i), '--', 'Color', border_color, 'LineWidth', 0.8);
end

%% =========================
%  Axis styling
%  =========================
title(ax1, sprintf('BLP Signals with Detected Event Windows: %.3f-%.3f s', t1, t2));
xlim(ax1, [t1, t2]);
ylim(ax1, [ytick_pos(1) - 1, ytick_pos(end) + 1]);
set(ax1, 'YTick', ytick_pos, 'YTickLabel', ytick_labels);
set(ax1, 'Box', 'off');
xlabel(ax1, 'Time (s)');

%% =========================
%  Legend
%  =========================
legend_handles = gobjects(n_bands, 1);
legend_labels = cell(n_bands, 1);

for b = 1:n_bands
    legend_handles(b) = plot(ax1, nan, nan, 'Color', band_colors(b, :), 'LineWidth', event_linewidth + 0.6);

    if isfield(R, 'params') && isfield(R.params, 'band_labels') && numel(R.params.band_labels) >= b
        legend_labels{b} = R.params.band_labels{b};
    else
        legend_labels{b} = sprintf('Band %d', b);
    end
end

legend(ax1, legend_handles, legend_labels, 'Location', 'eastoutside');
drawnow;

end


function R = load_event_results(cfg, output_root, event_input)
% Load saved event-detection results from a struct or file.

if isstruct(event_input)
    if isfield(event_input, 'DetectResults')
        R = event_input;
        return;
    end
    error('event_input struct must contain DetectResults.');
end

if ischar(event_input) || isstring(event_input)
    S = load(char(event_input));
    if ~isfield(S, 'R')
        error('The event-result file does not contain variable R.');
    end
    R = S.R;
    return;
end

search_dir = fullfile(output_root, cfg.file_stem, 'event_detection');
pattern = fullfile(search_dir, [cfg.file_stem, '_bandpass_events_*.mat']);
L = dir(pattern);

if isempty(L)
    error('No event-result file matching %s was found.', pattern);
end

if numel(L) > 1
    is_short_name = false(numel(L), 1);
    for i = 1:numel(L)
        is_short_name(i) = ~isempty(regexp(L(i).name, ...
            ['^', regexptranslate('escape', cfg.file_stem), '_bandpass_events_[0-9]+bands\.mat$'], 'once'));
    end

    if sum(is_short_name) == 1
        L = L(is_short_name);
    else
        names = cell(numel(L), 1);
        for i = 1:numel(L)
            names{i} = fullfile(L(i).folder, L(i).name);
        end
        error('Multiple event-result files found. Please pass event_input explicitly:\n%s', strjoin(names, newline));
    end
end

S = load(fullfile(L(1).folder, L(1).name));
if ~isfield(S, 'R')
    error('The event-result file does not contain variable R.');
end

R = S.R;
end


function t = build_global_time_axis(session_lengths, session_dx)
% Build a global time axis by concatenating sessions in time.

n_sessions = numel(session_lengths);
t_cells = cell(n_sessions, 1);

for k = 1:n_sessions
    n = double(session_lengths(k));
    dx = double(session_dx(k));

    t_local = (0:n-1) * dx;

    if k == 1
        t_global = t_local;
    else
        t_prev = t_cells{k-1};
        t_global = t_local + t_prev(end) + dx;
    end

    t_cells{k} = t_global;
end

t = cat(2, t_cells{:});
t = double(t(:));
end


function colors = default_band_colors(n_bands)
% Default colors for event-band overlays.

base = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880; ...
    0.4940, 0.1840, 0.5560; ...
    0.9290, 0.6940, 0.1250];

if n_bands <= size(base, 1)
    colors = base(1:n_bands, :);
else
    colors = lines(n_bands);
end
end
