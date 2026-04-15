function hfig = plot_blp_segment_with_events(base_plot_cache)
% Plot raw BLP traces with detected event windows overlaid in color.

if nargin < 1 || isempty(base_plot_cache)
    error('base_plot_cache must be provided. Build it first with build_blp_plot_window_cache.');
end

if ~isfield(base_plot_cache, 'show_events') || ~base_plot_cache.show_events
    error('base_plot_cache does not contain preloaded event results.');
end

t1 = double(base_plot_cache.time_range_sec(1));
t2 = double(base_plot_cache.time_range_sec(2));
t_seg = base_plot_cache.t_seg;
x_disp = base_plot_cache.x_disp;
offsets = base_plot_cache.offsets;
n_channels = base_plot_cache.n_channels;
R = base_plot_cache.R_events;
n_bands = base_plot_cache.n_event_bands;
idx_plot_1 = base_plot_cache.idx_plot_1;
idx_plot_2 = base_plot_cache.idx_plot_2;
band_colors = base_plot_cache.band_colors;
ytick_pos = base_plot_cache.ytick_pos;
ytick_labels = base_plot_cache.ytick_labels;
border_t_raw = base_plot_cache.border_t_raw;

plot_settings = base_plot_cache.plot_settings;
trace_linewidth = plot_settings.trace_linewidth;
event_linewidth = plot_settings.event_linewidth;
border_color = plot_settings.border_color;

hfig = figure('Color', 'w', 'Position', [100, 100, 1100, 460]);
ax1 = axes('Parent', hfig);
disableDefaultInteractivity(ax1);
hold(ax1, 'on');

for c = 1:n_channels
    plot(ax1, t_seg, x_disp(:, c) + offsets(c), 'k', 'LineWidth', trace_linewidth);
end

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

for i = 1:numel(border_t_raw)
    xline(ax1, border_t_raw(i), '--', 'Color', border_color, 'LineWidth', 0.8);
end

title(ax1, sprintf('BLP Signals with Detected Event Windows: %.3f-%.3f s', t1, t2));
xlim(ax1, [t1, t2]);
ylim(ax1, [ytick_pos(1) - 1, ytick_pos(end) + 1]);
set(ax1, 'YTick', ytick_pos, 'YTickLabel', ytick_labels);
set(ax1, 'Box', 'off');
xlabel(ax1, 'Time (s)');

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
