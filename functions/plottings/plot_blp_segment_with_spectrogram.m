function hfig = plot_blp_segment_with_spectrogram(base_plot_cache, cmap_in, clim)
% Plot raw BLP traces and region spectrograms for one preloaded window.

if nargin < 1 || isempty(base_plot_cache)
    error('base_plot_cache must be provided. Build it first with build_blp_plot_window_cache.');
end

if nargin < 2 || isempty(cmap_in)
    cmap_in = parula(256);
end

if nargin < 3
    clim = [];
end

t1 = double(base_plot_cache.time_range_sec(1));
t2 = double(base_plot_cache.time_range_sec(2));
t_seg = base_plot_cache.t_seg;
x_disp = base_plot_cache.x_disp;
offsets = base_plot_cache.offsets;
n_channels = base_plot_cache.n_channels;
ytick_pos = base_plot_cache.ytick_pos;
ytick_labels = base_plot_cache.ytick_labels;
border_t_raw = base_plot_cache.border_t_raw;
show_events = base_plot_cache.show_events;
R_events = base_plot_cache.R_events;
n_event_bands = base_plot_cache.n_event_bands;
idx_plot_1 = base_plot_cache.idx_plot_1;
idx_plot_2 = base_plot_cache.idx_plot_2;
band_colors = base_plot_cache.band_colors;
show_consensus = isfield(base_plot_cache, 'show_consensus') && base_plot_cache.show_consensus;
consensus_windows = [];
consensus_colors = [];
consensus_face_alpha = 0.2;
if show_consensus
    consensus_windows = base_plot_cache.consensus_windows;
    consensus_colors = base_plot_cache.consensus_colors;
    consensus_face_alpha = base_plot_cache.consensus_face_alpha;
end

if isempty(base_plot_cache.spec_seg)
    error('base_plot_cache does not contain spectrogram data.');
end

spec_seg = base_plot_cache.spec_seg;
t_spec = base_plot_cache.t_spec;
f_spec = base_plot_cache.f_spec;
regions = base_plot_cache.regions;
freq_lim = base_plot_cache.freq_lim;
spec_global_clim = [];
if isfield(base_plot_cache, 'spec_global_clim')
    spec_global_clim = base_plot_cache.spec_global_clim;
end

if isempty(clim) && ~isempty(spec_global_clim)
    clim = spec_global_clim;
end

trace_linewidth = base_plot_cache.plot_settings.trace_linewidth;

title_font_size = 22;
panel_title_font_size = 18;
axis_label_font_size = 16;
tick_font_size = 14;
legend_font_size = 18;
colorbar_font_size = 14;

n_regions = size(spec_seg, 3);
hfig = figure( ...
    'Color', 'w', ...
    'Position', [100, 100, 1000, 260 + 220 * n_regions], ...
    'MenuBar', 'none', ...
    'ToolBar', 'none');
tl = tiledlayout(n_regions + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_all = gobjects(n_regions + 1, 1);

ax1 = nexttile(tl, 1);
ax_all(1) = ax1;
disableDefaultInteractivity(ax1);
hold(ax1, 'on');

y_min = ytick_pos(1) - 1;
y_max = ytick_pos(end) + 1;

if show_consensus && ~isempty(consensus_windows)
    for i = 1:numel(consensus_windows)
        g1 = double(consensus_windows(i).win_global(1));
        g2 = double(consensus_windows(i).win_global(2));

        g1 = max(g1, idx_plot_1);
        g2 = min(g2, idx_plot_2);
        if g2 < g1
            continue;
        end

        t_start = t_seg(g1 - idx_plot_1 + 1);
        t_end = t_seg(g2 - idx_plot_1 + 1);

        state_code = double(consensus_windows(i).state_code);
        if state_code >= 1 && state_code <= size(consensus_colors, 1)
            patch_color = consensus_colors(state_code, :);
        else
            patch_color = [0.85, 0.85, 0.85];
        end

        patch(ax1, ...
            [t_start, t_end, t_end, t_start], ...
            [y_min, y_min, y_max, y_max], ...
            patch_color, ...
            'FaceAlpha', consensus_face_alpha, ...
            'EdgeColor', patch_color, ...
            'LineWidth', 1.0, ...
            'Clipping', 'on');
    end
end

for c = 1:n_channels
    plot(ax1, t_seg, x_disp(:, c) + offsets(c), 'k', 'LineWidth', trace_linewidth);
end

if show_events
    for c = 1:n_channels
        for b = 1:n_event_bands
            win = R_events.DetectResults{b, c}.event_win;
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
                    'LineWidth', max(trace_linewidth, 1.0));
            end
        end
    end
end

for i = 1:numel(border_t_raw)
    xline(ax1, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end

if show_events
    if show_consensus
        trace_title_lines = { ...
            'Original Local Field Potential Signals', ...
            sprintf('Detected Events + Consensus-State Windows | %.3f-%.3f s', t1, t2)};
    else
        trace_title_lines = { ...
            'Original Local Field Potential Signals', ...
            sprintf('Detected Events | %.3f-%.3f s', t1, t2)};
    end
else
    if show_consensus
        trace_title_lines = { ...
            'Original Local Field Potential Signals', ...
            sprintf('Consensus-State Windows | %.3f-%.3f s', t1, t2)};
    else
        trace_title_lines = { ...
            'Original Local Field Potential Signals', ...
            sprintf('%.3f-%.3f s', t1, t2)};
    end
end
title(ax1, trace_title_lines, 'FontSize', title_font_size, 'FontWeight', 'bold');
xlim(ax1, [t1, t2]);
ylim(ax1, [y_min, y_max]);
set(ax1, 'YTick', ytick_pos, 'YTickLabel', ytick_labels);
set(ax1, 'Box', 'off');
set(ax1, 'FontSize', tick_font_size);

if n_regions >= 1
    set(ax1, 'XTickLabel', []);
end

if show_events || show_consensus
    legend_handles = gobjects(0, 1);
    legend_labels = cell(0, 1);

    if show_events
        for b = 1:n_event_bands
            legend_handles(end+1, 1) = plot(ax1, nan, nan, 'Color', band_colors(b, :), 'LineWidth', 1.4); %#ok<AGROW>

            if isfield(R_events, 'params') && isfield(R_events.params, 'band_labels') && numel(R_events.params.band_labels) >= b
                legend_labels{end+1, 1} = sprintf('event: %s', R_events.params.band_labels{b}); %#ok<AGROW>
            else
                legend_labels{end+1, 1} = sprintf('event: band %d', b); %#ok<AGROW>
            end
        end
    end

    if show_consensus && ~isempty(base_plot_cache.consensus_state_catalog)
        state_catalog = base_plot_cache.consensus_state_catalog(:);
        n_states = min(numel(state_catalog), size(consensus_colors, 1));
        for s = 1:n_states
            legend_handles(end+1, 1) = patch(ax1, ... %#ok<AGROW>
                [nan, nan, nan, nan], ...
                [nan, nan, nan, nan], ...
                consensus_colors(s, :), ...
                'FaceAlpha', consensus_face_alpha, ...
                'EdgeColor', consensus_colors(s, :), ...
                'LineWidth', 1.0);
            legend_labels{end+1, 1} = sprintf('consensus: %s', char(string(state_catalog(s).label))); %#ok<AGROW>
        end
    end

    lgd = legend(ax1, legend_handles, legend_labels, 'Location', 'northeastoutside', 'Box', 'off');
    set(lgd, 'FontSize', legend_font_size);
    set(lgd, 'ItemTokenSize', [22, 18]);
end

for r = 1:n_regions
    ax = nexttile(tl, r + 1);
    ax_all(r + 1) = ax;
    disableDefaultInteractivity(ax);

    imagesc(ax, t_spec, f_spec, spec_seg(:, :, r));
    set(ax, 'YDir', 'normal');
    local_apply_colormap(ax, cmap_in);

    if ~isempty(clim)
        if size(clim, 1) == 1
            set(ax, 'CLim', clim);
        elseif size(clim, 1) >= r
            set(ax, 'CLim', clim(r, :));
        else
            error('clim must be 1x2 or N_regions x 2.');
        end
    end

    hold(ax, 'on');
    for i = 1:numel(border_t_raw)
        xline(ax, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end

    title(ax, sprintf('Power Spectrogram (%s)', regions{r}), ...
        'FontSize', panel_title_font_size, 'FontWeight', 'bold');
    ylabel(ax, 'Frequency (Hz)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
    xlim(ax, [t1, t2]);
    ylim(ax, freq_lim);
    set(ax, 'Box', 'off');
    set(ax, 'FontSize', tick_font_size);

    if r < n_regions
        set(ax, 'XTickLabel', []);
    else
        xlabel(ax, 'Time (s)', 'FontSize', axis_label_font_size, 'FontWeight', 'bold');
    end

    cb = colorbar(ax);
    cb.FontSize = colorbar_font_size;
end

linkaxes(ax_all, 'x');
drawnow;
end


function local_apply_colormap(ax, cmap_in)
if ischar(cmap_in) || isstring(cmap_in)
    cmap_name = char(cmap_in);
    colormap(ax, feval(cmap_name, 256));
elseif isa(cmap_in, 'function_handle')
    colormap(ax, cmap_in(256));
elseif isnumeric(cmap_in)
    colormap(ax, cmap_in);
else
    error('Unsupported colormap input.');
end
end
