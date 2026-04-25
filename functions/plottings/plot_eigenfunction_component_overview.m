function [fig, plot_info] = plot_eigenfunction_component_overview(result, plot_cfg)
%PLOT_EIGENFUNCTION_COMPONENT_OVERVIEW Plot eigenfunction heatmap and temporal components.

if nargin < 1 || isempty(result)
    error('result must be provided.');
end

if nargin < 2
    plot_cfg = struct();
end

plot_cfg = local_apply_defaults(plot_cfg, result);

X_time_by_mode = result.data.efun_feature_time_by_mode;
x_axis = result.input.time_axis(:);

if isempty(x_axis)
    x_axis = (1:size(X_time_by_mode, 1)).';
end

window_idx = local_resolve_window_idx(plot_cfg.window_idx, size(X_time_by_mode, 1));
x_plot = x_axis(window_idx);

n_modes_plot = min(plot_cfg.max_modes, size(X_time_by_mode, 2));
heatmap_data = X_time_by_mode(window_idx, 1:n_modes_plot).';

[comp_raw, comp_plot, comp_source_name] = local_get_component_series(result, plot_cfg);
comp_raw = comp_raw(window_idx, :);
comp_plot = comp_plot(window_idx, :);
n_comp = size(comp_plot, 2);

fig = figure( ...
    'Color', plot_cfg.background_color, ...
    'Position', plot_cfg.figure_position, ...
    'Name', plot_cfg.title, ...
    'Visible', plot_cfg.figure_visible, ...
    'NumberTitle', 'off');

tiled = tiledlayout(fig, n_comp + 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_heat = nexttile(tiled, [2, 1]);
imagesc(ax_heat, x_plot, 1:n_modes_plot, heatmap_data);
set(ax_heat, 'YDir', 'reverse');
colormap(ax_heat, parula(256));
local_style_axes(ax_heat, plot_cfg);
ylabel(ax_heat, 'Mode');
title(ax_heat, plot_cfg.title, 'Color', plot_cfg.text_color, 'Interpreter', 'none');
cb = colorbar(ax_heat);
cb.Color = plot_cfg.text_color;
cb.Label.String = sprintf('Eigenfunction (%s)', result.feature.variant);
cb.Label.Color = plot_cfg.text_color;
set(ax_heat, 'XTickLabel', []);

component_axes = gobjects(n_comp, 1);
for i = 1:n_comp
    ax = nexttile(tiled);
    component_axes(i) = ax;
    hold(ax, 'on');

    if plot_cfg.plot_raw_under_smooth && ~isempty(comp_raw)
        plot(ax, x_plot, comp_raw(:, i), ...
            'Color', plot_cfg.raw_line_color, ...
            'LineWidth', 0.9);
    end

    plot(ax, x_plot, comp_plot(:, i), ...
        'Color', plot_cfg.component_line_color, ...
        'LineWidth', 1.2);

    local_add_event_patches(ax, plot_cfg.event_windows, plot_cfg.event_colors, plot_cfg.event_alpha);
    uistack(findobj(ax, 'Type', 'Line'), 'top');

    local_style_axes(ax, plot_cfg);
    ylabel(ax, sprintf('Comp %d', i));
    if i == n_comp
        xlabel(ax, local_get_xlabel(result));
    else
        set(ax, 'XTickLabel', []);
    end

    if i == 1
        title(ax, sprintf('Temporal Components (%s)', comp_source_name), ...
            'Color', plot_cfg.text_color, 'Interpreter', 'none');
    end
end

linkaxes([ax_heat; component_axes], 'x');
xlim(ax_heat, [x_plot(1), x_plot(end)]);

plot_info = struct();
plot_info.window_idx = window_idx;
plot_info.component_source = comp_source_name;
plot_info.n_modes_plot = n_modes_plot;
plot_info.save_path = '';

if plot_cfg.save_figure
    if exist(plot_cfg.save_dir, 'dir') ~= 7
        mkdir(plot_cfg.save_dir);
    end
    save_path = local_build_save_path(plot_cfg, result);
    exportgraphics(fig, save_path, 'Resolution', 200);
    plot_info.save_path = save_path;
end
end


function plot_cfg = local_apply_defaults(plot_cfg, result)
if ~isfield(plot_cfg, 'max_modes') || isempty(plot_cfg.max_modes)
    plot_cfg.max_modes = min(50, size(result.data.efun_feature_time_by_mode, 2));
end

if ~isfield(plot_cfg, 'window_idx')
    plot_cfg.window_idx = [];
end

if ~isfield(plot_cfg, 'component_source') || isempty(plot_cfg.component_source)
    plot_cfg.component_source = 'smooth_if_available';
end

if ~isfield(plot_cfg, 'plot_raw_under_smooth') || isempty(plot_cfg.plot_raw_under_smooth)
    plot_cfg.plot_raw_under_smooth = true;
end

if ~isfield(plot_cfg, 'figure_position') || isempty(plot_cfg.figure_position)
    plot_cfg.figure_position = [120, 80, 1180, 820];
end

if ~isfield(plot_cfg, 'figure_visible') || isempty(plot_cfg.figure_visible)
    if isfield(plot_cfg, 'save_figure') && isequal(plot_cfg.save_figure, true)
        plot_cfg.figure_visible = 'off';
    else
        plot_cfg.figure_visible = 'on';
    end
end

if ~isfield(plot_cfg, 'background_color') || isempty(plot_cfg.background_color)
    plot_cfg.background_color = [0, 0, 0];
end

if ~isfield(plot_cfg, 'axes_color') || isempty(plot_cfg.axes_color)
    plot_cfg.axes_color = plot_cfg.background_color;
end

if ~isfield(plot_cfg, 'grid_color') || isempty(plot_cfg.grid_color)
    plot_cfg.grid_color = [0.75, 0.75, 0.75];
end

if ~isfield(plot_cfg, 'text_color') || isempty(plot_cfg.text_color)
    plot_cfg.text_color = [1, 1, 1];
end

if ~isfield(plot_cfg, 'raw_line_color') || isempty(plot_cfg.raw_line_color)
    plot_cfg.raw_line_color = [0.55, 0.55, 0.55];
end

if ~isfield(plot_cfg, 'component_line_color') || isempty(plot_cfg.component_line_color)
    plot_cfg.component_line_color = [0.96, 0.96, 0.96];
end

if ~isfield(plot_cfg, 'event_windows') || isempty(plot_cfg.event_windows)
    plot_cfg.event_windows = zeros(0, 3);
end

if size(plot_cfg.event_windows, 2) == 2
    plot_cfg.event_windows(:, 3) = 1;
end

if ~isfield(plot_cfg, 'event_colors') || isempty(plot_cfg.event_colors)
    plot_cfg.event_colors = [
        0.96, 0.84, 0.62;
        0.70, 0.80, 0.97;
        0.88, 0.72, 0.80;
        0.85, 0.72, 0.95];
end

if ~isfield(plot_cfg, 'event_alpha') || isempty(plot_cfg.event_alpha)
    plot_cfg.event_alpha = 0.38;
end

if ~isfield(plot_cfg, 'title') || isempty(plot_cfg.title)
    plot_cfg.title = 'Eigenfunctions And Temporal Components';
end

if ~isfield(plot_cfg, 'save_figure') || isempty(plot_cfg.save_figure)
    plot_cfg.save_figure = false;
end

if ~isfield(plot_cfg, 'save_dir') || isempty(plot_cfg.save_dir)
    if isfield(result, 'cfg') && isfield(result.cfg, 'save') && isfield(result.cfg.save, 'dir')
        plot_cfg.save_dir = result.cfg.save.dir;
    else
        plot_cfg.save_dir = pwd;
    end
end

if ~isfield(plot_cfg, 'save_tag')
    plot_cfg.save_tag = 'ov';
end
end


function window_idx = local_resolve_window_idx(window_idx, T)
if isempty(window_idx)
    window_idx = (1:T).';
    return;
end

window_idx = window_idx(:);
window_idx = window_idx(window_idx >= 1 & window_idx <= T);
window_idx = unique(window_idx, 'stable');

if isempty(window_idx)
    error('plot_cfg.window_idx does not overlap the available sample range.');
end
end


function [comp_raw, comp_plot, source_name] = local_get_component_series(result, plot_cfg)
comp_raw = result.core.temporal_components_time_by_comp;
comp_plot = comp_raw;
source_name = 'raw';

has_smooth = isfield(result, 'summary') && ...
    isfield(result.summary, 'temporal_components_smooth_time_by_comp') && ...
    ~isempty(result.summary.temporal_components_smooth_time_by_comp);

switch lower(plot_cfg.component_source)
    case 'raw'
        source_name = 'raw';

    case 'smooth'
        if ~has_smooth
            error('No smoothed temporal components are available in result.summary.');
        end
        comp_plot = result.summary.temporal_components_smooth_time_by_comp;
        source_name = 'smooth';

    case 'smooth_if_available'
        if has_smooth
            comp_plot = result.summary.temporal_components_smooth_time_by_comp;
            source_name = 'smooth';
        end

    otherwise
        error('Unknown plot_cfg.component_source = %s.', plot_cfg.component_source);
end
end


function local_style_axes(ax, plot_cfg)
set(ax, ...
    'Color', plot_cfg.axes_color, ...
    'XColor', plot_cfg.text_color, ...
    'YColor', plot_cfg.text_color, ...
    'GridColor', plot_cfg.grid_color, ...
    'MinorGridColor', plot_cfg.grid_color, ...
    'LineWidth', 0.8, ...
    'FontSize', 11);
grid(ax, 'on');
box(ax, 'off');
end


function local_add_event_patches(ax, event_windows, event_colors, event_alpha)
if isempty(event_windows)
    return;
end

y_lim = ylim(ax);
for i = 1:size(event_windows, 1)
    event_id = event_windows(i, 3);
    color_idx = mod(max(1, round(event_id)) - 1, size(event_colors, 1)) + 1;
    x1 = event_windows(i, 1);
    x2 = event_windows(i, 2);
    patch(ax, [x1, x2, x2, x1], [y_lim(1), y_lim(1), y_lim(2), y_lim(2)], ...
        event_colors(color_idx, :), ...
        'FaceAlpha', event_alpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
ylim(ax, y_lim);
end


function xlabel_text = local_get_xlabel(result)
if isfield(result.input, 'dt') && ~isempty(result.input.dt)
    xlabel_text = 'Time (s)';
else
    xlabel_text = 'Sample Index';
end
end


function save_path = local_build_save_path(plot_cfg, result)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {'efun_comp', lower(result.meta.path_kind), lower(result.meta.feature_variant)};

if isfield(plot_cfg, 'save_tag') && ~isempty(plot_cfg.save_tag)
    pieces{end + 1} = plot_cfg.save_tag;
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(plot_cfg.save_dir, sprintf('%s__%s.png', filename, timestamp));
end
