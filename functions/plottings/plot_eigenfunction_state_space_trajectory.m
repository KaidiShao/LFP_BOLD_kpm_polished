function [fig, plot_info] = plot_eigenfunction_state_space_trajectory(result, plot_cfg)
%PLOT_EIGENFUNCTION_STATE_SPACE_TRAJECTORY Plot 3D state-space trajectory in four colorings.

if nargin < 1 || isempty(result)
    error('result must be provided.');
end

if nargin < 2
    plot_cfg = struct();
end

plot_cfg = local_apply_defaults(plot_cfg, result);
[~, traj_plot, source_name] = local_get_component_series(result, plot_cfg);

if size(traj_plot, 2) < 3
    error('At least 3 temporal components are required to plot a 3D state-space trajectory.');
end

full_window_idx = local_resolve_window_idx(plot_cfg.window_idx, size(traj_plot, 1));
[window_idx, downsample_info] = local_apply_plot_downsample(full_window_idx, plot_cfg);
local_print(plot_cfg, ['[state-space] Using %d/%d samples ', ...
    '(downsample_step=%s, max_plot_points=%s).\n'], ...
    numel(window_idx), numel(full_window_idx), ...
    local_num_to_text(plot_cfg.downsample_step), ...
    local_num_to_text(plot_cfg.max_plot_points));
traj_plot = traj_plot(window_idx, :);

x = traj_plot(:, 1);
y = traj_plot(:, 2);
z = traj_plot(:, 3);

if isfield(result.input, 'time_axis') && ~isempty(result.input.time_axis)
    time_color = result.input.time_axis(window_idx);
    time_label = 'Time';
else
    time_color = window_idx;
    time_label = 'Sample';
end

fig = figure( ...
    'Color', plot_cfg.background_color, ...
    'Position', plot_cfg.figure_position, ...
    'Name', plot_cfg.title, ...
    'Visible', plot_cfg.figure_visible, ...
    'NumberTitle', 'off');

tiled = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

titles = {'Color = Dim 1', 'Color = Dim 2', 'Color = Dim 3', sprintf('Color = %s', time_label)};
color_values = {x, y, z, time_color};
cmaps = {plot_cfg.value_colormap, plot_cfg.value_colormap, plot_cfg.value_colormap, plot_cfg.time_colormap};

axes_list = gobjects(4, 1);
for i = 1:4
    local_print(plot_cfg, '[state-space] Drawing panel %d/4...\n', i);
    ax = nexttile(tiled);
    axes_list(i) = ax;
    local_plot_colored_trajectory(ax, x, y, z, color_values{i}, plot_cfg);
    colormap(ax, cmaps{i});
    local_style_axes(ax, plot_cfg);
    xlabel(ax, 'Dim 1');
    ylabel(ax, 'Dim 2');
    zlabel(ax, 'Dim 3');
    title(ax, titles{i}, 'Color', plot_cfg.text_color, 'Interpreter', 'none');
    view(ax, plot_cfg.view_angle(1), plot_cfg.view_angle(2));
    cb = colorbar(ax);
    cb.Color = plot_cfg.text_color;
end

sgtitle(tiled, sprintf('%s (%s)', plot_cfg.title, source_name), ...
    'Color', plot_cfg.text_color, 'Interpreter', 'none');

linkprop(axes_list, {'XLim', 'YLim', 'ZLim'});

plot_info = struct();
plot_info.window_idx = window_idx;
plot_info.n_full_window_samples = numel(full_window_idx);
plot_info.n_plotted_samples = numel(window_idx);
plot_info.downsample_step_estimate = downsample_info.step_estimate;
plot_info.requested_downsample_step = downsample_info.requested_step;
plot_info.max_plot_points_applied = downsample_info.cap_applied;
plot_info.component_source = source_name;
plot_info.save_path = '';

if plot_cfg.save_figure
    if exist(plot_cfg.save_dir, 'dir') ~= 7
        mkdir(plot_cfg.save_dir);
    end
    save_path = local_build_save_path(plot_cfg, result);
    local_print(plot_cfg, '[state-space] Saving figure: %s\n', save_path);
    exportgraphics(fig, save_path, 'Resolution', plot_cfg.figure_resolution);
    local_print(plot_cfg, '[state-space] Saved figure.\n');
    plot_info.save_path = save_path;
end
end


function plot_cfg = local_apply_defaults(plot_cfg, result)
if ~isfield(plot_cfg, 'window_idx')
    plot_cfg.window_idx = [];
end

if ~isfield(plot_cfg, 'component_source') || isempty(plot_cfg.component_source)
    plot_cfg.component_source = 'smooth_if_available';
end

if ~isfield(plot_cfg, 'figure_position') || isempty(plot_cfg.figure_position)
    plot_cfg.figure_position = [120, 100, 1040, 860];
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
    plot_cfg.grid_color = [0.82, 0.82, 0.82];
end

if ~isfield(plot_cfg, 'text_color') || isempty(plot_cfg.text_color)
    plot_cfg.text_color = [1, 1, 1];
end

if ~isfield(plot_cfg, 'line_width') || isempty(plot_cfg.line_width)
    plot_cfg.line_width = 1.15;
end

if ~isfield(plot_cfg, 'downsample_step') || isempty(plot_cfg.downsample_step)
    plot_cfg.downsample_step = 10;
end

if ~isfield(plot_cfg, 'max_plot_points') || isempty(plot_cfg.max_plot_points)
    plot_cfg.max_plot_points = Inf;
end

if ~isfield(plot_cfg, 'view_angle') || isempty(plot_cfg.view_angle)
    plot_cfg.view_angle = [37, 24];
end

if ~isfield(plot_cfg, 'time_colormap') || isempty(plot_cfg.time_colormap)
    plot_cfg.time_colormap = parula(256);
end

if ~isfield(plot_cfg, 'value_colormap') || isempty(plot_cfg.value_colormap)
    plot_cfg.value_colormap = turbo(256);
end

if ~isfield(plot_cfg, 'title') || isempty(plot_cfg.title)
    plot_cfg.title = 'State-Space Trajectory';
end

if ~isfield(plot_cfg, 'save_figure') || isempty(plot_cfg.save_figure)
    plot_cfg.save_figure = false;
end

if ~isfield(plot_cfg, 'figure_resolution') || isempty(plot_cfg.figure_resolution)
    plot_cfg.figure_resolution = 180;
end

if ~isfield(plot_cfg, 'verbose') || isempty(plot_cfg.verbose)
    plot_cfg.verbose = true;
end

if ~isfield(plot_cfg, 'save_dir') || isempty(plot_cfg.save_dir)
    if isfield(result, 'cfg') && isfield(result.cfg, 'save') && isfield(result.cfg.save, 'dir')
        plot_cfg.save_dir = result.cfg.save.dir;
    else
        plot_cfg.save_dir = pwd;
    end
end

if ~isfield(plot_cfg, 'save_tag')
    plot_cfg.save_tag = 'ss';
end
end


function [traj_raw, traj_plot, source_name] = local_get_component_series(result, plot_cfg)
traj_raw = result.core.temporal_components_time_by_comp;
traj_plot = traj_raw;
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
        traj_plot = result.summary.temporal_components_smooth_time_by_comp;
        source_name = 'smooth';

    case 'smooth_if_available'
        if has_smooth
            traj_plot = result.summary.temporal_components_smooth_time_by_comp;
            source_name = 'smooth';
        end

    otherwise
        error('Unknown plot_cfg.component_source = %s.', plot_cfg.component_source);
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


function [window_idx, info] = local_apply_plot_downsample(window_idx, plot_cfg)
n_full = numel(window_idx);
max_points = double(plot_cfg.max_plot_points);
requested_step = double(plot_cfg.downsample_step);

info = struct();
info.requested_step = requested_step;
info.cap_applied = false;
info.step_estimate = 1;

if isscalar(requested_step) && isfinite(requested_step) && requested_step > 1
    requested_step = max(1, round(requested_step));
    sample_pos = 1:requested_step:n_full;
    if sample_pos(end) ~= n_full
        sample_pos(end + 1) = n_full;
    end
    window_idx = window_idx(sample_pos);
    info.step_estimate = n_full / max(1, numel(window_idx));
end

if isscalar(max_points) && isfinite(max_points) && max_points > 0 && ...
        numel(window_idx) > max_points
    max_points = max(2, round(max_points));
    sample_pos = unique(round(linspace(1, numel(window_idx), max_points)), 'stable');
    window_idx = window_idx(sample_pos);
    info.cap_applied = true;
    info.step_estimate = n_full / max(1, numel(window_idx));
end
end


function local_plot_colored_trajectory(ax, x, y, z, c, plot_cfg)
surface(ax, ...
    [x.'; x.'], ...
    [y.'; y.'], ...
    [z.'; z.'], ...
    [c.'; c.'], ...
    'FaceColor', 'none', ...
    'EdgeColor', 'interp', ...
    'LineWidth', plot_cfg.line_width);
end


function local_style_axes(ax, plot_cfg)
set(ax, ...
    'Color', plot_cfg.axes_color, ...
    'XColor', plot_cfg.text_color, ...
    'YColor', plot_cfg.text_color, ...
    'ZColor', plot_cfg.text_color, ...
    'GridColor', plot_cfg.grid_color, ...
    'MinorGridColor', plot_cfg.grid_color, ...
    'LineWidth', 0.8, ...
    'FontSize', 11);
grid(ax, 'on');
box(ax, 'off');
end


function local_print(plot_cfg, varargin)
if isfield(plot_cfg, 'verbose') && isequal(plot_cfg.verbose, true)
    fprintf(varargin{:});
end
end


function txt = local_num_to_text(value)
if isempty(value)
    txt = '[]';
elseif isscalar(value) && isfinite(double(value))
    txt = sprintf('%g', double(value));
else
    txt = char(string(value));
end
end


function save_path = local_build_save_path(plot_cfg, result)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {'efun_ss', lower(result.meta.path_kind), lower(result.meta.feature_variant)};

if isfield(plot_cfg, 'save_tag') && ~isempty(plot_cfg.save_tag)
    pieces{end + 1} = plot_cfg.save_tag;
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(plot_cfg.save_dir, sprintf('%s__%s.png', filename, timestamp));
end
