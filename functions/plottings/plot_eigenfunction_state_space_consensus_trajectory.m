function [fig, plot_info] = plot_eigenfunction_state_space_consensus_trajectory(result, consensus_input, plot_cfg)
%PLOT_EIGENFUNCTION_STATE_SPACE_CONSENSUS_TRAJECTORY Color trajectory by consensus state code.

if nargin < 1 || isempty(result)
    error('result must be provided.');
end

if nargin < 2
    consensus_input = [];
end

if nargin < 3
    plot_cfg = struct();
end

[state_code_by_time, state_catalog, source_consensus_file] = ...
    local_resolve_consensus_input(consensus_input, plot_cfg);

plot_cfg = local_apply_defaults(plot_cfg, result, state_catalog);
[~, traj_plot, source_name] = local_get_component_series(result, plot_cfg);

if size(traj_plot, 2) < 3
    error('At least 3 temporal components are required to plot a 3D state-space trajectory.');
end

T = size(traj_plot, 1);
[state_code_by_traj, state_time_idx] = local_align_state_codes(state_code_by_time, T, plot_cfg);

full_window_idx = local_resolve_window_idx(plot_cfg.window_idx, T);
[window_idx, downsample_info] = local_apply_plot_downsample(full_window_idx, plot_cfg);
local_print(plot_cfg, ['[consensus-state-space] Using %d/%d samples ', ...
    '(downsample_step=%s, max_plot_points=%s).\n'], ...
    numel(window_idx), numel(full_window_idx), ...
    local_num_to_text(plot_cfg.downsample_step), ...
    local_num_to_text(plot_cfg.max_plot_points));
traj_plot = traj_plot(window_idx, :);
state_code_plot = double(state_code_by_traj(window_idx));
state_time_idx = state_time_idx(window_idx);

x = traj_plot(:, 1);
y = traj_plot(:, 2);
z = traj_plot(:, 3);

code_ticks = local_build_code_ticks(state_code_plot, state_catalog, plot_cfg.baseline_code);
code_labels = local_build_code_tick_labels(code_ticks, state_catalog, plot_cfg);
state_cmap = local_build_state_colormap(max(code_ticks), plot_cfg);

fig = figure( ...
    'Color', plot_cfg.background_color, ...
    'Position', plot_cfg.figure_position, ...
    'Name', plot_cfg.title, ...
    'Visible', plot_cfg.figure_visible, ...
    'NumberTitle', 'off');

ax = axes(fig);
hold(ax, 'on');

local_print(plot_cfg, '[consensus-state-space] Drawing trajectory surface...\n');
trajectory_surface = local_plot_consensus_surface(ax, x, y, z, ...
    state_code_plot, plot_cfg);

if plot_cfg.show_state_number_labels
    local_print(plot_cfg, '[consensus-state-space] Adding state labels...\n');
    local_add_state_number_labels(ax, x, y, z, state_code_plot, code_ticks, ...
        state_cmap, plot_cfg);
end

local_style_axes(ax, plot_cfg);
xlabel(ax, 'Dim 1');
ylabel(ax, 'Dim 2');
zlabel(ax, 'Dim 3');
title(ax, sprintf('%s (%s)', plot_cfg.title, source_name), ...
    'Color', plot_cfg.text_color, 'Interpreter', 'none');
view(ax, plot_cfg.view_angle(1), plot_cfg.view_angle(2));

colormap(ax, state_cmap);
clim(ax, [-0.5, max(code_ticks) + 0.5]);
cb = colorbar(ax);
cb.Color = plot_cfg.text_color;
cb.Ticks = code_ticks;
cb.TickLabels = code_labels;
cb.TickLabelInterpreter = 'none';
cb.Label.String = 'Consensus state code';
cb.Label.Color = plot_cfg.text_color;

plot_info = struct();
plot_info.window_idx = window_idx;
plot_info.n_full_window_samples = numel(full_window_idx);
plot_info.n_plotted_samples = numel(window_idx);
plot_info.downsample_step_estimate = downsample_info.step_estimate;
plot_info.requested_downsample_step = downsample_info.requested_step;
plot_info.max_plot_points_applied = downsample_info.cap_applied;
plot_info.state_time_idx = state_time_idx;
plot_info.state_code_by_sample = state_code_plot(:);
plot_info.code_ticks = code_ticks(:);
plot_info.code_tick_labels = code_labels(:);
plot_info.component_source = source_name;
plot_info.source_consensus_file = source_consensus_file;
plot_info.graphics_mode = 'single_surface';
plot_info.trajectory_surface = trajectory_surface;
plot_info.save_path = '';

if plot_cfg.save_figure
    if exist(plot_cfg.save_dir, 'dir') ~= 7
        mkdir(plot_cfg.save_dir);
    end
    save_path = local_build_save_path(plot_cfg, result);
    local_print(plot_cfg, '[consensus-state-space] Saving figure: %s\n', save_path);
    exportgraphics(fig, save_path, 'Resolution', plot_cfg.figure_resolution);
    local_print(plot_cfg, '[consensus-state-space] Saved figure.\n');
    plot_info.save_path = save_path;
end
end


function [state_code_by_time, state_catalog, source_consensus_file] = ...
    local_resolve_consensus_input(consensus_input, plot_cfg)
state_catalog = struct([]);
source_consensus_file = '';

if isempty(consensus_input) && isfield(plot_cfg, 'consensus_file') && ...
        ~isempty(plot_cfg.consensus_file)
    consensus_input = plot_cfg.consensus_file;
end

if isempty(consensus_input) && isfield(plot_cfg, 'state_code_by_time') && ...
        ~isempty(plot_cfg.state_code_by_time)
    consensus_input = plot_cfg.state_code_by_time;
end

if isempty(consensus_input)
    error('Provide consensus_input, plot_cfg.consensus_file, or plot_cfg.state_code_by_time.');
end

if isnumeric(consensus_input) || islogical(consensus_input)
    state_code_by_time = uint8(consensus_input(:));
    source_consensus_file = '[state code vector input]';
    return;
end

if ischar(consensus_input) || isstring(consensus_input)
    source_consensus_file = char(consensus_input);
    S = load(source_consensus_file);
    if ~isfield(S, 'C')
        error('Consensus file %s does not contain variable C.', source_consensus_file);
    end
    consensus_input = S.C;
end

if isstruct(consensus_input)
    if ~isfield(consensus_input, 'state_code_by_time') || ...
            isempty(consensus_input.state_code_by_time)
        error('Consensus struct must contain nonempty state_code_by_time.');
    end

    state_code_by_time = uint8(consensus_input.state_code_by_time(:));
    if isfield(consensus_input, 'state_catalog')
        state_catalog = consensus_input.state_catalog(:);
    end
    if isempty(source_consensus_file) && isfield(consensus_input, 'save_file')
        source_consensus_file = consensus_input.save_file;
    elseif isempty(source_consensus_file)
        source_consensus_file = '[consensus struct input]';
    end
    return;
end

error('Unsupported consensus_input type.');
end


function plot_cfg = local_apply_defaults(plot_cfg, result, state_catalog)
if ~isfield(plot_cfg, 'window_idx')
    plot_cfg.window_idx = [];
end

if ~isfield(plot_cfg, 'component_source') || isempty(plot_cfg.component_source)
    plot_cfg.component_source = 'smooth_if_available';
end

if ~isfield(plot_cfg, 'state_start_idx') || isempty(plot_cfg.state_start_idx)
    plot_cfg.state_start_idx = 1;
end

if ~isfield(plot_cfg, 'baseline_code') || isempty(plot_cfg.baseline_code)
    plot_cfg.baseline_code = 0;
end

if ~isfield(plot_cfg, 'baseline_label') || isempty(plot_cfg.baseline_label)
    plot_cfg.baseline_label = 'baseline';
end

if ~isfield(plot_cfg, 'state_catalog') || isempty(plot_cfg.state_catalog)
    plot_cfg.state_catalog = state_catalog;
end

if ~isfield(plot_cfg, 'figure_position') || isempty(plot_cfg.figure_position)
    plot_cfg.figure_position = [140, 90, 1120, 860];
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
    plot_cfg.grid_color = [0.78, 0.78, 0.78];
end

if ~isfield(plot_cfg, 'text_color') || isempty(plot_cfg.text_color)
    plot_cfg.text_color = [1, 1, 1];
end

if ~isfield(plot_cfg, 'baseline_line_width') || isempty(plot_cfg.baseline_line_width)
    plot_cfg.baseline_line_width = 0.75;
end

if ~isfield(plot_cfg, 'event_line_width') || isempty(plot_cfg.event_line_width)
    plot_cfg.event_line_width = 2.0;
end

if ~isfield(plot_cfg, 'single_sample_marker_size') || isempty(plot_cfg.single_sample_marker_size)
    plot_cfg.single_sample_marker_size = 18;
end

if ~isfield(plot_cfg, 'patch_line_width') || isempty(plot_cfg.patch_line_width)
    plot_cfg.patch_line_width = 1.15;
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

if ~isfield(plot_cfg, 'state_colormap') || isempty(plot_cfg.state_colormap)
    plot_cfg.state_colormap = [];
end

if ~isfield(plot_cfg, 'show_state_number_labels') || isempty(plot_cfg.show_state_number_labels)
    plot_cfg.show_state_number_labels = true;
end

if ~isfield(plot_cfg, 'label_baseline') || isempty(plot_cfg.label_baseline)
    plot_cfg.label_baseline = false;
end

if ~isfield(plot_cfg, 'state_label_marker_size') || isempty(plot_cfg.state_label_marker_size)
    plot_cfg.state_label_marker_size = 13;
end

if ~isfield(plot_cfg, 'title') || isempty(plot_cfg.title)
    plot_cfg.title = 'State-Space Trajectory By Consensus State';
end

if ~isfield(plot_cfg, 'save_figure') || isempty(plot_cfg.save_figure)
    plot_cfg.save_figure = false;
end

if ~isfield(plot_cfg, 'figure_resolution') || isempty(plot_cfg.figure_resolution)
    plot_cfg.figure_resolution = 200;
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
    plot_cfg.save_tag = 'ssc';
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


function [state_code_by_traj, state_time_idx] = local_align_state_codes(state_code_by_time, T, plot_cfg)
state_start_idx = double(plot_cfg.state_start_idx);
if ~isscalar(state_start_idx) || state_start_idx < 1 || state_start_idx ~= round(state_start_idx)
    error('plot_cfg.state_start_idx must be a positive integer.');
end

state_time_idx = (state_start_idx:(state_start_idx + T - 1)).';
if state_time_idx(end) > numel(state_code_by_time)
    error(['Consensus state vector length (%d) is shorter than requested trajectory ', ...
        'state index range [%d, %d]. Adjust plot_cfg.state_start_idx or window_idx.'], ...
        numel(state_code_by_time), state_time_idx(1), state_time_idx(end));
end

state_code_by_traj = state_code_by_time(state_time_idx);
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


function h = local_plot_consensus_surface(ax, x, y, z, state_code_plot, plot_cfg)
if numel(x) < 2
    h = scatter3(ax, x, y, z, plot_cfg.single_sample_marker_size, ...
        state_code_plot(:), 'filled', ...
        'MarkerEdgeColor', plot_cfg.text_color);
    return;
end

h = surface(ax, ...
    [x(:).'; x(:).'], ...
    [y(:).'; y(:).'], ...
    [z(:).'; z(:).'], ...
    [state_code_plot(:).'; state_code_plot(:).'], ...
    'FaceColor', 'none', ...
    'EdgeColor', 'flat', ...
    'CDataMapping', 'scaled', ...
    'LineWidth', plot_cfg.patch_line_width);
end


function code_ticks = local_build_code_ticks(state_code_plot, state_catalog, baseline_code)
code_ticks = unique([double(baseline_code); double(state_code_plot(:))], 'stable');

if ~isempty(state_catalog) && isfield(state_catalog, 'code')
    catalog_codes = double([state_catalog.code]).';
    code_ticks = unique([double(baseline_code); catalog_codes(:); code_ticks(:)], 'stable');
end

code_ticks = sort(code_ticks(:)).';
end


function code_labels = local_build_code_tick_labels(code_ticks, state_catalog, plot_cfg)
code_labels = cell(numel(code_ticks), 1);

for i = 1:numel(code_ticks)
    code = code_ticks(i);
    if code == double(plot_cfg.baseline_code)
        label = sprintf('%d %s', code, plot_cfg.baseline_label);
    else
        label = sprintf('%d state', code);
        if ~isempty(state_catalog) && isfield(state_catalog, 'code') && isfield(state_catalog, 'label')
            catalog_codes = double([state_catalog.code]);
            idx = find(catalog_codes == code, 1, 'first');
            if ~isempty(idx)
                label = sprintf('%d %s', code, char(string(state_catalog(idx).label)));
            end
        end
    end
    code_labels{i} = label;
end
end


function cmap = local_build_state_colormap(max_code, plot_cfg)
if ~isempty(plot_cfg.state_colormap)
    cmap = plot_cfg.state_colormap;
else
    cmap = [ ...
        0.34, 0.36, 0.38;  % 0 baseline
        0.86, 0.74, 0.42;  % 1 theta
        0.45, 0.64, 0.90;  % 2 gamma
        0.86, 0.55, 0.75;  % 3 ripple
        0.63, 0.53, 0.88;  % 4 theta-gamma
        0.78, 0.45, 0.60]; % 5 sharp-wave-ripple
end

needed_rows = max(1, max_code + 1);
if size(cmap, 1) < needed_rows
    extra = lines(needed_rows - size(cmap, 1));
    cmap = [cmap; extra];
end

cmap = cmap(1:needed_rows, :);
end


function color = local_color_for_code(cmap, code)
idx = max(1, round(double(code)) + 1);
idx = min(idx, size(cmap, 1));
color = cmap(idx, :);
end


function local_add_state_number_labels(ax, x, y, z, state_code_plot, code_ticks, cmap, plot_cfg)
for i = 1:numel(code_ticks)
    code = code_ticks(i);
    if code == double(plot_cfg.baseline_code) && ~plot_cfg.label_baseline
        continue;
    end

    idx = find(state_code_plot == code);
    if isempty(idx)
        continue;
    end

    x0 = local_nanmedian(x(idx));
    y0 = local_nanmedian(y(idx));
    z0 = local_nanmedian(z(idx));
    if ~all(isfinite([x0, y0, z0]))
        continue;
    end

    text(ax, x0, y0, z0, sprintf('%d', code), ...
        'Color', plot_cfg.text_color, ...
        'BackgroundColor', local_color_for_code(cmap, code), ...
        'EdgeColor', plot_cfg.text_color, ...
        'FontWeight', 'bold', ...
        'FontSize', plot_cfg.state_label_marker_size, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Margin', 2);
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
axis(ax, 'vis3d');
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
pieces = {'efun_ssc', lower(result.meta.path_kind), ...
    lower(result.meta.feature_variant)};

if isfield(plot_cfg, 'save_tag') && ~isempty(plot_cfg.save_tag)
    pieces{end + 1} = plot_cfg.save_tag;
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(plot_cfg.save_dir, sprintf('%s__%s.png', filename, timestamp));
end
