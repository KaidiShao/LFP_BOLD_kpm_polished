function [fig, plot_info] = plot_eigenfunction_spectrum_diagnostics(result, plot_cfg)
%PLOT_EIGENFUNCTION_SPECTRUM_DIAGNOSTICS Plot spectrum-path geometry and clustering.

if nargin < 1 || isempty(result)
    error('result must be provided.');
end

if ~isfield(result, 'aux') || ~isfield(result.aux, 'spectrum')
    error('result.aux.spectrum is required for spectrum diagnostics.');
end

if nargin < 2
    plot_cfg = struct();
end

plot_cfg = local_apply_defaults(plot_cfg, result);

spec = result.aux.spectrum;
D = spec.distance_mode_by_mode;
Z = spec.embedding_mode_by_dim;
labels = spec.cluster_labels_mode(:);
N = size(Z, 1);

if plot_cfg.use_bilinear_evalues && ~isempty(result.data.evalues_bilinear)
    lambda_plot = result.data.evalues_bilinear(:);
    xlab = 'Real parts';
    ylab = 'Imag.parts';
    plot_unit_circle = false;
else
    lambda_plot = result.data.evalues_discrete(:);
    xlab = 'Real parts';
    ylab = 'Imag.parts';
    plot_unit_circle = plot_cfg.draw_unit_circle;
end

cluster_colors = local_build_cluster_colors(max(labels), plot_cfg.cluster_colormap);

fig = figure( ...
    'Color', plot_cfg.background_color, ...
    'Position', plot_cfg.figure_position, ...
    'Name', plot_cfg.title, ...
    'NumberTitle', 'off');

tiled = tiledlayout(fig, 1, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

ax1 = nexttile(tiled);
imagesc(ax1, 1:size(D, 2), 1:size(D, 1), D);
axis(ax1, 'square');
set(ax1, 'YDir', 'normal');
colormap(ax1, plot_cfg.distance_colormap);
local_style_axes(ax1, plot_cfg);
title(ax1, 'Correlation-based distance', 'Color', plot_cfg.text_color, 'Interpreter', 'none');
xlabel(ax1, 'Eigenfunction #');
ylabel(ax1, 'Eigenfunction #');
cb1 = colorbar(ax1);
cb1.Color = plot_cfg.text_color;

ax2 = nexttile(tiled);
local_scatter_embedding(ax2, Z, (1:N).', plot_cfg.marker_size);
colormap(ax2, plot_cfg.embedding_index_colormap);
local_style_axes(ax2, plot_cfg);
title(ax2, sprintf('%dD %s embedding', size(Z, 2), upper(result.meta.method)), ...
    'Color', plot_cfg.text_color, 'Interpreter', 'none');
local_label_embedding_axes(ax2, size(Z, 2));
cb2 = colorbar(ax2);
cb2.Color = plot_cfg.text_color;

ax3 = nexttile(tiled);
local_scatter_embedding(ax3, Z, labels, plot_cfg.marker_size);
colormap(ax3, cluster_colors);
caxis(ax3, [1, max(labels)]);
local_style_axes(ax3, plot_cfg);
title(ax3, 'Clustering in reduced space', 'Color', plot_cfg.text_color, 'Interpreter', 'none');
local_label_embedding_axes(ax3, size(Z, 2));
cb3 = colorbar(ax3);
cb3.Color = plot_cfg.text_color;
cb3.Ticks = 1:max(labels);

ax4 = nexttile(tiled);
scatter(ax4, real(lambda_plot), imag(lambda_plot), ...
    plot_cfg.marker_size, labels, 'filled', ...
    'MarkerEdgeColor', [0.92, 0.92, 0.92], ...
    'LineWidth', 0.4);
hold(ax4, 'on');
if plot_unit_circle
    theta = linspace(0, 2 * pi, 400);
    plot(ax4, cos(theta), sin(theta), '--', 'Color', [0.85, 0.85, 0.85], 'LineWidth', 0.8);
end
hold(ax4, 'off');
colormap(ax4, cluster_colors);
caxis(ax4, [1, max(labels)]);
local_style_axes(ax4, plot_cfg);
title(ax4, 'Clustering in spectrum', 'Color', plot_cfg.text_color, 'Interpreter', 'none');
xlabel(ax4, xlab);
ylabel(ax4, ylab);
axis(ax4, 'equal');
cb4 = colorbar(ax4);
cb4.Color = plot_cfg.text_color;
cb4.Ticks = 1:max(labels);

plot_info = struct();
plot_info.save_path = '';

if plot_cfg.save_figure
    if exist(plot_cfg.save_dir, 'dir') ~= 7
        mkdir(plot_cfg.save_dir);
    end
    save_path = local_build_save_path(plot_cfg, result);
    exportgraphics(fig, save_path, 'Resolution', 220);
    plot_info.save_path = save_path;
end
end


function plot_cfg = local_apply_defaults(plot_cfg, result)
if ~isfield(plot_cfg, 'figure_position') || isempty(plot_cfg.figure_position)
    plot_cfg.figure_position = [80, 120, 1260, 360];
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

if ~isfield(plot_cfg, 'distance_colormap') || isempty(plot_cfg.distance_colormap)
    plot_cfg.distance_colormap = turbo(256);
end

if ~isfield(plot_cfg, 'embedding_index_colormap') || isempty(plot_cfg.embedding_index_colormap)
    plot_cfg.embedding_index_colormap = parula(256);
end

if ~isfield(plot_cfg, 'cluster_colormap')
    plot_cfg.cluster_colormap = [];
end

if ~isfield(plot_cfg, 'marker_size') || isempty(plot_cfg.marker_size)
    plot_cfg.marker_size = 28;
end

if ~isfield(plot_cfg, 'draw_unit_circle') || isempty(plot_cfg.draw_unit_circle)
    plot_cfg.draw_unit_circle = true;
end

if ~isfield(plot_cfg, 'use_bilinear_evalues') || isempty(plot_cfg.use_bilinear_evalues)
    plot_cfg.use_bilinear_evalues = false;
end

if ~isfield(plot_cfg, 'title') || isempty(plot_cfg.title)
    plot_cfg.title = 'Spectrum Path Diagnostics';
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
    plot_cfg.save_tag = 'spectrum_diag';
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


function local_scatter_embedding(ax, Z, cval, marker_size)
if size(Z, 2) >= 3
    scatter3(ax, Z(:, 1), Z(:, 2), Z(:, 3), marker_size, cval, ...
        'filled', 'MarkerEdgeColor', [0.92, 0.92, 0.92], 'LineWidth', 0.35);
    view(ax, 3);
else
    scatter(ax, Z(:, 1), Z(:, 2), marker_size, cval, ...
        'filled', 'MarkerEdgeColor', [0.92, 0.92, 0.92], 'LineWidth', 0.35);
end
end


function local_label_embedding_axes(ax, n_dim)
xlabel(ax, 'Dim 1');
ylabel(ax, 'Dim 2');
if n_dim >= 3
    zlabel(ax, 'Dim 3');
end
end


function cmap = local_build_cluster_colors(k, cluster_colormap)
if ~isempty(cluster_colormap)
    cmap = cluster_colormap;
    if size(cmap, 1) >= k
        cmap = cmap(1:k, :);
        return;
    end
end

base = lines(max(k, 1));
if size(base, 1) >= k
    cmap = base(1:k, :);
else
    cmap = interp1(linspace(0, 1, size(base, 1)), base, linspace(0, 1, k));
end
end


function save_path = local_build_save_path(plot_cfg, result)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {'eigenfunction_spectrum', lower(result.meta.feature_variant), lower(result.meta.method)};

if isfield(plot_cfg, 'save_tag') && ~isempty(plot_cfg.save_tag)
    pieces{end + 1} = plot_cfg.save_tag; %#ok<AGROW>
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(plot_cfg.save_dir, sprintf('%s__%s.png', filename, timestamp));
end
