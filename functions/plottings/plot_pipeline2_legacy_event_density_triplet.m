function fig = plot_pipeline2_legacy_event_density_triplet(E_pipeline2, E_mevt, E_sevt, plot_title, output_file)
%PLOT_PIPELINE2_LEGACY_EVENT_DENSITY_TRIPLET Plot current and legacy densities.

if nargin < 4 || isempty(plot_title)
    plot_title = 'Pipeline2 And Legacy PL Event Density';
end
if nargin < 5
    output_file = '';
end

fig = figure('Color', 'w', 'Position', [100, 100, 1500, 1050], 'Visible', 'off');
tl = tiledlayout(fig, 3, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, plot_title, 'Interpreter', 'none', 'FontWeight', 'bold');

plot_pipeline2_panel(nexttile(tl, 1), E_pipeline2);
plot_legacy_panel(nexttile(tl, 2), E_mevt, 'mevt\_pl');
plot_legacy_panel(nexttile(tl, 3), E_sevt, 'sevt\_pl');

if ~isempty(output_file)
    out_dir = fileparts(output_file);
    if ~isempty(out_dir) && exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
    exportgraphics(fig, output_file, 'Resolution', 220);
end
end


function plot_pipeline2_panel(ax, E)
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'off');

if isempty(E)
    plot_missing_panel(ax, 'pipeline2');
    return;
end

t = double(E.t_centers(:));
Y = double(E.smoothed_density_mean);
if size(Y, 1) ~= numel(t) && size(Y, 2) == numel(t)
    Y = Y.';
end

labels = cellstr(E.band_labels(:));
if numel(labels) < size(Y, 2)
    for k = numel(labels)+1:size(Y, 2)
        labels{k, 1} = sprintf('band%02d', k); %#ok<AGROW>
    end
end

colors = lines(max(size(Y, 2), 1));
for k = 1:size(Y, 2)
    plot(ax, t, Y(:, k), ...
        'LineWidth', 1.2, ...
        'Color', colors(k, :), ...
        'DisplayName', labels{k});
end

title(ax, 'pipeline2 current detector', 'Interpreter', 'none');
xlabel(ax, 'time (s)');
ylabel(ax, sprintf('mean density / %.3gs bin (Hz)', double(E.bin_sec)));
if size(Y, 2) > 0
    legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
end
set_axis_limits(ax, t);
end


function plot_legacy_panel(ax, E, panel_name)
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'off');

if isempty(E)
    plot_missing_panel(ax, panel_name);
    return;
end

t = double(E.t_centers(:));
Y = double(E.smoothed_density_by_type);
if size(Y, 1) ~= numel(t) && size(Y, 2) == numel(t)
    Y = Y.';
end

labels = cellstr(E.event_labels(:));
n_types = size(Y, 2);
colors = lines(max(n_types, 1));

for k = 1:n_types
    label = labels{k};
    if isfield(E, 'total_event_count') && numel(E.total_event_count) >= k
        label = sprintf('%s (n=%d)', label, round(double(E.total_event_count(k))));
    end
    plot(ax, t, Y(:, k), ...
        'LineWidth', 1.2, ...
        'Color', colors(k, :), ...
        'DisplayName', label);
end

title(ax, panel_name, 'Interpreter', 'tex');
xlabel(ax, 'time (s)');
ylabel(ax, sprintf('density / %.3gs bin (Hz)', double(E.bin_sec)));
if n_types > 0
    legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
end
set_axis_limits(ax, t);
end


function plot_missing_panel(ax, panel_name)
text(ax, 0.5, 0.5, 'missing', ...
    'Units', 'normalized', ...
    'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', ...
    'FontSize', 14, ...
    'Color', [0.4, 0.4, 0.4]);
title(ax, panel_name, 'Interpreter', 'tex');
xlabel(ax, 'time (s)');
ylabel(ax, 'event density (Hz)');
grid(ax, 'on');
box(ax, 'off');
end


function set_axis_limits(ax, t)
if ~isempty(t)
    xlim(ax, [min(t), max(t)]);
end
ylim_current = ylim(ax);
if ylim_current(1) > 0
    ylim(ax, [0, ylim_current(2)]);
end
end
