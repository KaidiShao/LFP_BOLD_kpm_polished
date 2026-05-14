function fig = plot_legacy_blp_event_density_pair(E_mevt, E_sevt, plot_title, output_file)
%PLOT_LEGACY_BLP_EVENT_DENSITY_PAIR Plot mevt_pl and sevt_pl density together.

if nargin < 3 || isempty(plot_title)
    plot_title = 'Legacy PL Event Density';
end
if nargin < 4
    output_file = '';
end

fig = figure('Color', 'w', 'Position', [100, 100, 1400, 800], 'Visible', 'off');
tl = tiledlayout(fig, 2, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, plot_title, 'Interpreter', 'none', 'FontWeight', 'bold');

plot_one_panel(nexttile(tl, 1), E_mevt, 'mevt\_pl');
plot_one_panel(nexttile(tl, 2), E_sevt, 'sevt\_pl');

if ~isempty(output_file)
    out_dir = fileparts(output_file);
    if ~isempty(out_dir) && exist(out_dir, 'dir') ~= 7
        mkdir(out_dir);
    end
    exportgraphics(fig, output_file, 'Resolution', 220);
end
end


function plot_one_panel(ax, E, panel_name)
axes(ax); %#ok<LAXES>
hold(ax, 'on');
grid(ax, 'on');
box(ax, 'off');

if isempty(E)
    text(ax, 0.5, 0.5, 'missing', ...
        'Units', 'normalized', ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontSize', 14, ...
        'Color', [0.4, 0.4, 0.4]);
    title(ax, panel_name, 'Interpreter', 'tex');
    xlabel(ax, 'time (s)');
    ylabel(ax, 'event density (Hz)');
    return;
end

t = double(E.t_centers(:));
Y = double(E.smoothed_density_by_type);
labels = cellstr(E.event_labels(:));
n_types = size(Y, 2);
colors = lines(max(n_types, 1));

for k = 1:n_types
    plot(ax, t, Y(:, k), ...
        'LineWidth', 1.2, ...
        'Color', colors(k, :), ...
        'DisplayName', labels{k});
end

title_text = sprintf('%s | n = %s', panel_name, mat2str(double(E.total_event_count(:).')));
title(ax, title_text, 'Interpreter', 'tex');
xlabel(ax, 'time (s)');
ylabel(ax, sprintf('density / %.3gs bin (Hz)', double(E.bin_sec)));

if n_types > 0
    legend(ax, 'Location', 'eastoutside', 'Interpreter', 'none');
end

xlim(ax, [min(t), max(t)]);
ylim_current = ylim(ax);
if ylim_current(1) > 0
    ylim(ax, [0, ylim_current(2)]);
end
end
