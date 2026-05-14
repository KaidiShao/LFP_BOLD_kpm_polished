function [overview_png, overview_fig] = export_spkt_residual_cross_correlation_overview(result, cfg, run_info, figure_dir, params)
%EXPORT_SPKT_RESIDUAL_CROSS_CORRELATION_OVERVIEW Export one SPKT overview figure.

[~, result_stem] = fileparts(result.save_paths.main_mat);
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end
overview_png = fullfile(figure_dir, [result_stem, '_overview.png']);
overview_fig = fullfile(figure_dir, [result_stem, '_overview.fig']);

if params.spkt_cross_skip_existing && exist(overview_png, 'file') == 2
    return;
end

fig = local_plot_overview(result, cfg, run_info, params);
drawnow;
exportgraphics(fig, overview_png, 'Resolution', params.spkt_cross_overview_resolution);
if params.save_fig
    savefig(fig, overview_fig);
else
    overview_fig = '';
end
if params.close_figures && isgraphics(fig)
    close(fig);
end
end


function fig = local_plot_overview(result, cfg, run_info, params)
T = result.pooled_corr_table;
if isempty(T)
    error('pooled_corr_table is empty for %s.', run_info.run_name);
end

features_to_plot = params.spkt_cross_overview_features;
if numel(features_to_plot) < 2
    features_to_plot = [features_to_plot, features_to_plot(1)];
end
features_to_plot = features_to_plot(1:2);

channel_labels = cellstr(string(result.channel_table.channel_label));
mode_labels = compose('m%d', double(result.mode_table.mode_rank));
n_channels = height(result.channel_table);
n_modes = height(result.mode_table);

plot_mats = cell(numel(features_to_plot), 1);
for i = 1:numel(features_to_plot)
    plot_mats{i} = local_build_corr_matrix( ...
        T, features_to_plot{i}, result.channel_table, result.mode_table);
end

all_vals = [];
for i = 1:numel(plot_mats)
    vals = plot_mats{i};
    all_vals = [all_vals; vals(isfinite(vals))]; %#ok<AGROW>
end

if isempty(all_vals)
    clim_use = [-1, 1];
else
    cmax = max(abs(all_vals));
    if cmax == 0
        cmax = 1;
    end
    clim_use = [-cmax, cmax];
end

top_n = min(params.spkt_cross_overview_top_n, height(T));
top_tbl = T(1:top_n, :);
top_labels = strings(top_n, 1);
for i = 1:top_n
    top_labels(i) = sprintf('%s | %s | m%d', ...
        char(top_tbl.channel_label(i)), ...
        char(top_tbl.residual_feature(i)), ...
        double(top_tbl.mode_rank(i)));
end

fig = figure('Color', 'w', 'Position', [100, 100, 1500, 900]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(features_to_plot)
    nexttile(i);
    imagesc(plot_mats{i});
    axis tight;
    set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
        'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
        'TickLabelInterpreter', 'none');
    xtickangle(0);
    xlabel('Residual Mode Rank');
    ylabel('SPKT Channel');
    title(sprintf('Pooled Cross-Correlation Heatmap: %s', strrep(features_to_plot{i}, '_', '\_')));
    colormap(gca, parula(256));
    colorbar;
    caxis(clim_use);
end

nexttile([1, 2]);
if top_n > 0
    barh(1:top_n, top_tbl.corr, 'FaceColor', [0.18 0.46 0.71], 'EdgeColor', 'none');
    set(gca, 'YTick', 1:top_n, 'YTickLabel', top_labels, ...
        'TickLabelInterpreter', 'none', 'YDir', 'reverse');
else
    text(0.5, 0.5, 'No finite pooled cross-correlations found', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Correlation');
ylabel('Top Pairs');
title(sprintf('Top %d Pooled SPKT-Residual Cross-Correlations', top_n));
grid on;
xline(0, 'k-');

sgtitle(sprintf('%s | %s | full-time SPKT vs deconv residual cross-correlation', ...
    cfg.file_stem, run_info.run_name), 'FontWeight', 'bold', 'Interpreter', 'none');
end


function M = local_build_corr_matrix(T, feature_name, channel_table, mode_table)
n_channels = height(channel_table);
n_modes = height(mode_table);
M = nan(n_channels, n_modes);

mask = strcmp(string(T.residual_feature), string(feature_name));
Tf = T(mask, :);

for i = 1:height(Tf)
    ch_idx = find(channel_table.channel_index == Tf.channel_index(i), 1, 'first');
    mode_idx = find(mode_table.mode_rank == Tf.mode_rank(i), 1, 'first');
    if isempty(ch_idx) || isempty(mode_idx)
        continue;
    end
    M(ch_idx, mode_idx) = Tf.corr(i);
end
end
