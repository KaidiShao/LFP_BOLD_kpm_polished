function [overview_png, overview_fig] = export_mua_residual_cross_correlation_overview(result, cfg, figure_dir, params)
%EXPORT_MUA_RESIDUAL_CROSS_CORRELATION_OVERVIEW Export one MUA overview figure.

[~, result_stem] = fileparts(result.save_paths.main_mat);
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end
overview_png = fullfile(figure_dir, [result_stem, '_overview.png']);
overview_fig = fullfile(figure_dir, [result_stem, '_overview.fig']);

if params.mua_cross_skip_existing && exist(overview_png, 'file') == 2
    return;
end

fig = local_plot_overview(result, cfg);
drawnow;
exportgraphics(fig, overview_png, 'Resolution', params.mua_cross_overview_resolution);
if params.save_fig
    savefig(fig, overview_fig);
else
    overview_fig = '';
end
if params.close_figures && isgraphics(fig)
    close(fig);
end
end


function fig = local_plot_overview(R, cfg)
T = R.pooled_corr_table;
if isempty(T)
    error('pooled_corr_table is empty in %s', R.save_paths.main_mat);
end

pairings_to_plot = {'abs_abs', 'raw_real', 'raw_imag'};
pairing_titles = containers.Map( ...
    {'abs_abs', 'raw_real', 'raw_imag'}, ...
    {'|MUA proxy| vs |residual|', 'MUA proxy vs Re(residual)', 'MUA proxy vs Im(residual)'});

channel_labels = cellstr(string(R.channel_table.channel_label));
mode_labels = compose('m%d', double(R.mode_table.mode_rank));
n_channels = height(R.channel_table);
n_modes = height(R.mode_table);

plot_mats = cell(numel(pairings_to_plot), 1);
for i = 1:numel(pairings_to_plot)
    plot_mats{i} = local_build_corr_matrix(T, pairings_to_plot{i}, R.channel_table, R.mode_table);
end

all_vals = [];
for i = 1:numel(plot_mats)
    vals = plot_mats{i};
    all_vals = [all_vals; vals(isfinite(vals))]; %#ok<AGROW>
end

if isempty(all_vals)
    clim = [-1, 1];
else
    cmax = max(abs(all_vals));
    if cmax == 0
        cmax = 1;
    end
    clim = [-cmax, cmax];
end

top_n = min(20, height(T));
top_tbl = T(1:top_n, :);
top_labels = strings(top_n, 1);
for i = 1:top_n
    top_labels(i) = sprintf('%s | %s | m%d', ...
        char(top_tbl.channel_label(i)), ...
        char(top_tbl.pairing_label(i)), ...
        double(top_tbl.mode_rank(i)));
end

fig = figure('Color', 'w', 'Position', [80, 80, 1600, 980]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(pairings_to_plot)
    nexttile(i);
    imagesc(plot_mats{i});
    axis tight;
    set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
        'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
        'TickLabelInterpreter', 'none');
    xtickangle(0);
    xlabel('Residual Mode Rank');
    ylabel('MUA-Proxy Channel');
    title(sprintf('Pooled Cross-Correlation: %s', pairing_titles(pairings_to_plot{i})));
    colormap(gca, parula(256));
    colorbar;
    caxis(clim);
end

nexttile(4);
if top_n > 0
    barh(1:top_n, top_tbl.corr, 'FaceColor', [0.17 0.48 0.72], 'EdgeColor', 'none');
    set(gca, 'YTick', 1:top_n, 'YTickLabel', top_labels, ...
        'TickLabelInterpreter', 'none', 'YDir', 'reverse');
else
    text(0.5, 0.5, 'No finite pooled cross-correlations found', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Correlation');
ylabel('Top Pairs');
title(sprintf('Top %d Pooled MUA-Residual Cross-Correlations', top_n));
grid on;
xline(0, 'k-');

sgtitle(sprintf('%s MUA-Proxy vs LFP Residual Cross-Correlation (last BLP band)\nBand index %d of %d', ...
    cfg.file_stem, double(R.source.blp_band_index), double(R.source.n_blp_bands)), ...
    'FontWeight', 'bold');
end


function M = local_build_corr_matrix(T, pairing_name, channel_table, mode_table)
n_channels = height(channel_table);
n_modes = height(mode_table);
M = nan(n_channels, n_modes);

mask = strcmp(string(T.pairing_label), string(pairing_name));
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
