this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

result_dir = 'E:\DataPons_processed\e10gb1\spike_residual_comparison';
result_file = '';

if isempty(result_file)
    L = dir(fullfile(result_dir, '*_koopman_residual_edmd_top*.mat'));
    if isempty(L)
        error('No comparison MAT files were found in %s', result_dir);
    end
    [~, order] = sort([L.datenum], 'descend');
    L = L(order);
    result_file = fullfile(L(1).folder, L(1).name);
end

S = load(result_file, 'result');
if ~isfield(S, 'result')
    error('Result file does not contain variable "result": %s', result_file);
end

R = S.result;
T = R.pooled_corr_table;

if isempty(T)
    error('pooled_corr_table is empty in %s', result_file);
end

features_to_plot = {'abs_mean', 'abs_rms'};
channel_labels = cellstr(string(R.channel_table.channel_label));
mode_labels = compose('m%d', double(R.mode_table.mode_rank));
n_channels = height(R.channel_table);
n_modes = height(R.mode_table);

plot_mats = cell(numel(features_to_plot), 1);
for i = 1:numel(features_to_plot)
    plot_mats{i} = local_build_corr_matrix(T, features_to_plot{i}, R.channel_table, R.mode_table);
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
        char(top_tbl.residual_feature(i)), ...
        double(top_tbl.mode_rank(i)));
end

hfig = figure('Color', 'w', 'Position', [100, 100, 1500, 900]);
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
    ylabel('Spike Channel');
    title(sprintf('Pooled Correlation Heatmap: %s', strrep(features_to_plot{i}, '_', '\_')));
    colormap(gca, parula(256));
    colorbar;
    caxis(clim);
end

nexttile([1, 2]);
barh(1:top_n, top_tbl.corr, 'FaceColor', [0.18 0.46 0.71], 'EdgeColor', 'none');
set(gca, 'YTick', 1:top_n, 'YTickLabel', top_labels, 'TickLabelInterpreter', 'none', ...
    'YDir', 'reverse');
xlabel('Correlation');
ylabel('Top Pairs');
title(sprintf('Top %d Pooled Spike-Residual Correlations', top_n));
grid on;
xline(0, 'k-');

sgtitle('E10gb1 Spike vs LFP Residual Correlation', 'FontWeight', 'bold');

[~, result_stem] = fileparts(result_file);
png_file = fullfile(result_dir, [result_stem '_overview.png']);
fig_file = fullfile(result_dir, [result_stem '_overview.fig']);

exportgraphics(hfig, png_file, 'Resolution', 220);
savefig(hfig, fig_file);

fprintf('Saved correlation overview PNG:\n  %s\n', png_file);
fprintf('Saved correlation overview FIG:\n  %s\n', fig_file);


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
