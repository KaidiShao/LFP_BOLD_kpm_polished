this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

result_dir = 'E:\DataPons_processed\e10gb1\spike_residual_comparison';
result_file = '';

if isempty(result_file)
    L = dir(fullfile(result_dir, '*_spk_all_koopman_residual_edmd_top*.mat'));
    if isempty(L)
        L = dir(fullfile(result_dir, '*_koopman_residual_edmd_top*.mat'));
    end
    if isempty(L)
        error('No comparison MAT files were found in %s', result_dir);
    end
    [~, order] = sort([L.datenum], 'descend');
    L = L(order);
    result_file = fullfile(L(1).folder, L(1).name);
end

[~, result_stem] = fileparts(result_file);

params = struct();
params.features = {'abs_rms', 'abs_mean'};
params.max_lag_bins = 80;
params.top_n_rows = 200;
params.save_results = true;
params.save_dir = result_dir;
params.save_tag = sprintf('%s_lagged_pm%dbins', result_stem, params.max_lag_bins);
params.verbose = true;

out = compute_spike_residual_lagged_correlation(result_file, params);

feature_primary = 'abs_rms';
feature_secondary = 'abs_mean';
channel_labels = cellstr(string(out.channel_table.channel_label));
mode_labels = compose('m%d', double(out.mode_table.mode_rank));
n_channels = height(out.channel_table);
n_modes = height(out.mode_table);

peak_abs_primary = out.lag_corr.([feature_primary, '_peak']).peak_abs_corr;
peak_lag_primary = out.lag_corr.([feature_primary, '_peak']).peak_lag_sec;
peak_abs_secondary = out.lag_corr.([feature_secondary, '_peak']).peak_abs_corr;
peak_lag_secondary = out.lag_corr.([feature_secondary, '_peak']).peak_lag_sec;

mag_vals = [peak_abs_primary(:); peak_abs_secondary(:)];
mag_vals = mag_vals(isfinite(mag_vals));
if isempty(mag_vals)
    mag_clim = [0, 1];
else
    mag_clim = [0, max(mag_vals)];
    if mag_clim(2) == 0
        mag_clim(2) = 1;
    end
end

lag_vals = [peak_lag_primary(:); peak_lag_secondary(:)];
lag_vals = lag_vals(isfinite(lag_vals));
if isempty(lag_vals)
    lag_clim = [-max(abs(out.lag_sec)), max(abs(out.lag_sec))];
else
    lag_lim = max(abs(lag_vals));
    if lag_lim == 0
        lag_lim = max(abs(out.lag_sec));
    end
    lag_clim = [-lag_lim, lag_lim];
end

hfig = figure('Color', 'w', 'Position', [80, 60, 1850, 980]);
tiledlayout(2, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile(1);
imagesc(peak_abs_primary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak |corr|: %s', strrep(feature_primary, '_', '\_')));
colormap(gca, parula(256));
colorbar;
caxis(mag_clim);

nexttile(2);
imagesc(peak_lag_primary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak lag (s): %s', strrep(feature_primary, '_', '\_')));
colormap(gca, local_blue_white_red(256));
colorbar;
caxis(lag_clim);

nexttile(3, [2 1]);
local_plot_top_curves(out, feature_primary, 6);
title(sprintf('Top lag curves: %s', strrep(feature_primary, '_', '\_')));

nexttile(4);
imagesc(peak_abs_secondary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak |corr|: %s', strrep(feature_secondary, '_', '\_')));
colormap(gca, parula(256));
colorbar;
caxis(mag_clim);

nexttile(5);
imagesc(peak_lag_secondary);
axis tight;
set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
    'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(0);
xlabel('Residual Mode Rank');
ylabel('Spike Channel');
title(sprintf('Peak lag (s): %s', strrep(feature_secondary, '_', '\_')));
colormap(gca, local_blue_white_red(256));
colorbar;
caxis(lag_clim);

sgtitle(sprintf(['E10gb1 Lagged Spike vs LFP Residual Correlation\n', ...
    'Positive lag means spike leads residual | max lag = %.2f s'], max(abs(out.lag_sec))), ...
    'FontWeight', 'bold');

png_file = fullfile(result_dir, [params.save_tag, '_overview.png']);
fig_file = fullfile(result_dir, [params.save_tag, '_overview.fig']);

exportgraphics(hfig, png_file, 'Resolution', 220);
savefig(hfig, fig_file);

fprintf('Saved lagged-correlation overview PNG:\n  %s\n', png_file);
fprintf('Saved lagged-correlation overview FIG:\n  %s\n', fig_file);


function local_plot_top_curves(out, feature_name, top_n)
T = out.peak_table;
mask = strcmp(string(T.feature), string(feature_name)) & isfinite(T.peak_abs_corr);
Tf = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
top_n = min(top_n, height(Tf));

hold on;
grid on;
box on;

if top_n == 0
    text(0.5, 0.5, 'No finite lagged correlations found', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
    xlabel('Lag (s)');
    ylabel('Correlation');
    return;
end

curve_cube = out.lag_corr.(feature_name);
colors = lines(top_n);

for i = 1:top_n
    ch_row = find(out.channel_table.channel_index == Tf.channel_index(i), 1, 'first');
    mode_row = find(out.mode_table.mode_rank == Tf.mode_rank(i), 1, 'first');
    if isempty(ch_row) || isempty(mode_row)
        continue;
    end

    curve = squeeze(curve_cube(ch_row, mode_row, :));
    plot(out.lag_sec, curve, 'LineWidth', 1.6, 'Color', colors(i, :), ...
        'DisplayName', sprintf('%s | m%d | peak %.3fs', ...
        char(Tf.channel_label(i)), double(Tf.mode_rank(i)), double(Tf.peak_lag_sec(i))));
    plot(double(Tf.peak_lag_sec(i)), double(Tf.peak_corr(i)), 'o', ...
        'MarkerSize', 7, 'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), ...
        'HandleVisibility', 'off');
end

xline(0, 'k--', 'LineWidth', 1.0);
yline(0, 'k-', 'LineWidth', 0.8);
xlabel('Lag (s)');
ylabel('Correlation');
legend('Location', 'eastoutside', 'Interpreter', 'none');
end


function cmap = local_blue_white_red(n)
if nargin < 1
    n = 256;
end
half_n = ceil(n / 2);
blue = [0.13, 0.35, 0.75];
white = [1, 1, 1];
red = [0.76, 0.16, 0.12];

left = [linspace(blue(1), white(1), half_n).', ...
    linspace(blue(2), white(2), half_n).', ...
    linspace(blue(3), white(3), half_n).'];
right = [linspace(white(1), red(1), n - half_n + 1).', ...
    linspace(white(2), red(2), n - half_n + 1).', ...
    linspace(white(3), red(3), n - half_n + 1).'];
cmap = [left; right(2:end, :)];
end
