function paths = plot_bold_efun_density_cross_correlation_summary(out, params)
%PLOT_BOLD_EFUN_DENSITY_CROSS_CORRELATION_SUMMARY Plot BOLD-density xcorr.

if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);
if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

paths = struct('summary_png', '', 'top_curves_png', '', 'top_overlay_png', '');

fig = figure('Color', 'w', 'Position', [80, 80, 1500, 820]);
tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

T = out.peak_table;
if isempty(T)
    text(0.5, 0.5, 'No finite correlations found', ...
        'HorizontalAlignment', 'center');
else
    nexttile;
    local_plot_peak_bars(T, 'peak_abs_corr');
    title('Max |cross-correlation|');

    nexttile;
    local_plot_peak_bars(T, 'peak_lag_sec');
    title('Optimal lag (s)');
end

sgtitle(sprintf('BOLD eigenfunction vs density cross-correlation | positive lag: %s', ...
    out.positive_lag_definition), 'Interpreter', 'none');
paths.summary_png = fullfile(params.save_dir, [params.save_tag, '_summary.png']);
exportgraphics(fig, paths.summary_png, 'Resolution', params.resolution);

fig2 = figure('Color', 'w', 'Position', [100, 80, 1400, 900]);
local_plot_top_lag_curves(out, params.top_n);
paths.top_curves_png = fullfile(params.save_dir, [params.save_tag, '_top_lag_curves.png']);
exportgraphics(fig2, paths.top_curves_png, 'Resolution', params.resolution);

fig3 = figure('Color', 'w', 'Position', [120, 80, 1500, 950]);
local_plot_top_signal_overlays(out, params.top_n);
paths.top_overlay_png = fullfile(params.save_dir, [params.save_tag, '_top_signal_overlay.png']);
exportgraphics(fig3, paths.top_overlay_png, 'Resolution', params.resolution);

if params.close_figures
    close(fig);
    close(fig2);
    close(fig3);
end
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'save_dir', pwd);
params = local_set_default(params, 'save_tag', 'bold_efun_density_xcorr');
params = local_set_default(params, 'top_n', 5);
params = local_set_default(params, 'resolution', 220);
params = local_set_default(params, 'close_figures', true);
end


function local_plot_peak_bars(T, value_name)
Tf = T(isfinite(T.(value_name)), :);
if isempty(Tf)
    text(0.5, 0.5, 'No finite values', 'HorizontalAlignment', 'center');
    return;
end
Tf = sortrows(Tf, 'peak_abs_corr', 'descend');
Tf = Tf(1:min(height(Tf), 30), :);
labels = compose('%s | %s | d%d/m%d', ...
    string(Tf.density_name), string(Tf.bold_feature), ...
    Tf.density_index, Tf.bold_mode_index);
bar(Tf.(value_name));
grid on;
set(gca, 'XTick', 1:height(Tf), 'XTickLabel', labels, ...
    'TickLabelInterpreter', 'none');
xtickangle(70);
ylabel(strrep(value_name, '_', '\_'));
end


function local_plot_top_lag_curves(out, top_n)
T = out.peak_table;
if isempty(T) || isempty(out.source_results)
    text(0.5, 0.5, 'No finite correlations found', ...
        'HorizontalAlignment', 'center');
    return;
end
Tf = T(isfinite(T.peak_abs_corr), :);
Tf = sortrows(Tf, 'peak_abs_corr', 'descend');
Tf = Tf(1:min(height(Tf), top_n), :);

tiledlayout(numel(out.source_results), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i_src = 1:numel(out.source_results)
    sr = local_get_source_result(out.source_results, i_src);
    nexttile;
    hold on;
    grid on;
    mask = strcmp(string(Tf.density_name), string(sr.best_candidate.name)) & ...
        strcmp(string(Tf.bold_feature), string(sr.best_feature.name));
    Tsrc = Tf(mask, :);
    colors = lines(max(1, height(Tsrc)));
    for i = 1:height(Tsrc)
        d = Tsrc.density_index(i);
        m = Tsrc.bold_mode_index(i);
        curve = squeeze(sr.best_maps.corr_cube(d, m, :));
        plot(out.lag_sec, curve, 'LineWidth', 1.6, 'Color', colors(i, :), ...
            'DisplayName', sprintf('%s d%d/m%d peak %.3f lag %.1fs', ...
            char(sr.best_candidate.name), d, m, ...
            Tsrc.peak_corr(i), Tsrc.peak_lag_sec(i)));
        plot(Tsrc.peak_lag_sec(i), Tsrc.peak_corr(i), 'o', ...
            'MarkerFaceColor', colors(i, :), 'MarkerEdgeColor', colors(i, :), ...
            'HandleVisibility', 'off');
    end
    xline(0, 'k--');
    yline(0, 'k-');
    xlabel('Lag (s), positive = density leads BOLD');
    ylabel('corr');
    title(sprintf('%s | condition: %s | feature: %s', ...
        sr.best_candidate.name, sr.best_candidate.field_used, sr.best_feature.name), ...
        'Interpreter', 'none');
    if height(Tsrc) > 0
        legend('Location', 'eastoutside', 'Interpreter', 'none');
    end
    hold off;
end
end


function local_plot_top_signal_overlays(out, top_n)
T = out.peak_table;
if isempty(T) || isempty(out.source_results)
    text(0.5, 0.5, 'No finite correlations found', ...
        'HorizontalAlignment', 'center');
    return;
end
Tf = T(isfinite(T.peak_abs_corr), :);
Tf = sortrows(Tf, 'peak_abs_corr', 'descend');
Tf = Tf(1:min(height(Tf), top_n), :);
if isempty(Tf)
    text(0.5, 0.5, 'No finite top correlations found', ...
        'HorizontalAlignment', 'center');
    return;
end

t = (0:(out.T - 1)).' * out.dt;
tiledlayout(height(Tf), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
for i = 1:height(Tf)
    nexttile;
    sr = local_find_source_result(out, Tf.density_name(i), Tf.bold_feature(i));
    if isempty(sr)
        continue;
    end
    d = Tf.density_index(i);
    m = Tf.bold_mode_index(i);
    density = sr.best_candidate.X(:, d);
    bold = sr.best_bold_matrix(:, m);
    lag_bins = round(Tf.peak_lag_sec(i) / out.dt);

    if lag_bins > 0
        ix = 1:(numel(t) - lag_bins);
        iy = (1 + lag_bins):numel(t);
    elseif lag_bins < 0
        a = -lag_bins;
        ix = (1 + a):numel(t);
        iy = 1:(numel(t) - a);
    else
        ix = 1:numel(t);
        iy = ix;
    end

    plot(t(ix), density(ix), 'Color', [0.1 0.35 0.75], 'LineWidth', 1.0);
    hold on;
    plot(t(ix), bold(iy), 'Color', [0.75 0.2 0.1], 'LineWidth', 1.0);
    yline(0, 'k-');
    grid on;
    xlabel('Time (s)');
    ylabel('z-score');
    title(sprintf('%s | %s | density %d / BOLD mode %d | corr %.3f | lag %.1fs', ...
        char(Tf.density_name(i)), char(Tf.bold_feature(i)), d, m, ...
        Tf.peak_corr(i), Tf.peak_lag_sec(i)), 'Interpreter', 'none');
    legend({'density(t)', 'BOLD(t+lag)'}, 'Location', 'eastoutside');
    hold off;
end
end


function sr = local_find_source_result(out, density_name, feature_name)
sr = [];
for i = 1:numel(out.source_results)
    item = local_get_source_result(out.source_results, i);
    cand = item.best_candidate;
    feat = item.best_feature;
    if strcmp(string(cand.name), string(density_name)) && ...
            strcmp(string(feat.name), string(feature_name))
        sr = item;
        return;
    end
end
end


function item = local_get_source_result(source_results, idx)
if iscell(source_results)
    item = source_results{idx};
else
    item = source_results(idx);
end
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end
