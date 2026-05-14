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

fig = figure('Color', 'w', 'Position', [80, 80, 1700, 900]);
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
local_export_png(fig, paths.summary_png, params.resolution);

fig2 = figure('Color', 'w', 'Position', [100, 80, 1800, 1000]);
local_plot_top_lag_curves(out, params.top_n);
paths.top_curves_png = fullfile(params.save_dir, [params.save_tag, '_top_lag_curves.png']);
local_export_png(fig2, paths.top_curves_png, params.resolution);

overlay_height = max(1300, 260 * max(1, params.top_n));
fig3 = figure('Color', 'w', 'Position', [120, 80, 2100, overlay_height]);
local_plot_top_signal_overlays(out, params.top_n);
paths.top_overlay_png = fullfile(params.save_dir, [params.save_tag, '_top_signal_overlay.png']);
local_export_png(fig3, paths.top_overlay_png, params.resolution);

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

    [ix, iy] = local_lag_indices(numel(t), lag_bins);
    valid_pair = local_valid_lag_pair(out, ix, iy, density, bold);
    density_plot = density(ix);
    bold_plot = bold(iy);
    density_plot(~valid_pair) = NaN;
    bold_plot(~valid_pair) = NaN;

    if Tf.peak_corr(i) < 0
        bold_plot = -bold_plot;
        bold_label = '-BOLD(t+lag)';
        flip_note = ' | BOLD sign flipped for display';
    else
        bold_label = 'BOLD(t+lag)';
        flip_note = '';
    end

    plot(t(ix), density_plot, 'Color', [0.1 0.35 0.75], 'LineWidth', 1.0);
    hold on;
    plot(t(ix), bold_plot, 'Color', [0.75 0.2 0.1], 'LineWidth', 1.0);
    local_draw_session_borders(out);
    yline(0, 'k-');
    grid on;
    xlabel('Time (s)');
    ylabel('z-score');
    title(sprintf('%s | %s | density %d / BOLD mode %d | corr %.3f | lag %.1fs%s', ...
        char(Tf.density_name(i)), char(Tf.bold_feature(i)), d, m, ...
        Tf.peak_corr(i), Tf.peak_lag_sec(i), flip_note), 'Interpreter', 'none');
    legend({'density(t)', bold_label}, 'Location', 'eastoutside');
    hold off;
end
end


function [ix, iy] = local_lag_indices(n, lag_bins)
if lag_bins > 0
    ix = (1:(n - lag_bins)).';
    iy = ((1 + lag_bins):n).';
elseif lag_bins < 0
    a = -lag_bins;
    ix = ((1 + a):n).';
    iy = (1:(n - a)).';
else
    ix = (1:n).';
    iy = ix;
end
end


function valid_pair = local_valid_lag_pair(out, ix, iy, density, bold)
base_mask = true(out.T, 1);
if isfield(out, 'border_mask') && ~isempty(out.border_mask)
    base_mask = logical(out.border_mask(:));
end
if isfield(out, 'session')
    session = out.session;
else
    session = struct();
end
session_id = local_session_id_vector(out.T, session);
valid_pair = base_mask(ix) & base_mask(iy) & ...
    (session_id(ix) == session_id(iy)) & ...
    isfinite(density(ix)) & isfinite(bold(iy));
end


function local_draw_session_borders(out)
if ~isfield(out, 'session') || isempty(out.session)
    return;
end
session = out.session;
border_idx = [];
if isfield(session, 'border_idx') && ~isempty(session.border_idx)
    border_idx = double(session.border_idx(:));
elseif isfield(session, 'session_end_idx') && numel(session.session_end_idx) > 1
    border_idx = double(session.session_end_idx(1:end-1));
end
border_idx = border_idx(isfinite(border_idx) & border_idx >= 1 & border_idx <= out.T);
for i_border = 1:numel(border_idx)
    xline((border_idx(i_border) - 1) * out.dt, '--', ...
        'Color', [0.25 0.25 0.25], 'LineWidth', 0.9, ...
        'HandleVisibility', 'off');
end
end


function ids = local_session_id_vector(T, session)
ids = ones(T, 1);
if ~isstruct(session) || ~isfield(session, 'session_start_idx') || ...
        isempty(session.session_start_idx) || ~isfield(session, 'session_end_idx') || ...
        isempty(session.session_end_idx)
    return;
end
ids(:) = NaN;
starts = double(session.session_start_idx(:));
ends = double(session.session_end_idx(:));
if isfield(session, 'session_ids') && ~isempty(session.session_ids)
    session_ids = double(session.session_ids(:));
else
    session_ids = (1:numel(starts)).';
end
for i = 1:numel(starts)
    a = max(1, starts(i));
    b = min(T, ends(i));
    if a <= b
        ids(a:b) = session_ids(i);
    end
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


function local_export_png(fig, png_file, resolution)
[png_dir, ~, ~] = fileparts(png_file);
if exist(png_dir, 'dir') ~= 7
    mkdir(png_dir);
end

tmp_dir = fullfile(tempdir, 'koopman_png_shortpath');
if exist(tmp_dir, 'dir') ~= 7
    mkdir(tmp_dir);
end
tmp_png = [tempname(tmp_dir), '.png'];
cleanup_tmp = onCleanup(@() local_delete_if_exists(tmp_png));

exportgraphics(fig, tmp_png, 'Resolution', resolution);
[ok, msg] = copyfile(tmp_png, png_file, 'f');
if ~ok
    error('plot_bold_xcorr:CopyPngFailed', ...
        'Unable to copy exported PNG to target path:\n%s\n%s', png_file, msg);
end

clear cleanup_tmp
end


function local_delete_if_exists(path_in)
if exist(path_in, 'file') == 2
    try
        delete(path_in);
    catch
    end
end
end
