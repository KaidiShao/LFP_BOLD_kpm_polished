function plot_e10gb1_observable_separated_subprocess_probe()
% Split the E10gb1 standardized-csplit P8/P10 probe by BOLD observable.
%
% The goal is to test whether the current slow-variable / intrinsic-trigger
% interpretation holds separately for each BOLD observable mode.

repo_root = fileparts(fileparts(mfilename('fullpath')));
result_dir = fullfile(repo_root, 'results', ...
    'e10gb1_standardized_csplit_p8_p10_probe');
input_csv = fullfile(result_dir, 'all_top_rows.csv');
fig_dir = fullfile('E:\DataPons_processed', 'summary_figures', ...
    'pipeline11_current_analysis_summary', ...
    'p8_p10_e10gb1_standardized_csplit_probe', ...
    'observable_separated_subprocess_probe');
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

T = readtable(input_csv, 'TextType', 'string', 'Delimiter', ',');
T.observable_mode = map_observable_mode(T.observable);
T.observable_class = map_observable_class(T.observable_mode);
T.two_subprocess_label = fillmissing(T.two_subprocess_label, 'constant', "unlabelled");

top_n = 50;
dimred_summary = summarize_dimred_by_observable(T, top_n);
raw_summary = summarize_raw_timescale_by_observable(T, top_n);
top5_dimred = summarize_top_dimred_rows(T, 5);

writetable(dimred_summary, fullfile(result_dir, ...
    'observable_separated_dimred_label_summary_top50.csv'));
writetable(raw_summary, fullfile(result_dir, ...
    'observable_separated_raw_timescale_summary_top50.csv'));
writetable(top5_dimred, fullfile(result_dir, ...
    'observable_separated_dimred_top5_rows.csv'));

plot_dimred_label_fraction(dimred_summary, fig_dir, top_n);
plot_dimred_score_heatmaps(dimred_summary, fig_dir, top_n);
plot_raw_timescale_heatmaps(raw_summary, fig_dir, top_n);

fprintf('Wrote observable-separated probe to:\n  %s\n', fig_dir);
end


function summary = summarize_dimred_by_observable(T, top_n)
labels = ["theta_strict", "theta_relaxed", ...
    "ripple_gamma_no_pure_theta", "other", "unlabelled"];
rows = {};
for analysis = ["P8", "P10"]
    for family = ["efun", "deconv_efun"]
        for obs = observable_order()
            mask = T.analysis == analysis & ...
                T.feature_family == family & ...
                T.observable_mode == obs & ...
                T.source_kind == "dimred" & ...
                isfinite(T.peak_abs_corr);
            Tg = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
            if height(Tg) > top_n
                Tg = Tg(1:top_n, :);
            end
            n = height(Tg);
            counts = zeros(1, numel(labels));
            for i = 1:numel(labels)
                counts(i) = nnz(Tg.two_subprocess_label == labels(i));
            end
            if n > 0
                fractions = counts ./ n;
                median_score = median(Tg.peak_abs_corr, 'omitnan');
                max_score = max(Tg.peak_abs_corr, [], 'omitnan');
                top1_label = Tg.two_subprocess_label(1);
                top1_score = Tg.peak_abs_corr(1);
            else
                fractions = nan(size(counts));
                median_score = NaN;
                max_score = NaN;
                top1_label = "";
                top1_score = NaN;
            end
            row = {char(analysis), char(family), char(obs), ...
                char(observable_class(obs)), n, median_score, max_score, ...
                char(top1_label), top1_score};
            for i = 1:numel(labels)
                row{end + 1} = counts(i); %#ok<AGROW>
                row{end + 1} = fractions(i); %#ok<AGROW>
            end
            rows(end + 1, :) = row; %#ok<AGROW>
        end
    end
end
var_names = {'analysis', 'feature_family', 'observable_mode', ...
    'observable_class', 'n_top_rows', 'median_peak_abs_corr', ...
    'max_peak_abs_corr', 'top1_label', 'top1_peak_abs_corr'};
for i = 1:numel(labels)
    var_names{end + 1} = sprintf('%s_count', labels(i)); %#ok<AGROW>
    var_names{end + 1} = sprintf('%s_fraction', labels(i)); %#ok<AGROW>
end
summary = cell2table(rows, 'VariableNames', var_names);
end


function summary = summarize_raw_timescale_by_observable(T, top_n)
rows = {};
for analysis = ["P8", "P10"]
    for family = ["efun", "deconv_efun"]
        for obs = observable_order()
            mask = T.analysis == analysis & ...
                T.feature_family == family & ...
                T.observable_mode == obs & ...
                T.source_kind == "raw" & ...
                isfinite(T.peak_abs_corr);
            Tg = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
            if height(Tg) > top_n
                Tg = Tg(1:top_n, :);
            end
            tau = Tg.timescale_sec_preferred;
            tau = tau(isfinite(tau) & tau > 0);
            if isempty(tau)
                median_tau = NaN;
                mean_tau = NaN;
                fast_frac = NaN;
                slow_frac = NaN;
            else
                median_tau = median(tau, 'omitnan');
                mean_tau = mean(tau, 'omitnan');
                fast_frac = mean(tau < 0.05, 'omitnan');
                slow_frac = mean(tau > 1.0, 'omitnan');
            end
            if height(Tg) > 0
                median_score = median(Tg.peak_abs_corr, 'omitnan');
                max_score = max(Tg.peak_abs_corr, [], 'omitnan');
                top1_raw_idx = Tg.raw_efun_index(1);
                top1_tau = Tg.timescale_sec_preferred(1);
            else
                median_score = NaN;
                max_score = NaN;
                top1_raw_idx = NaN;
                top1_tau = NaN;
            end
            rows(end + 1, :) = {char(analysis), char(family), char(obs), ...
                char(observable_class(obs)), height(Tg), median_score, max_score, ...
                median_tau, mean_tau, fast_frac, slow_frac, ...
                top1_raw_idx, top1_tau}; %#ok<AGROW>
        end
    end
end
summary = cell2table(rows, 'VariableNames', { ...
    'analysis', 'feature_family', 'observable_mode', 'observable_class', ...
    'n_top_rows', 'median_peak_abs_corr', 'max_peak_abs_corr', ...
    'median_timescale_sec', 'mean_timescale_sec', ...
    'fraction_tau_lt_50ms', 'fraction_tau_gt_1s', ...
    'top1_raw_efun_index', 'top1_timescale_sec'});
end


function top_rows = summarize_top_dimred_rows(T, top_n)
rows = {};
for analysis = ["P8", "P10"]
    for family = ["efun", "deconv_efun"]
        for obs = observable_order()
            mask = T.analysis == analysis & ...
                T.feature_family == family & ...
                T.observable_mode == obs & ...
                T.source_kind == "dimred" & ...
                isfinite(T.peak_abs_corr);
            Tg = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
            Tg = Tg(1:min(top_n, height(Tg)), :);
            for i = 1:height(Tg)
                rows(end + 1, :) = {char(analysis), char(family), char(obs), ...
                    char(observable_class(obs)), i, char(Tg.context_id(i)), ...
                    char(Tg.bold_feature(i)), double(Tg.bold_mode_index(i)), ...
                    char(Tg.density_name(i)), char(Tg.blp_method(i)), ...
                    double(Tg.blp_k(i)), double(Tg.density_index(i)), ...
                    char(Tg.two_subprocess_label(i)), ...
                    double(Tg.peak_abs_corr(i)), double(Tg.peak_corr(i)), ...
                    double(Tg.peak_lag_sec(i))}; %#ok<AGROW>
            end
        end
    end
end
top_rows = cell2table(rows, 'VariableNames', { ...
    'analysis', 'feature_family', 'observable_mode', 'observable_class', ...
    'rank', 'context_id', 'bold_feature', 'bold_mode_index', ...
    'density_name', 'blp_method', 'blp_k', 'density_index', ...
    'label', 'peak_abs_corr', 'peak_corr', 'peak_lag_sec'});
end


function plot_dimred_label_fraction(S, fig_dir, top_n)
labels = ["theta_strict", "theta_relaxed", ...
    "ripple_gamma_no_pure_theta", "other", "unlabelled"];
colors = [
    44, 123, 182
    171, 217, 233
    215, 48, 39
    150, 150, 150
    220, 220, 220] ./ 255;

fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100, 100, 1600, 1050]);
tl = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('E10gb1 standardized csplit dimred label fraction by BOLD observable | top%d', top_n), ...
    'FontWeight', 'bold');
panels = {
    "P8", "efun"
    "P8", "deconv_efun"
    "P10", "efun"
    "P10", "deconv_efun"};
for p = 1:size(panels, 1)
    ax = nexttile(tl);
    analysis = panels{p, 1};
    family = panels{p, 2};
    rows = S(string(S.analysis) == analysis & string(S.feature_family) == family, :);
    M = nan(numel(observable_order()), numel(labels));
    for i_obs = 1:numel(observable_order())
        obs = observable_order();
        row = rows(string(rows.observable_mode) == obs(i_obs), :);
        if isempty(row)
            continue;
        end
        for i_label = 1:numel(labels)
            M(i_obs, i_label) = row.(sprintf('%s_fraction', labels(i_label)));
        end
    end
    bar(ax, M, 'stacked');
    colormap(ax, colors);
    for i = 1:numel(labels)
        ax.Children(numel(labels) - i + 1).FaceColor = colors(i, :);
    end
    ylim(ax, [0, 1]);
    set(ax, 'XTick', 1:numel(observable_order()), ...
        'XTickLabel', observable_order(), 'XTickLabelRotation', 25);
    ylabel(ax, 'fraction');
    title(ax, sprintf('%s | %s', analysis, family), 'Interpreter', 'none');
    grid(ax, 'on');
    if p == 2
        legend(ax, cellstr(labels), 'Location', 'eastoutside', 'Interpreter', 'none');
    end
end
exportgraphics(fig, fullfile(fig_dir, ...
    'observable_separated_dimred_label_fraction_top50.png'), 'Resolution', 180);
close(fig);
end


function plot_dimred_score_heatmaps(S, fig_dir, top_n)
score_field = 'median_peak_abs_corr';
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100, 100, 1350, 800]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('E10gb1 dimred xcorr strength by BOLD observable | top%d median | standardized csplit', top_n), ...
    'FontWeight', 'bold');
for i_analysis = 1:2
    analysis = ["P8", "P10"];
    ax = nexttile(tl);
    obs_order = observable_order();
    M = nan(2, numel(observable_order()));
    for i_family = 1:2
        families = ["efun", "deconv_efun"];
        for i_obs = 1:numel(observable_order())
            row = S(string(S.analysis) == analysis(i_analysis) & ...
                string(S.feature_family) == families(i_family) & ...
                string(S.observable_mode) == obs_order(i_obs), :);
            if ~isempty(row)
                M(i_family, i_obs) = row.(score_field);
            end
        end
    end
    imagesc(ax, M);
    clim(ax, [0, max(M(:), [], 'omitnan')]);
    colorbar(ax);
    set(ax, 'XTick', 1:numel(observable_order()), ...
        'XTickLabel', observable_order(), 'XTickLabelRotation', 25, ...
        'YTick', 1:2, 'YTickLabel', {'efun', 'deconv_efun'});
    title(ax, analysis(i_analysis));
    for r = 1:size(M, 1)
        for c = 1:size(M, 2)
            if isfinite(M(r, c))
                text(ax, c, r, sprintf('%.3f', M(r, c)), ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
        end
    end
end
exportgraphics(fig, fullfile(fig_dir, ...
    'observable_separated_dimred_xcorr_strength_top50.png'), 'Resolution', 180);
close(fig);
end


function plot_raw_timescale_heatmaps(S, fig_dir, top_n)
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100, 100, 1350, 800]);
tl = tiledlayout(fig, 1, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('E10gb1 raw BLP efun timescale by BOLD observable | top%d median | standardized csplit', top_n), ...
    'FontWeight', 'bold');
for i_analysis = 1:2
    analysis = ["P8", "P10"];
    ax = nexttile(tl);
    obs_order = observable_order();
    M = nan(2, numel(observable_order()));
    for i_family = 1:2
        families = ["efun", "deconv_efun"];
        for i_obs = 1:numel(observable_order())
            row = S(string(S.analysis) == analysis(i_analysis) & ...
                string(S.feature_family) == families(i_family) & ...
                string(S.observable_mode) == obs_order(i_obs), :);
            if ~isempty(row)
                M(i_family, i_obs) = row.median_timescale_sec;
            end
        end
    end
    imagesc(ax, log10(M));
    colorbar(ax);
    set(ax, 'XTick', 1:numel(observable_order()), ...
        'XTickLabel', observable_order(), 'XTickLabelRotation', 25, ...
        'YTick', 1:2, 'YTickLabel', {'efun', 'deconv_efun'});
    title(ax, analysis(i_analysis));
    for r = 1:size(M, 1)
        for c = 1:size(M, 2)
            if isfinite(M(r, c))
                text(ax, c, r, format_tau(M(r, c)), ...
                    'HorizontalAlignment', 'center', 'FontWeight', 'bold');
            end
        end
    end
end
exportgraphics(fig, fullfile(fig_dir, ...
    'observable_separated_raw_timescale_top50.png'), 'Resolution', 180);
close(fig);
end


function modes = map_observable_mode(run_tags)
modes = strings(size(run_tags));
for i = 1:numel(run_tags)
    switch lower(char(run_tags(i)))
        case 'pv_gsvd100'
            modes(i) = "global_svd100";
        case 'pv_gsvd100_ds'
            modes(i) = "gsvd100_ds";
        case 'pv_hp100'
            modes(i) = "HP_svd100";
        case {'pv_roi', 'pv_roi_mean'}
            modes(i) = "roi_mean";
        otherwise
            modes(i) = string(run_tags(i));
    end
end
end


function out = map_observable_class(modes)
out = strings(size(modes));
for i = 1:numel(modes)
    out(i) = observable_class(modes(i));
end
end


function cls = observable_class(mode)
if lower(string(mode)) == "hp_svd100"
    cls = "local_HP";
else
    cls = "global_state";
end
end


function modes = observable_order()
modes = ["global_svd100", "gsvd100_ds", "HP_svd100", "roi_mean"];
end


function s = format_tau(tau)
if ~isfinite(tau)
    s = '';
elseif tau < 0.001
    s = sprintf('%.0fus', tau * 1e6);
elseif tau < 1
    s = sprintf('%.1fms', tau * 1e3);
else
    s = sprintf('%.1fs', tau);
end
end
