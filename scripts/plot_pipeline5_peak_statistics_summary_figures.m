%PLOT_PIPELINE5_PEAK_STATISTICS_SUMMARY_FIGURES
% Draw flat summary figures from pipeline5 peak-statistics CSV files.
%
% Optional caller variables:
%   processed_root  - default E:\DataPons_processed
%   dataset_stems   - default {'e10gb1','e10fV1','e10gh1','f12m01','e10gw1'}
%   peak_stage      - default pipeline5_eigenfunction_peaks_by_state_maxabs
%   summary_dir     - default <processed_root>\summary_figures\pipeline5_peak_statistics_maxabs
%   max_figures     - default [] for all

set(groot, 'defaultFigureVisible', 'off');

if ~exist('processed_root', 'var') || isempty(processed_root)
    processed_root = 'E:\DataPons_processed';
end

if ~exist('dataset_stems', 'var') || isempty(dataset_stems)
    dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01', 'e10gw1'};
end

if ~exist('peak_stage', 'var') || isempty(peak_stage)
    peak_stage = 'pipeline5_eigenfunction_peaks_by_state_maxabs';
end

if ~exist('summary_dir', 'var') || isempty(summary_dir)
    summary_dir = fullfile(processed_root, 'summary_figures', ...
        'pipeline5_peak_statistics_maxabs');
end

if ~exist('max_figures', 'var')
    max_figures = [];
end

if exist(summary_dir, 'dir') ~= 7
    mkdir(summary_dir);
end

stats_files = local_find_stats_files(processed_root, dataset_stems, peak_stage);
if ~isempty(max_figures)
    stats_files = stats_files(1:min(numel(stats_files), max_figures));
end

fprintf('Peak statistics summary dir: %s\n', summary_dir);
fprintf('Discovered %d peak-stats CSV files.\n', numel(stats_files));

index_rows = table();
n_ok = 0;
n_failed = 0;
for i = 1:numel(stats_files)
    source_csv = stats_files(i).path;
    dataset = stats_files(i).dataset;
    variant = stats_files(i).variant;
    method_tag = stats_files(i).method_tag;
    out_name = sprintf('%s__%s__%s__peak_stats_summary.png', ...
        dataset, variant, method_tag);
    out_path = fullfile(summary_dir, out_name);

    fprintf('[%03d/%03d] %s / %s / %s\n', ...
        i, numel(stats_files), dataset, variant, method_tag);

    status = "OK";
    message = "";
    try
        local_plot_one_stats(source_csv, out_path, dataset, variant, method_tag);
        n_ok = n_ok + 1;
    catch ME
        status = "FAILED";
        message = string(ME.message);
        n_failed = n_failed + 1;
        fprintf('  FAILED: %s\n', ME.message);
    end

    one = table( ...
        string(dataset), string(variant), string(method_tag), ...
        string(source_csv), string(out_path), status, message, ...
        'VariableNames', {'dataset', 'variant', 'method_tag', ...
        'source_csv', 'figure_png', 'status', 'message'});
    index_rows = [index_rows; one]; %#ok<AGROW>
end

index_csv = fullfile(summary_dir, 'index.csv');
writetable(index_rows, index_csv);
fprintf('Done. OK: %d | failed: %d\n', n_ok, n_failed);
fprintf('Wrote index: %s\n', index_csv);


function stats_files = local_find_stats_files(processed_root, dataset_stems, peak_stage)
stats_files = struct('dataset', {}, 'variant', {}, 'method_tag', {}, 'path', {});
for d = 1:numel(dataset_stems)
    dataset = char(string(dataset_stems{d}));
    root = fullfile(processed_root, dataset, peak_stage);
    if exist(root, 'dir') ~= 7
        continue;
    end

    variant_dirs = dir(root);
    variant_dirs = variant_dirs([variant_dirs.isdir]);
    variant_dirs = variant_dirs(~ismember({variant_dirs.name}, {'.', '..'}));
    for v = 1:numel(variant_dirs)
        variant = variant_dirs(v).name;
        variant_path = fullfile(root, variant);
        method_dirs = dir(variant_path);
        method_dirs = method_dirs([method_dirs.isdir]);
        method_dirs = method_dirs(~ismember({method_dirs.name}, {'.', '..'}));
        for m = 1:numel(method_dirs)
            method_tag = method_dirs(m).name;
            method_path = fullfile(variant_path, method_tag);
            files = dir(fullfile(method_path, '*_peaks_stats.csv'));
            for f = 1:numel(files)
                stats_files(end + 1).dataset = dataset; %#ok<AGROW>
                stats_files(end).variant = variant;
                stats_files(end).method_tag = method_tag;
                stats_files(end).path = fullfile(files(f).folder, files(f).name);
            end
        end
    end
end

if isempty(stats_files)
    return;
end

keys = strings(numel(stats_files), 1);
for i = 1:numel(stats_files)
    keys(i) = sprintf('%s|%s|%s|%s', stats_files(i).dataset, ...
        stats_files(i).variant, stats_files(i).method_tag, stats_files(i).path);
end
[~, order] = sort(keys);
stats_files = stats_files(order);
end


function local_plot_one_stats(source_csv, out_path, dataset, variant, method_tag)
T = readtable(source_csv, 'TextType', 'string');
required_cols = {'state_code', 'state_label', 'component_idx', ...
    'event_mean_peak', 'mean_event_minus_baseline', ...
    'cohen_d_paired_vs_baseline', 'q_vs_baseline_paired_ttest_two_sided'};
for c = 1:numel(required_cols)
    if ~ismember(required_cols{c}, T.Properties.VariableNames)
        error('Missing required column: %s', required_cols{c});
    end
end

state_codes = unique(T.state_code(:).', 'stable');
component_idx = unique(T.component_idx(:).', 'stable');
component_idx = sort(component_idx);
state_labels = strings(numel(state_codes), 1);
for s = 1:numel(state_codes)
    first_idx = find(T.state_code == state_codes(s), 1, 'first');
    state_labels(s) = T.state_label(first_idx);
end

mean_peak = local_stats_matrix(T, 'event_mean_peak', state_codes, component_idx);
delta_peak = local_stats_matrix(T, 'mean_event_minus_baseline', state_codes, component_idx);
cohen_d = local_stats_matrix(T, 'cohen_d_paired_vs_baseline', state_codes, component_idx);
q_value = local_stats_matrix(T, 'q_vs_baseline_paired_ttest_two_sided', state_codes, component_idx);
neg_log_q = -log10(max(q_value, 1e-16));
neg_log_q(~isfinite(neg_log_q)) = NaN;

fig = figure('Color', 'w', 'Position', [100, 80, 1500, 900], ...
    'Visible', 'off', 'Name', 'Pipeline5 peak statistics', 'NumberTitle', 'off');
layout = tiledlayout(fig, 2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');
title(layout, sprintf('%s / %s / %s', dataset, variant, method_tag), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

local_heatmap_tile(nexttile(layout), mean_peak, state_labels, component_idx, ...
    'Event mean peak', 'sequential', []);
local_heatmap_tile(nexttile(layout), delta_peak, state_labels, component_idx, ...
    'Mean event minus baseline', 'diverging', q_value);
local_heatmap_tile(nexttile(layout), cohen_d, state_labels, component_idx, ...
    'Cohen d paired vs baseline', 'diverging', q_value);
local_heatmap_tile(nexttile(layout), neg_log_q, state_labels, component_idx, ...
    '-log10 q paired t-test', 'q', []);

exportgraphics(fig, out_path, 'Resolution', 180);
close(fig);
end


function M = local_stats_matrix(T, column_name, state_codes, component_idx)
M = nan(numel(state_codes), numel(component_idx));
for r = 1:numel(state_codes)
    for c = 1:numel(component_idx)
        mask = T.state_code == state_codes(r) & T.component_idx == component_idx(c);
        if any(mask)
            vals = T.(column_name)(mask);
            M(r, c) = vals(1);
        end
    end
end
end


function local_heatmap_tile(ax, M, state_labels, component_idx, title_text, mode, q_value)
imagesc(ax, M, 'AlphaData', ~isnan(M));
axis(ax, 'tight');
set(ax, 'YTick', 1:numel(state_labels), 'YTickLabel', cellstr(state_labels), ...
    'XTick', 1:numel(component_idx), ...
    'XTickLabel', compose('c%d', component_idx), ...
    'TickLabelInterpreter', 'none', 'FontSize', 9);
xlabel(ax, 'Component');
ylabel(ax, 'Event');
title(ax, title_text, 'Interpreter', 'none');
colorbar(ax);

switch mode
    case 'diverging'
        colormap(ax, local_blue_white_red(256));
        finite_vals = M(isfinite(M));
        if isempty(finite_vals)
            clim(ax, [-1, 1]);
        else
            max_abs = max(abs(finite_vals));
            if max_abs == 0
                max_abs = 1;
            end
            clim(ax, [-max_abs, max_abs]);
        end
    case 'q'
        colormap(ax, flipud(hot(256)));
        finite_vals = M(isfinite(M));
        if isempty(finite_vals)
            clim(ax, [0, 1]);
        else
            clim(ax, [0, max(1, min(16, max(finite_vals)))]);
        end
    otherwise
        colormap(ax, turbo(256));
end

for r = 1:size(M, 1)
    for c = 1:size(M, 2)
        if ~isfinite(M(r, c))
            continue;
        end
        if strcmp(mode, 'q')
            label = sprintf('%.1f', M(r, c));
        else
            label = sprintf('%.2f', M(r, c));
        end
        text(ax, c, r, label, 'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', 'FontSize', 8, ...
            'Color', local_annotation_color(M, M(r, c), mode));
        if ~isempty(q_value) && isfinite(q_value(r, c)) && q_value(r, c) <= 0.05
            text(ax, c + 0.32, r - 0.30, '*', 'HorizontalAlignment', 'center', ...
                'VerticalAlignment', 'middle', 'FontSize', 12, ...
                'FontWeight', 'bold', 'Color', [0.05, 0.05, 0.05]);
        end
    end
end
end


function color = local_annotation_color(M, value, mode)
finite_vals = M(isfinite(M));
if isempty(finite_vals)
    color = [0.05, 0.05, 0.05];
    return;
end
switch mode
    case 'diverging'
        max_abs = max(abs(finite_vals));
        if max_abs == 0
            intensity = 0.5;
        else
            intensity = abs(value) / max_abs;
        end
    otherwise
        min_val = min(finite_vals);
        max_val = max(finite_vals);
        if max_val == min_val
            intensity = 0.5;
        else
            intensity = (value - min_val) / (max_val - min_val);
        end
end
if intensity > 0.58
    color = [1, 1, 1];
else
    color = [0.05, 0.05, 0.05];
end
end


function cmap = local_blue_white_red(n)
if nargin < 1 || isempty(n)
    n = 256;
end
half = floor(n / 2);
blue = [40, 92, 166] ./ 255;
white = [247, 247, 247] ./ 255;
red = [188, 65, 52] ./ 255;

lower = [linspace(blue(1), white(1), half)', ...
    linspace(blue(2), white(2), half)', ...
    linspace(blue(3), white(3), half)'];
upper = [linspace(white(1), red(1), n - half)', ...
    linspace(white(2), red(2), n - half)', ...
    linspace(white(3), red(3), n - half)'];
cmap = [lower; upper];
end
