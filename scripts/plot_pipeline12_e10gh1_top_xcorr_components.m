function plot_pipeline12_e10gh1_top_xcorr_components()
% Plot E10gH1 Pipeline12 top-xcorr dimred BLP components.
%
% Outputs:
%   results/pipeline12_standardized_csplit_consistency/e10gh1/top_xcorr_component_figures/
%
% The detailed figures show, for each P8/P10 x efun/deconv x BOLD observable
% context, the top-N dimred BLP density components and the matched BOLD-side
% target shifted by the peak lag. The gallery figures deduplicate dimred BLP
% components across contexts.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'functions')));

dataset = "e10gh1";
top_n_detail = 5;
top_n_gallery = 10;

result_dir = fullfile(repo_root, 'results', ...
    'pipeline12_standardized_csplit_consistency', char(dataset));
input_csv = fullfile(result_dir, 'p12_dimred_top_rows_with_labels.csv');
fig_dir = fullfile(result_dir, 'top_xcorr_component_figures');
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

T = readtable(input_csv, 'TextType', 'string', 'Delimiter', ',');
T = T(isfinite(T.peak_abs_corr), :);

detail_csv = fullfile(fig_dir, 'e10gh1_dimred_top_xcorr_detail_rows.csv');
if exist(detail_csv, 'file') == 2
    fprintf('Using existing detail figures listed in %s\n', detail_csv);
else
    plot_context_details(T, fig_dir, top_n_detail);
end
plot_unique_component_gallery(T, fig_dir, top_n_gallery);
end


function plot_context_details(T, fig_dir, top_n)
groups = unique(T(:, {'pipeline', 'feature_family', 'observable_mode'}), 'rows');
summary_rows = {};

for i_group = 1:height(groups)
    pipeline = groups.pipeline(i_group);
    family = groups.feature_family(i_group);
    observable = groups.observable_mode(i_group);
    mask = T.pipeline == pipeline & ...
        T.feature_family == family & ...
        T.observable_mode == observable & ...
        T.rank_within_group <= top_n;
    Tg = sortrows(T(mask, :), 'rank_within_group', 'ascend');
    if isempty(Tg)
        continue;
    end

    [png_file, rows_i] = plot_one_context(Tg, pipeline, family, observable, fig_dir);
    summary_rows = [summary_rows; rows_i]; %#ok<AGROW>
    fprintf('Wrote %s\n', png_file);
end

if ~isempty(summary_rows)
    S = cell2table(summary_rows, 'VariableNames', { ...
        'pipeline', 'feature_family', 'observable_mode', 'rank', ...
        'density_name', 'density_file', 'density_index', ...
        'blp_method', 'blp_k', 'process_label', ...
        'bold_feature', 'target_index', 'peak_abs_corr', ...
        'peak_corr', 'peak_lag_sec', 'target_file', 'figure_file'});
    writetable(S, fullfile(fig_dir, 'e10gh1_dimred_top_xcorr_detail_rows.csv'));
end
end


function [png_file, summary_rows] = plot_one_context(Tg, pipeline, family, observable, fig_dir)
n = height(Tg);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [80, 80, 1900, 270 * n + 170]);
tl = tiledlayout(fig, n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('E10gH1 dimred BLP top-xcorr components | %s | %s | %s', ...
    pipeline, family, observable), 'Interpreter', 'none', 'FontWeight', 'bold');

summary_rows = cell(n, 17);
for i = 1:n
    row = Tg(i, :);
    density_file = char(row.density_file);
    density_idx = double(row.density_index);
    target_feature = string(row.bold_feature);
    target_idx = resolve_target_index(row);
    lag_sec = double(row.peak_lag_sec);
    peak_corr = double(row.peak_corr);
    peak_abs = double(row.peak_abs_corr);
    label = string(row.two_subprocess_label);

    [density_z, density_t] = load_density_series(density_file, density_idx);
    [target_z, target_t, target_file] = load_target_series(row);
    [t_plot, x_plot, y_plot] = align_by_lag( ...
        density_z, target_z, density_t, target_t, lag_sec);

    nexttile(tl);
    hold on;
    plot(t_plot / 60, x_plot, 'Color', [0.08, 0.32, 0.80], 'LineWidth', 1.05);
    plot(t_plot / 60, y_plot, 'Color', [0.82, 0.18, 0.10], 'LineWidth', 1.05);
    yline(0, '-', 'Color', [0.72, 0.72, 0.72], 'LineWidth', 0.7);
    grid on;
    box off;
    if ~isempty(t_plot)
        xlim([min(t_plot / 60), max(t_plot / 60)]);
    end
    ylabel('z');
    if i == n
        xlabel('time (min)');
    end
    legend({'BLP dimred density', 'BOLD-side target shifted by lag'}, ...
        'Location', 'northeastoutside', 'Box', 'off');
    title(sprintf(['#%d | |r|=%.3f, r=%.3f, lag=%+.1fs | %s k%s comp%d | %s | %s idx %d\n' ...
        '%s'], ...
        double(row.rank_within_group), peak_abs, peak_corr, lag_sec, ...
        row.blp_method, string(row.blp_k), density_idx, label, ...
        target_feature, target_idx, density_file), ...
        'Interpreter', 'none', 'FontSize', 9);

    summary_rows(i, :) = { ...
        char(pipeline), char(family), char(observable), ...
        double(row.rank_within_group), char(row.density_name), ...
        density_file, density_idx, char(row.blp_method), ...
        double(row.blp_k), char(label), char(target_feature), ...
        target_idx, peak_abs, peak_corr, lag_sec, target_file, ''};
end

safe_name = sprintf('detail_top05__%s__%s__%s.png', pipeline, family, observable);
safe_name = regexprep(char(safe_name), '[^A-Za-z0-9_.-]', '_');
png_file = fullfile(fig_dir, safe_name);
exportgraphics(fig, png_file, 'Resolution', 170);
close(fig);
for i = 1:n
    summary_rows{i, 17} = png_file;
end
end


function plot_unique_component_gallery(T, fig_dir, top_n)
Tg = T(T.rank_within_group <= top_n, :);
if isempty(Tg)
    return;
end

[G, keys] = findgroups(Tg(:, {'density_file', 'density_index'}));
n_groups = max(G);
rows = cell(n_groups, 12);

for i = 1:n_groups
    Ti = Tg(G == i, :);
    Ti = sortrows(Ti, 'peak_abs_corr', 'descend');
    first = Ti(1, :);
    contexts = strings(height(Ti), 1);
    for j = 1:height(Ti)
        contexts(j) = sprintf('%s/%s/%s/r%d', ...
            Ti.pipeline(j), Ti.feature_family(j), ...
            Ti.observable_mode(j), double(Ti.rank_within_group(j)));
    end
    rows(i, :) = { ...
        char(first.density_name), char(first.density_file), ...
        double(first.density_index), char(first.blp_method), ...
        double(first.blp_k), char(first.two_subprocess_label), ...
        height(Ti), max(Ti.peak_abs_corr), mean(Ti.peak_abs_corr, 'omitnan'), ...
        strjoin(contexts, '; '), '', ''};
end

S = cell2table(rows, 'VariableNames', { ...
    'density_name', 'density_file', 'density_index', ...
    'blp_method', 'blp_k', 'process_label', ...
    'topN_occurrence_count', 'max_peak_abs_corr', ...
    'mean_peak_abs_corr', 'contexts', 'figure_file', 'page'});
S.density_name = table_col_to_string(S.density_name);
S.density_file = table_col_to_string(S.density_file);
S.density_index = table_col_to_double(S.density_index);
S.blp_method = table_col_to_string(S.blp_method);
S.blp_k = table_col_to_double(S.blp_k);
S.process_label = table_col_to_string(S.process_label);
S.topN_occurrence_count = table_col_to_double(S.topN_occurrence_count);
S.max_peak_abs_corr = table_col_to_double(S.max_peak_abs_corr);
S.mean_peak_abs_corr = table_col_to_double(S.mean_peak_abs_corr);
S.contexts = table_col_to_string(S.contexts);
S.figure_file = strings(height(S), 1);
S.page = strings(height(S), 1);
S = sortrows(S, {'topN_occurrence_count', 'max_peak_abs_corr'}, {'descend', 'descend'});

per_page = 12;
n_pages = ceil(height(S) / per_page);
for page = 1:n_pages
    i1 = (page - 1) * per_page + 1;
    i2 = min(page * per_page, height(S));
    Sp = S(i1:i2, :);
    fig = figure('Visible', 'off', 'Color', 'w', ...
        'Units', 'pixels', 'Position', [80, 80, 1900, 235 * height(Sp) + 150]);
    tl = tiledlayout(fig, height(Sp), 1, 'TileSpacing', 'compact', 'Padding', 'compact');
    title(tl, sprintf('E10gH1 unique dimred BLP components appearing in top%d xcorr | page %d/%d', ...
        top_n, page, n_pages), 'Interpreter', 'none', 'FontWeight', 'bold');

    for i = 1:height(Sp)
        [x, t] = load_density_series(char(Sp.density_file(i)), double(Sp.density_index(i)));
        nexttile(tl);
        plot(t / 60, x, 'Color', [0.08, 0.32, 0.80], 'LineWidth', 1.0);
        yline(0, '-', 'Color', [0.75, 0.75, 0.75], 'LineWidth', 0.7);
        grid on;
        box off;
        ylabel('z');
        if i == height(Sp)
            xlabel('time (min)');
        end
        title(sprintf('%02d | %s k%02d comp%d | %s | n=%d, max |r|=%.3f, mean |r|=%.3f', ...
            i1 + i - 1, Sp.blp_method(i), double(Sp.blp_k(i)), ...
            double(Sp.density_index(i)), Sp.process_label(i), ...
            double(Sp.topN_occurrence_count(i)), ...
            double(Sp.max_peak_abs_corr(i)), ...
            double(Sp.mean_peak_abs_corr(i))), ...
            'Interpreter', 'none', 'FontSize', 9);
    end

    png_file = fullfile(fig_dir, sprintf('gallery_unique_dimred_components_top%02d_page%02d.png', top_n, page));
    exportgraphics(fig, png_file, 'Resolution', 170);
    close(fig);
    S.figure_file(i1:i2) = string(png_file);
    S.page(i1:i2) = string(sprintf('page%02d', page));
    fprintf('Wrote %s\n', png_file);
end

writetable(S, fullfile(fig_dir, sprintf('e10gh1_unique_dimred_components_top%02d.csv', top_n)));
end


function s = table_col_to_string(col)
if iscell(col)
    s = string(col);
else
    s = string(col);
end
end


function v = table_col_to_double(col)
if iscell(col)
    v = cellfun(@scalar_double, col);
else
    v = double(col);
end
v = v(:);
end


function [x, t] = load_density_series(density_file, density_idx)
persistent density_cache
if isempty(density_cache)
    density_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
key = char(density_file);
if isKey(density_cache, key)
    cached = density_cache(key);
    X = cached.X;
    t = cached.t;
else
    S = load_mat_file_with_short_path(key, 'D', 'density', 't_centers');
    if isfield(S, 'D') && isfield(S.D, 'density_time_by_mode')
        X = double(S.D.density_time_by_mode);
        t = local_time_from_density_struct(S.D, size(X, 1));
    elseif isfield(S, 'D') && isfield(S.D, 'density_time_by_component')
        X = double(S.D.density_time_by_component);
        t = local_time_from_density_struct(S.D, size(X, 1));
    elseif isfield(S, 'density')
        X = double(S.density);
        if isfield(S, 't_centers') && numel(S.t_centers) == size(X, 1)
            t = double(S.t_centers(:));
        else
            t = (0:size(X, 1)-1)' * 2;
        end
    else
        error('Could not find density matrix in %s', key);
    end
    density_cache(key) = struct('X', X, 't', t);
end
if density_idx < 1 || density_idx > size(X, 2)
    error('Density index %d is outside 1:%d for %s', density_idx, size(X, 2), key);
end
x = zscore_vector(X(:, density_idx));
end


function t = local_time_from_density_struct(D, n)
if isfield(D, 't_centers') && numel(D.t_centers) == n
    t = double(D.t_centers(:));
elseif isfield(D, 't_start') && numel(D.t_start) == n
    t = double(D.t_start(:));
else
    t = (0:n-1)' * 2;
end
end


function [y, t, target_file] = load_target_series(row)
persistent target_cache
if isempty(target_cache)
    target_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

pipeline = char(row.pipeline);
target_feature = string(row.bold_feature);
target_idx = resolve_target_index(row);

switch upper(pipeline)
    case 'P8'
        target_file = resolve_p8_bold_post(string(row.observable_mode));
        cache_key = sprintf('P8|%s|%s', target_file, target_feature);
        if isKey(target_cache, cache_key)
            cached = target_cache(cache_key);
            Y = cached.Y;
            t = cached.t;
        else
            S = load_mat_file_with_short_path(target_file, 'BOLD_POST');
            [Y, t] = build_p8_feature_matrix(S.BOLD_POST, target_feature);
            target_cache(cache_key) = struct('Y', Y, 't', t);
        end

    case 'P10'
        target_file = resolve_p10_dimred_result(row);
        cache_key = sprintf('P10|%s|%s', target_file, target_feature);
        if isKey(target_cache, cache_key)
            cached = target_cache(cache_key);
            Y = cached.Y;
            t = cached.t;
        else
            S = load_mat_file_with_short_path(target_file, 'result');
            result = S.result;
            C = double(result.core.temporal_components_time_by_comp);
            if contains(lower(char(target_feature)), "component_abs")
                Y = abs(C);
            else
                Y = real(C);
            end
            t = resolve_result_time(result, size(Y, 1));
            target_cache(cache_key) = struct('Y', Y, 't', t);
        end

    otherwise
        error('Unsupported pipeline: %s', pipeline);
end

if target_idx < 1 || target_idx > size(Y, 2)
    error('Target index %d is outside 1:%d for %s', target_idx, size(Y, 2), target_file);
end
y = zscore_vector(Y(:, target_idx));
end


function idx = resolve_target_index(row)
idx = NaN;
if ismember('bold_component_index', row.Properties.VariableNames)
    idx = scalar_double(row.bold_component_index);
end
if ~isfinite(idx) || idx < 1
    idx = scalar_double(row.bold_mode_index);
end
end


function x = scalar_double(value)
if isstring(value) || ischar(value)
    x = str2double(value);
elseif iscell(value)
    x = str2double(string(value{1}));
else
    x = double(value);
end
if isempty(x)
    x = NaN;
else
    x = x(1);
end
end


function file = resolve_p8_bold_post(observable)
root = fullfile('E:\DataPons_processed', 'e10gh1', ...
    'pipeline7_bold_reskoopnet_postprocessing');
switch lower(char(observable))
    case 'global_svd100'
        patterns = {'*global_svd100_untilstale_20260512*global_svd100_bold_post.mat', ...
            '*global_svd100*bold_post.mat'};
    case 'gsvd100_ds'
        patterns = {'*gsvd100_ds_tune_lr3em6_reg0p01*lowreg*bold_post.mat', ...
            '*gsvd100_ds*bold_post.mat'};
    case 'hp_svd100'
        patterns = {'*HP_svd100_tune_lr1em5_reg0p03*bold_post.mat', ...
            '*HP_svd100*bold_post.mat', '*hp_svd100*bold_post.mat'};
    case 'roi_mean'
        patterns = {'*roi_mean_tune_lr1em6_reg0p003*lowreg*bold_post.mat', ...
            '*roi_mean*bold_post.mat'};
    otherwise
        patterns = {sprintf('*%s*bold_post.mat', char(observable))};
end
file = resolve_latest_matching_file(root, patterns);
end


function file = resolve_p10_dimred_result(row)
source_csv = char(row.source_csv);
p10_run_dir = fileparts(fileparts(source_csv));
[~, p10_tag] = fileparts(p10_run_dir);
parts = split(string(p10_tag), "__");
if numel(parts) >= 3
    run_tag = parts(1);
    p10_feature = parts(2);
    p10_method_tag = parts(3);
else
    run_tag = observable_to_p7_run_tag(string(row.observable_mode));
    p10_feature = p10_feature_from_row(row);
    p10_method_tag = string(row.p10_method_k);
end

root = fullfile('E:\DataPons_processed', 'e10gh1', ...
    'pipeline9_bold_eigenfunction_reduction', ...
    char(run_tag), char(p10_feature), char(p10_method_tag), 'mat');
patterns = {sprintf('e10gh1_%s_%s_%s.mat', ...
    char(run_tag), char(p10_feature), char(p10_method_tag)), '*.mat'};
file = resolve_latest_matching_file(root, patterns);
end


function run_tag = observable_to_p7_run_tag(observable)
switch lower(char(observable))
    case 'global_svd100'
        run_tag = "pv_gsvd100";
    case 'gsvd100_ds'
        run_tag = "pv_gsvd100_ds";
    case 'hp_svd100'
        run_tag = "pv_hp100";
    case 'roi_mean'
        run_tag = "pv_roi";
    otherwise
        run_tag = string(observable);
end
end


function p10_feature = p10_feature_from_row(row)
feature = string(row.bold_feature);
if startsWith(feature, "deconv")
    p10_feature = "deconv_real";
elseif startsWith(feature, "efun")
    p10_feature = "efun_real";
else
    switch lower(char(row.feature_family))
        case 'deconv_efun'
            p10_feature = "deconv_real";
        otherwise
            p10_feature = "efun_real";
    end
end
end


function file = resolve_latest_matching_file(root, patterns)
if exist(root, 'dir') ~= 7
    error('Missing root: %s', root);
end
for i = 1:numel(patterns)
    L = dir(fullfile(root, '**', patterns{i}));
    L = L(~[L.isdir]);
    if isempty(L)
        continue;
    end
    [~, idx] = max([L.datenum]);
    file = fullfile(L(idx).folder, L(idx).name);
    return;
end
error('No matching file under %s', root);
end


function [Y, t] = build_p8_feature_matrix(B, feature_name)
E = B.EDMD_outputs;
switch lower(char(feature_name))
    case 'efun_abs'
        Y = abs(E.efuns);
    case 'efun_real'
        Y = real(E.efuns);
    case 'deconv_abs'
        Y = abs(get_deconv_matrix(E));
    case 'deconv_real'
        Y = real(get_deconv_matrix(E));
    otherwise
        error('Unsupported P8 feature: %s', feature_name);
end
Y = double(Y);
if isfield(B, 'time_vec') && numel(B.time_vec) == size(Y, 1)
    t = double(B.time_vec(:));
elseif isfield(B, 'dt') && ~isempty(B.dt)
    t = (0:size(Y, 1)-1)' * double(B.dt(1));
elseif isfield(E, 'dt') && ~isempty(E.dt)
    t = (0:size(Y, 1)-1)' * double(E.dt(1));
else
    t = (0:size(Y, 1)-1)' * 2;
end
end


function U = get_deconv_matrix(E)
U = [];
if isfield(E, 'deconv_efuns') && isstruct(E.deconv_efuns)
    if isfield(E.deconv_efuns, 'u_sel') && ~isempty(E.deconv_efuns.u_sel)
        U = E.deconv_efuns.u_sel;
    elseif isfield(E.deconv_efuns, 'u_all') && ~isempty(E.deconv_efuns.u_all)
        U = E.deconv_efuns.u_all;
    end
end
if isempty(U)
    error('BOLD EDMD output has no deconv_efuns matrix.');
end
end


function t = resolve_result_time(result, n)
dt = 2;
if isfield(result, 'input') && isfield(result.input, 'dt') && ~isempty(result.input.dt)
    dt = double(result.input.dt(1));
elseif isfield(result, 'run_info') && isfield(result.run_info, 'dt') && ~isempty(result.run_info.dt)
    dt = double(result.run_info.dt(1));
end
if isfield(result, 'input') && isfield(result.input, 'time_axis') && ...
        numel(result.input.time_axis) == n
    t = double(result.input.time_axis(:));
else
    t = (0:n-1)' * dt;
end
end


function [t, x, y] = align_by_lag(x_raw, y_raw, tx, ty, lag_sec)
n = min([numel(x_raw), numel(y_raw), numel(tx), numel(ty)]);
x_raw = x_raw(1:n);
y_raw = y_raw(1:n);
tx = tx(1:n);
ty = ty(1:n); %#ok<NASGU>
dt = median(diff(tx(isfinite(tx))));
if ~isfinite(dt) || dt <= 0
    dt = 2;
end
lag_bins = round(lag_sec / dt);
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
t = tx(ix);
x = x_raw(ix);
y = y_raw(iy);
valid = isfinite(t) & isfinite(x) & isfinite(y);
t = t(valid);
x = x(valid);
y = y(valid);
end


function y = zscore_vector(x)
x = double(x(:));
mu = mean(x, 'omitnan');
sig = std(x, 0, 'omitnan');
if ~isfinite(sig) || sig == 0
    y = nan(size(x));
else
    y = (x - mu) ./ sig;
end
end
