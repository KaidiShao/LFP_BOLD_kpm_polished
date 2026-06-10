function plot_e10gb1_dimred_csplit_top5_xcorr_time_series()
% Plot top-5 dimred BLP density xcorr pairs as aligned time series.
%
% This uses the standardized csplit P8/P10 probe table. The goal is to see
% whether dimred BLP efun densities preserve the same efun-vs-deconv split
% seen in raw BLP efun densities.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'functions')));

result_dir = fullfile(repo_root, 'results', ...
    'e10gb1_standardized_csplit_p8_p10_probe');
input_csv = fullfile(result_dir, 'all_top_rows.csv');

fig_dir = fullfile('E:\DataPons_processed', 'summary_figures', ...
    'pipeline11_current_analysis_summary', ...
    'p8_p10_e10gb1_standardized_csplit_probe', ...
    'dimred_csplit_top5_xcorr_timeseries');
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

T = readtable(input_csv, 'TextType', 'string', 'Delimiter', ',');
T = T(T.source_kind == "dimred" & isfinite(T.peak_abs_corr), :);
groups = unique(T(:, {'analysis', 'feature_family'}), 'rows');

summary_rows = {};
for i_group = 1:height(groups)
    analysis = groups.analysis(i_group);
    family = groups.feature_family(i_group);
    mask = T.analysis == analysis & T.feature_family == family;
    Tg = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
    if isempty(Tg)
        continue;
    end
    Tg = Tg(1:min(5, height(Tg)), :);

    fprintf('Plotting dimred %s %s (%d rows)\n', analysis, family, height(Tg));
    [png_file, rows_i] = plot_one_group(Tg, analysis, family, fig_dir);
    summary_rows = [summary_rows; rows_i]; %#ok<AGROW>
    fprintf('  %s\n', png_file);
end

if ~isempty(summary_rows)
    S = cell2table(summary_rows, 'VariableNames', { ...
        'pipeline', 'feature_family', 'rank', 'observable', ...
        'bold_feature', 'bold_mode_index', 'density_name', ...
        'density_file', 'density_index', 'blp_method', 'blp_k', ...
        'label', 'peak_abs_corr', 'peak_corr', 'peak_lag_sec', ...
        'target_file', 'figure_file'});
    writetable(S, fullfile(result_dir, 'dimred_csplit_top5_xcorr_timeseries_rows.csv'));
end
end


function [png_file, summary_rows] = plot_one_group(Tg, analysis, family, fig_dir)
n = height(Tg);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100, 100, 1850, 285 * n + 175]);
tl = tiledlayout(fig, n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('E10gb1 dimred csplit top5 xcorr time series | %s | %s', ...
    analysis, family), 'Interpreter', 'none', 'FontWeight', 'bold');

summary_rows = cell(n, 17);
for i = 1:n
    row = Tg(i, :);
    density_file = char(row.density_file);
    density_idx = double(row.density_index);
    target_feature = string(row.bold_feature);
    target_idx = double(row.bold_mode_index);
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
    plot(t_plot / 60, x_plot, 'Color', [0.1, 0.35, 0.85], 'LineWidth', 1.0);
    plot(t_plot / 60, y_plot, 'Color', [0.85, 0.2, 0.1], 'LineWidth', 1.0);
    yline(0, '-', 'Color', [0.72, 0.72, 0.72], 'LineWidth', 0.7);
    grid on;
    box off;
    xlim([min(t_plot / 60), max(t_plot / 60)]);
    ylabel('z');
    if i == n
        xlabel('time (min)');
    end
    legend({'dimred BLP efun density', 'BOLD-side target shifted to peak lag'}, ...
        'Location', 'northeastoutside', 'Box', 'off');
    title(sprintf(['#%d  |r|=%.3f, r=%.3f, lag=%+.1fs | %s k%02d comp%d %s | %s idx %d\n' ...
        '%s'], ...
        i, peak_abs, peak_corr, lag_sec, row.blp_method, double(row.blp_k), ...
        density_idx, label, target_feature, target_idx, density_file), ...
        'Interpreter', 'none', 'FontSize', 9);

    summary_rows(i, :) = { ...
        char(analysis), char(family), i, char(row.observable), ...
        char(target_feature), target_idx, char(row.density_name), ...
        density_file, density_idx, char(row.blp_method), double(row.blp_k), ...
        char(label), peak_abs, peak_corr, lag_sec, target_file, ''};
end

safe_name = sprintf('dimred_csplit_top5_xcorr_timeseries__%s__%s.png', ...
    analysis, family);
safe_name = regexprep(char(safe_name), '[^A-Za-z0-9_.-]', '_');
png_file = fullfile(fig_dir, safe_name);
exportgraphics(fig, png_file, 'Resolution', 170);
close(fig);
for i = 1:n
    summary_rows{i, 17} = png_file;
end
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
        if isfield(S.D, 't_centers') && ~isempty(S.D.t_centers)
            t = double(S.D.t_centers(:));
        else
            t = [];
        end
    elseif isfield(S, 'density')
        X = double(S.density);
        if isfield(S, 't_centers') && ~isempty(S.t_centers)
            t = double(S.t_centers(:));
        else
            t = [];
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
if isempty(t)
    t = (0:numel(x)-1)' * 2;
end
end


function [y, t, target_file] = load_target_series(row)
persistent target_cache
if isempty(target_cache)
    target_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end

pipeline = char(row.analysis);
target_feature = string(row.bold_feature);
target_idx = double(row.bold_mode_index);
switch upper(pipeline)
    case 'P8'
        target_file = resolve_p8_bold_post(string(row.observable));
        cache_key = sprintf('P8|%s|%s', target_file, target_feature);
        if isKey(target_cache, cache_key)
            cached = target_cache(cache_key);
            Y = cached.Y;
            t = cached.t;
        else
            S = load_mat_file_with_short_path(target_file, 'BOLD_POST');
            B = S.BOLD_POST;
            [Y, t] = build_p8_feature_matrix(B, target_feature);
            target_cache(cache_key) = struct('Y', Y, 't', t);
        end

    case 'P10'
        [p10_feature, p10_method_tag] = resolve_p10_context(row);
        target_file = resolve_p10_dimred_result( ...
            string(row.observable), p10_feature, p10_method_tag);
        cache_key = sprintf('P10|%s|%s', target_file, target_feature);
        if isKey(target_cache, cache_key)
            cached = target_cache(cache_key);
            Y = cached.Y;
            t = cached.t;
        else
            S = load_mat_file_with_short_path(target_file, 'result');
            result = S.result;
            C = double(result.core.temporal_components_time_by_comp);
            if contains(lower(target_feature), "component_abs")
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


function [p10_feature, p10_method_tag] = resolve_p10_context(row)
p10_feature = string(row.p10_bold_feature);
p10_method_tag = string(row.p10_bold_method_tag);
if ismissing(p10_feature) || strlength(p10_feature) == 0 || ...
        ismissing(p10_method_tag) || strlength(p10_method_tag) == 0
    parts = split(string(row.context_id), "::");
    if numel(parts) >= 4
        p10_feature = parts(3);
        p10_method_tag = parts(4);
    end
end
if ismissing(p10_feature) || strlength(p10_feature) == 0 || ...
        ismissing(p10_method_tag) || strlength(p10_method_tag) == 0
    error('Cannot resolve P10 context from row context_id=%s', string(row.context_id));
end
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


function file = resolve_p8_bold_post(observable)
root = fullfile('E:\DataPons_processed', 'e10gb1', ...
    'pipeline7_bold_reskoopnet_postprocessing');
switch lower(char(observable))
    case 'pv_gsvd100'
        patterns = {'*global_svd100_untilstale_20260511*global_svd100_bold_post.mat'};
    case 'pv_gsvd100_ds'
        patterns = {'*gsvd100_ds_tune_lr3em6_reg0p01*lowreg*bold_post.mat', ...
            '*gsvd100_ds*bold_post.mat'};
    case {'pv_roi', 'pv_roi_mean'}
        patterns = {'*roi_mean_tune_lr3em6_reg0p01*lowreg*bold_post.mat', ...
            '*roi_mean*bold_post.mat'};
    case {'pv_hp100', 'pv_hp_svd100'}
        patterns = {'*HP_svd100*bold_post.mat', '*hp_svd100*bold_post.mat'};
    otherwise
        patterns = {sprintf('*%s*bold_post.mat', char(observable))};
end
file = resolve_latest_matching_file(root, patterns);
end


function file = resolve_p10_dimred_result(observable, p10_feature, p10_method_tag)
root = fullfile('E:\DataPons_processed', 'e10gb1', ...
    'pipeline9_bold_eigenfunction_reduction', ...
    char(observable), char(p10_feature), char(p10_method_tag), 'mat');
patterns = {sprintf('e10gb1_%s_%s_%s.mat', ...
    char(observable), char(p10_feature), char(p10_method_tag)), '*.mat'};
file = resolve_latest_matching_file(root, patterns);
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
