function plot_e10gb1_raw_csplit_top5_xcorr_time_series()
% Plot top-5 raw csplit density xcorr pairs as aligned time series.
%
% This is a diagnostic figure for the E10gb1 standardized-vs-nonstandard
% csplit probe. Each output figure shows the global top-5 xcorr rows for one
% condition x pipeline x feature-family group.

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(fullfile(repo_root, 'functions')));

result_dir = fullfile(repo_root, 'results', ...
    'e10gb1_standardized_csplit_p8_p10_probe');
long_csv = fullfile(result_dir, ...
    'raw_csplit_xcorr_distribution_std_vs_nonstandard_long.csv');

fig_dir = fullfile('E:\DataPons_processed', 'summary_figures', ...
    'pipeline11_current_analysis_summary', ...
    'p8_p10_e10gb1_standardized_csplit_probe', ...
    'raw_csplit_top5_xcorr_timeseries');
if exist(fig_dir, 'dir') ~= 7
    mkdir(fig_dir);
end

T = readtable(long_csv, 'TextType', 'string', 'Delimiter', ',');
groups = unique(T(:, {'condition', 'pipeline', 'feature_family'}), 'rows');

summary_rows = {};
for i_group = 1:height(groups)
    condition = groups.condition(i_group);
    pipeline = groups.pipeline(i_group);
    family = groups.feature_family(i_group);

    mask = T.condition == condition & ...
        T.pipeline == pipeline & ...
        T.feature_family == family & ...
        isfinite(T.peak_abs_corr);
    Tg = sortrows(T(mask, :), 'peak_abs_corr', 'descend');
    if isempty(Tg)
        continue;
    end
    Tg = Tg(1:min(5, height(Tg)), :);

    fprintf('Plotting %s %s %s (%d rows)\n', ...
        condition, pipeline, family, height(Tg));
    [png_file, rows_i] = plot_one_group(Tg, condition, pipeline, family, fig_dir);
    summary_rows = [summary_rows; rows_i]; %#ok<AGROW>
    fprintf('  %s\n', png_file);
end

if ~isempty(summary_rows)
    S = cell2table(summary_rows, 'VariableNames', { ...
        'condition', 'pipeline', 'feature_family', 'rank', ...
        'source_file', 'density_file', 'target_file', ...
        'density_index', 'target_feature', 'target_index', ...
        'peak_abs_corr', 'peak_corr', 'peak_lag_sec', ...
        'figure_file'});
    writetable(S, fullfile(result_dir, 'raw_csplit_top5_xcorr_timeseries_rows.csv'));
end
end


function [png_file, summary_rows] = plot_one_group(Tg, condition, pipeline, family, fig_dir)
n = height(Tg);
fig = figure('Visible', 'off', 'Color', 'w', ...
    'Units', 'pixels', 'Position', [100, 100, 1800, 280 * n + 170]);
tl = tiledlayout(fig, n, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
title(tl, sprintf('E10gb1 raw csplit top5 xcorr time series | %s | %s | %s', ...
    condition, pipeline, family), ...
    'Interpreter', 'none', 'FontWeight', 'bold');

summary_rows = cell(n, 14);
for i = 1:n
    row = Tg(i, :);
    source_csv = normalize_path(row.source_file);
    peak_table = readtable(source_csv, 'TextType', 'string', 'Delimiter', ',');
    rank_in_source = double(row.rank_in_source_csv);
    peak_row = peak_table(rank_in_source, :);

    density_file = normalize_path(peak_row.density_file);
    density_idx = double(peak_row.density_index);
    target_feature = string(peak_row.bold_feature);
    target_idx = double(peak_row.bold_mode_index);
    peak_corr = double(peak_row.peak_corr);
    peak_abs = double(peak_row.peak_abs_corr);
    lag_sec = double(peak_row.peak_lag_sec);

    [density_z, density_t, density_label] = load_density_series(density_file, density_idx);
    [target_z, target_t, target_label, target_file] = ...
        load_target_series(source_csv, pipeline, target_feature, target_idx);

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
    legend({'raw BLP efun density', 'BOLD-side target shifted to peak lag'}, ...
        'Location', 'northeastoutside', 'Box', 'off');
    title(sprintf(['#%d  |r|=%.3f, r=%.3f, lag=%+.1fs | %s idx %d vs %s idx %d\n' ...
        '%s'], ...
        i, peak_abs, peak_corr, lag_sec, density_label, density_idx, ...
        target_label, target_idx, source_csv), ...
        'Interpreter', 'none', 'FontSize', 9);

    summary_rows(i, :) = { ...
        char(condition), char(pipeline), char(family), i, ...
        source_csv, density_file, target_file, ...
        density_idx, char(target_feature), target_idx, ...
        peak_abs, peak_corr, lag_sec, ''};
end

safe_name = sprintf('raw_csplit_top5_xcorr_timeseries__%s__%s__%s.png', ...
    condition, pipeline, family);
safe_name = regexprep(char(safe_name), '[^A-Za-z0-9_.-]', '_');
png_file = fullfile(fig_dir, safe_name);
exportgraphics(fig, png_file, 'Resolution', 170);
close(fig);
for i = 1:n
    summary_rows{i, 14} = png_file;
end
end


function [x, t, label] = load_density_series(density_file, density_idx)
persistent density_cache
if isempty(density_cache)
    density_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
if isKey(density_cache, density_file)
    cached = density_cache(density_file);
    X = cached.X;
    t = cached.t;
else
    S = load_mat_file_with_short_path(density_file, 'D', 'density', 't_centers');
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
        error('Could not find density matrix in %s', density_file);
    end
    density_cache(density_file) = struct('X', X, 't', t);
end
if density_idx < 1 || density_idx > size(X, 2)
    error('Density index %d is outside 1:%d for %s', ...
        density_idx, size(X, 2), density_file);
end
x = zscore_vector(X(:, density_idx));
if isempty(t)
    t = (0:numel(x)-1)' * 2;
end
label = 'raw density';
end


function [y, t, label, target_file] = load_target_series(source_csv, pipeline, target_feature, target_idx)
target_key = sprintf('%s|%s|%s', char(pipeline), extract_run_dir(source_csv), char(target_feature));
persistent target_cache
if isempty(target_cache)
    target_cache = containers.Map('KeyType', 'char', 'ValueType', 'any');
end
if isKey(target_cache, target_key)
    cached = target_cache(target_key);
    Y = cached.Y;
    t = cached.t;
    label = cached.label;
    target_file = cached.target_file;
else
run_dir = extract_run_dir(source_csv);
switch upper(char(pipeline))
    case 'P8'
        run_tag = string(get_last_path_part(run_dir));
        target_file = resolve_p8_bold_post(run_tag);
        S = load_mat_file_with_short_path(target_file, 'BOLD_POST');
        if isfield(S, 'BOLD_POST')
            B = S.BOLD_POST;
        else
            error('BOLD_POST variable missing from %s', target_file);
        end
        [Y, t] = build_p8_feature_matrix(B, target_feature);
        label = 'BOLD ' + target_feature;

    case 'P10'
        run_tag = string(get_last_path_part(run_dir));
        target_file = resolve_p10_dimred_result(run_tag);
        S = load_mat_file_with_short_path(target_file, 'result');
        result = S.result;
        C = double(result.core.temporal_components_time_by_comp);
        if contains(lower(target_feature), "component_abs")
            Y = abs(C);
        else
            Y = real(C);
        end
        t = resolve_result_time(result, size(Y, 1));
        label = 'P9 ' + target_feature;

    otherwise
        error('Unsupported pipeline: %s', pipeline);
end
    target_cache(target_key) = struct( ...
        'Y', Y, 't', t, 'label', label, 'target_file', target_file);
end

if target_idx < 1 || target_idx > size(Y, 2)
    error('Target index %d is outside 1:%d for %s', ...
        target_idx, size(Y, 2), target_file);
end
y = zscore_vector(Y(:, target_idx));
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
ty = ty(1:n);

dt = median(diff(tx(isfinite(tx))));
if ~isfinite(dt) || dt <= 0
    dt = median(diff(ty(isfinite(ty))));
end
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


function file = resolve_p8_bold_post(run_tag)
root = fullfile('E:\DataPons_processed', 'e10gb1', ...
    'pipeline7_bold_reskoopnet_postprocessing');
switch lower(char(run_tag))
    case 'pv_gsvd100'
        patterns = {'*global_svd100_untilstale_20260511*global_svd100_bold_post.mat'};
    case 'pv_gsvd100_ds'
        patterns = {'*gsvd100_ds_tune_lr3em6_reg0p01*b*lowreg*bold_post.mat', ...
            '*gsvd100_ds_tune_lr3em6_reg0p01*lowreg*bold_post.mat', ...
            '*gsvd100_ds*bold_post.mat'};
    case {'pv_roi', 'pv_roi_mean'}
        patterns = {'*roi_mean_tune_lr3em6_reg0p01*lowreg*bold_post.mat', ...
            '*roi_mean*bold_post.mat'};
    case {'pv_hp100', 'pv_hp_svd100'}
        patterns = {'*HP_svd100*bold_post.mat', '*hp_svd100*bold_post.mat'};
    otherwise
        patterns = {sprintf('*%s*bold_post.mat', char(run_tag))};
end
file = resolve_latest_matching_file(root, patterns);
end


function file = resolve_p10_dimred_result(run_tag)
parts = split(string(run_tag), "__");
if numel(parts) < 3
    error('Cannot parse P10 run tag: %s', run_tag);
end
obs = parts(1);
feature = parts(2);
method_tag = parts(3);
root = fullfile('E:\DataPons_processed', 'e10gb1', ...
    'pipeline9_bold_eigenfunction_reduction', ...
    char(obs), char(feature), char(method_tag), 'mat');
patterns = {sprintf('e10gb1_%s_%s_%s.mat', char(obs), char(feature), char(method_tag)), ...
    '*.mat'};
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


function run_dir = extract_run_dir(source_csv)
source_csv = char(source_csv);
marker = [filesep 'feature' filesep];
idx = strfind(source_csv, marker);
if isempty(idx)
    error('Cannot extract run dir from %s', source_csv);
end
run_dir = source_csv(1:idx(1)-1);
end


function part = get_last_path_part(path_in)
parts = split(string(path_in), filesep);
parts = parts(parts ~= "");
part = char(parts(end));
end


function path_out = normalize_path(path_in)
path_out = char(string(path_in));
if startsWith(path_out, '/mnt/e/')
    path_out = ['E:\', strrep(extractAfter(path_out, strlength('/mnt/e/')), '/', '\')];
elseif startsWith(path_out, '/mnt/d/')
    path_out = ['D:\', strrep(extractAfter(path_out, strlength('/mnt/d/')), '/', '\')];
else
    path_out = strrep(path_out, '/', '\');
end
end
