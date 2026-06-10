function [fig, plot_info] = plot_bold_observable_efun_dimred_summary(bold_post_input, result_input, params)
%PLOT_BOLD_OBSERVABLE_EFUN_DIMRED_SUMMARY Plot observables, BOLD efuns, and P9 components.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(result_input)
    error('result_input is required.');
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

B = local_load_bold_post(bold_post_input);
result = local_load_result(result_input);

[X_obs, obs_labels, obs_source] = local_observable_matrix(B);
X_src = result.data.bold_feature_time_by_mode;
X_comp = local_component_matrix(result);
t = result.input.time_axis(:);

T = min([numel(t), size(X_obs, 1), size(X_src, 1), size(X_comp, 1)]);
if T < 2
    error('Not enough aligned samples for P9 summary plot.');
end
t = t(1:T);
X_obs = X_obs(1:T, :);
X_src = X_src(1:T, :);
X_comp = X_comp(1:T, :);

idx = local_plot_indices(T, params);
t_plot = t(idx);

obs_idx = local_select_columns(X_obs(idx, :), params.max_observables, params.select_observables_by);
src_idx = local_select_columns(X_src(idx, :), params.max_source_modes, params.select_modes_by);
comp_idx = 1:min(size(X_comp, 2), params.max_components);

Z_obs = local_zscore_columns(X_obs(idx, obs_idx));
Z_src = local_zscore_columns(X_src(idx, src_idx));
Z_comp = local_zscore_columns(X_comp(idx, comp_idx));

fig = figure('Color', 'w', 'Visible', params.figure_visible, ...
    'Position', params.figure_position);

ax1 = subplot(3, 1, 1);
imagesc(ax1, t_plot, 1:numel(obs_idx), Z_obs.');
set(ax1, 'YDir', 'normal');
title(ax1, sprintf('BOLD observables | %s | %d/%d shown', ...
    obs_source, numel(obs_idx), size(X_obs, 2)), 'Interpreter', 'none');
ylabel(ax1, 'observable');
local_apply_axis_style(ax1, params);
local_draw_session_borders(ax1, B, t_plot);

ax2 = subplot(3, 1, 2);
imagesc(ax2, t_plot, 1:numel(src_idx), Z_src.');
set(ax2, 'YDir', 'normal');
title(ax2, sprintf('%s | %s | %d/%d modes shown', ...
    result.feature.name, result.feature.normalization, ...
    numel(src_idx), size(X_src, 2)), 'Interpreter', 'none');
ylabel(ax2, 'mode');
local_apply_axis_style(ax2, params);
local_draw_session_borders(ax2, B, t_plot);

ax3 = subplot(3, 1, 3);
imagesc(ax3, t_plot, 1:numel(comp_idx), Z_comp.');
set(ax3, 'YDir', 'normal');
title(ax3, sprintf('Reduced BOLD eigenfunction components | %s/%s | k=%d', ...
    result.meta.path_kind, result.meta.method, size(X_comp, 2)), 'Interpreter', 'none');
ylabel(ax3, 'component');
xlabel(ax3, 'time (s)');
local_apply_axis_style(ax3, params);
local_draw_session_borders(ax3, B, t_plot);

linkaxes([ax1, ax2, ax3], 'x');
sgtitle(local_figure_title(B, result), 'Interpreter', 'none');

plot_info = struct();
plot_info.observable_source = obs_source;
plot_info.observable_indices = obs_idx(:);
plot_info.observable_labels = local_labels_at(obs_labels, obs_idx);
plot_info.source_mode_indices = src_idx(:);
plot_info.component_indices = comp_idx(:);
plot_info.window_idx = idx(:);
plot_info.time_start = t_plot(1);
plot_info.time_end = t_plot(end);
plot_info.png_file = '';
plot_info.fig_file = '';

if params.save_png || params.save_fig
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end
end
if params.save_png
    plot_info.png_file = fullfile(params.save_dir, [params.save_tag, '.png']);
    exportgraphics(fig, plot_info.png_file, 'Resolution', params.resolution);
end
if params.save_fig
    plot_info.fig_file = fullfile(params.save_dir, [params.save_tag, '.fig']);
    savefig(fig, plot_info.fig_file);
end
if params.close_figure
    close(fig);
end
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'max_observables', 40);
params = local_set_default(params, 'max_source_modes', 40);
params = local_set_default(params, 'max_components', 20);
params = local_set_default(params, 'max_samples', 2000);
params = local_set_default(params, 'window_start', 1);
params = local_set_default(params, 'window_idx', []);
params = local_set_default(params, 'select_observables_by', 'variance');
params = local_set_default(params, 'select_modes_by', 'order');
params = local_set_default(params, 'colormap', 'turbo');
params = local_set_default(params, 'figure_position', [80, 80, 1280, 820]);
params = local_set_default(params, 'figure_visible', 'off');
params = local_set_default(params, 'save_dir', pwd);
params = local_set_default(params, 'save_tag', 'bold_efun_dimred_summary');
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'resolution', 220);
params = local_set_default(params, 'close_figure', true);
end


function B = local_load_bold_post(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load_mat_file_with_short_path(file, 'BOLD_POST');
    if ~isfield(S, 'BOLD_POST')
        error('BOLD_POST variable missing in %s.', file);
    end
    B = S.BOLD_POST;
    B.source_file = file;
elseif isstruct(input)
    B = input;
else
    error('bold_post_input must be a path or struct.');
end
end


function result = local_load_result(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load_mat_file_with_short_path(file, 'result');
    if ~isfield(S, 'result')
        error('result variable missing in %s.', file);
    end
    result = S.result;
elseif isstruct(input)
    result = input;
else
    error('result_input must be a path or struct.');
end
end


function [X, labels, source] = local_observable_matrix(B)
X = [];
labels = {};
source = 'BOLD_POST.observable';

if isfield(B, 'observable') && isstruct(B.observable)
    [X, labels] = local_observable_from_struct(B.observable);
end

if isempty(X) && isfield(B, 'observable_file') && ~isempty(B.observable_file) && ...
        exist(B.observable_file, 'file') == 2
    source = B.observable_file;
    S = load_mat_file_with_short_path(B.observable_file);
    [X, labels] = local_observable_from_struct(S);
end

if isempty(X)
    error('Could not resolve BOLD observable matrix from BOLD_POST or observable_file.');
end

X = double(X);
if isempty(labels)
    labels = cellstr(compose("obs%03d", 1:size(X, 2)));
end
end


function [X, labels] = local_observable_from_struct(S)
X = [];
labels = {};
if isfield(S, 'O') && isstruct(S.O) && isfield(S.O, 'data') && ~isempty(S.O.data)
    X = S.O.data;
    labels = local_observable_labels(S.O);
elseif isfield(S, 'obs') && ~isempty(S.obs)
    X = S.obs;
    labels = local_observable_labels(S);
elseif isfield(S, 'data') && ~isempty(S.data)
    X = S.data;
    labels = local_observable_labels(S);
end
end


function labels = local_observable_labels(S)
labels = {};
if isfield(S, 'observable_labels') && ~isempty(S.observable_labels)
    labels = cellstr(string(S.observable_labels(:)));
elseif isfield(S, 'observable_info') && istable(S.observable_info) && ...
        ismember('observable_label', S.observable_info.Properties.VariableNames)
    labels = cellstr(string(S.observable_info.observable_label(:)));
elseif isfield(S, 'obs_info') && istable(S.obs_info) && ...
        ismember('observable_label', S.obs_info.Properties.VariableNames)
    labels = cellstr(string(S.obs_info.observable_label(:)));
end
end


function C = local_component_matrix(result)
if isfield(result, 'summary') && isfield(result.summary, 'temporal_components_smooth_time_by_comp') && ...
        ~isempty(result.summary.temporal_components_smooth_time_by_comp)
    C = result.summary.temporal_components_smooth_time_by_comp;
else
    C = result.core.temporal_components_time_by_comp;
end
end


function idx = local_plot_indices(T, params)
if ~isempty(params.window_idx)
    idx = unique(round(double(params.window_idx(:))), 'stable');
    idx = idx(idx >= 1 & idx <= T);
    if ~isempty(idx)
        return;
    end
end
start_idx = max(1, round(double(params.window_start)));
stop_idx = min(T, start_idx + max(1, round(double(params.max_samples))) - 1);
idx = (start_idx:stop_idx).';
end


function idx = local_select_columns(X, max_cols, mode_name)
n = size(X, 2);
max_cols = min(n, max(1, round(double(max_cols))));
switch lower(char(string(mode_name)))
    case {'variance', 'var'}
        score = var(double(X), 0, 1, 'omitnan');
        score(~isfinite(score)) = -Inf;
        [~, ord] = sort(score, 'descend');
        idx = ord(1:max_cols);
    case {'order', 'first', 'preserve'}
        idx = 1:max_cols;
    otherwise
        error('Unsupported column selection mode: %s.', mode_name);
end
idx = idx(:).';
end


function Z = local_zscore_columns(X)
X = double(X);
mu = mean(X, 1, 'omitnan');
sigma = std(X, 0, 1, 'omitnan');
sigma(~isfinite(sigma) | sigma <= 0) = 1;
Z = (X - mu) ./ sigma;
Z(~isfinite(Z)) = 0;
% Keep the observable traces un-clipped so large excursions remain visible.
end


function local_apply_axis_style(ax, params)
cmap = local_colormap(params.colormap);
colormap(ax, cmap);
colorbar(ax);
axis(ax, 'tight');
set(ax, 'TickDir', 'out', 'Box', 'off');
end


function cmap = local_colormap(name)
if isa(name, 'function_handle')
    cmap = name();
    return;
end
name = char(string(name));
if exist(name, 'file') == 2 || exist(name, 'builtin') == 5
    cmap = feval(name, 256);
else
    cmap = parula(256);
end
end


function local_draw_session_borders(ax, B, t_plot)
if ~isfield(B, 'session') || ~isstruct(B.session) || ...
        ~isfield(B.session, 'border_idx') || isempty(B.session.border_idx)
    return;
end
dt = local_get_field(B, 'dt', []);
if isempty(dt) || ~isfinite(dt) || dt <= 0
    return;
end
borders = double(B.session.border_idx(:));
border_t = (borders - 1) * double(dt);
border_t = border_t(border_t >= min(t_plot) & border_t <= max(t_plot));
for i = 1:numel(border_t)
    xline(ax, border_t(i), '-', 'Color', [0 0 0], 'LineWidth', 0.75);
end
end


function labels = local_labels_at(all_labels, idx)
labels = {};
if isempty(all_labels)
    return;
end
idx = idx(idx >= 1 & idx <= numel(all_labels));
labels = all_labels(idx);
end


function title_text = local_figure_title(B, result)
run_info = local_get_field(B, 'run_info', struct());
dataset = local_get_field(run_info, 'dataset_stem', '');
run_name = local_get_field(run_info, 'run_name', '');
title_text = sprintf('Pipeline 9 | %s | %s | %s | %s', ...
    dataset, run_name, result.feature.name, result.meta.method);
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
