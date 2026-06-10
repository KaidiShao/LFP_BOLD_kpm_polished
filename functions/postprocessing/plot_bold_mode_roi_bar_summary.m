function [fig, out] = plot_bold_mode_roi_bar_summary(bold_post_input, params)
%PLOT_BOLD_MODE_ROI_BAR_SUMMARY Plot ROI-averaged Koopman mode bar summaries.
%
% This helper reduces each selected Koopman mode into one value per ROI,
% then plots one horizontal bar chart per selected basis.
% The same observable-to-voxel mapping contract used by activation maps is
% reused here so voxel, roi_mean, svd/global_svd100, and multi-feature
% observables follow one consistent path.
%
% Required input:
%   bold_post_input : path to *_bold_post.mat, an in-memory BOLD_POST
%                     struct, or a prebuilt plot_ctx from
%                     build_bold_activation_plot_context
%
% Useful params:
%   selection_mode      : 'sorted' or 'raw'
%   basis_indices       : vector of Koopman basis indices to visualize
%   roi_value_mode      : 'mean_abs', 'abs_mean', 'real_mean', 'imag_mean'
%   roi_reduce          : 'mean', 'median', or 'sum'
%   mode_normalization  : 'range', 'maxabs', or 'none'
%   flip_roi_order      : true/false, flips ROI order top-to-bottom
%   layout_mode         : 'row' for 1 x n_basis, 'column' for n_basis x 1
%   show_roi_labels     : 'all', 'first', 'last', or 'none'
%   mode_title_texts    : optional cellstr to override per-panel titles
%   figure_title_text   : optional string to override the shared figure title
%   highlight_regions   : true/false, highlight selected ROI rows
%   highlight_exact_names / highlight_contains_tokens
%   show_cortical_subcortical_separator : true/false

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

if local_is_plot_context(bold_post_input)
    plot_ctx = bold_post_input;
    B = struct();
else
    [plot_ctx, B] = build_bold_activation_plot_context(bold_post_input, struct( ...
        'datapons_root', params.datapons_root, ...
        'roi_ts_file', params.roi_ts_file, ...
        'feature_reduce', params.feature_reduce));
end

[raw_indices, selection_labels] = local_resolve_raw_indices( ...
    plot_ctx.evalues, params.selection_mode, params.basis_indices);
if isempty(raw_indices)
    error('No Koopman basis indices were resolved for plotting.');
end

n_mode = numel(raw_indices);
n_roi = numel(plot_ctx.region_labels);
roi_values_native = nan(n_mode, n_roi);
map_infos = cell(n_mode, 1);

for i_mode = 1:n_mode
    raw_idx = raw_indices(i_mode);
    mode_obs = double(plot_ctx.kpm_modes(raw_idx, :));
    if ~plot_ctx.direct_voxel_mode
        source_values = mode_obs * plot_ctx.coeff_t;
    else
        source_values = mode_obs;
    end

    [voxel_values, map_info] = map_bold_source_values_to_voxels( ...
        source_values(:), plot_ctx, struct('feature_reduce', params.feature_reduce));
    map_infos{i_mode} = map_info;
    roi_values_native(i_mode, :) = local_reduce_voxels_to_regions( ...
        voxel_values(:), plot_ctx.voxel_region_idx(:), n_roi, ...
        params.roi_value_mode, params.roi_reduce);
end

roi_values_plot = local_normalize_mode_rows(roi_values_native, params.mode_normalization);
plot_order = 1:n_roi;
if params.flip_roi_order
    plot_order = fliplr(plot_order);
end
roi_labels_plot = cellstr(string(plot_ctx.region_labels(plot_order)));
roi_values_plot = roi_values_plot(:, plot_order);
highlight_spec = local_build_highlight_spec(roi_labels_plot, params);
separator_spec = local_build_group_separator_spec(roi_labels_plot, params);
roi_labels_display = local_format_roi_labels(roi_labels_plot, highlight_spec, params);

mode_colors = local_resolve_mode_colors(n_mode, params);
fig = figure('Color', 'w', 'Name', 'BOLD ROI-average Koopman mode summary');
if isprop(fig, 'ToolBar')
    fig.ToolBar = 'none';
end
if isprop(fig, 'MenuBar')
    fig.MenuBar = 'none';
end
[n_rows, n_cols, show_mode_normalized] = local_resolve_layout(n_mode, params.layout_mode, params.show_roi_labels);
[fig_width, fig_height] = local_resolve_figure_size(n_mode, n_roi, params.layout_mode);
fig.Position = [60 40 fig_width fig_height];
tl = tiledlayout(n_rows, n_cols, 'TileSpacing', 'compact', 'Padding', 'compact');

axes_list = gobjects(n_mode, 1);
for i_mode = 1:n_mode
    raw_idx = raw_indices(i_mode);
    eval_i = plot_ctx.evalues(raw_idx);
    axes_list(i_mode) = nexttile;
    ax = axes_list(i_mode);
    x_limits_i = local_x_limits(roi_values_plot(i_mode, :), params.mode_normalization, params.roi_value_mode);
    hold(ax, 'on');
    local_draw_highlight_rows(ax, x_limits_i, highlight_spec);
    barh(ax, 1:n_roi, roi_values_plot(i_mode, :), 0.78, ...
        'FaceColor', mode_colors(i_mode, :), ...
        'EdgeColor', 'none', ...
        'FaceAlpha', 0.92);
    local_draw_group_separator(ax, x_limits_i, separator_spec);
    hold(ax, 'off');
    ax.Box = 'off';
    ax.YDir = 'reverse';
    ax.YTick = 1:n_roi;
    ax.YTickLabel = roi_labels_display;
    ax.YGrid = 'on';
    ax.XGrid = 'on';
    ax.GridAlpha = 0.18;
    ax.Layer = 'top';
    ax.FontSize = params.axis_font_size;
    ax.YAxis.FontSize = params.roi_label_font_size;
    ax.TickLabelInterpreter = 'none';
    if isprop(ax, 'Toolbar') && ~isempty(ax.Toolbar)
        ax.Toolbar.Visible = 'off';
    end

    if ~local_should_show_roi_labels(show_mode_normalized, i_mode, n_mode)
        ax.YTickLabel = repmat({''}, n_roi, 1);
    end

    xlim(ax, x_limits_i);
    ylabel(ax, '');

    title(ax, local_mode_title_text(params, selection_labels{i_mode}, raw_idx, eval_i, i_mode), ...
        'Interpreter', 'tex', 'FontWeight', 'bold');
end

linkaxes(axes_list, 'y');
xlabel(tl, local_value_label(params), 'Interpreter', 'none');
sgtitle(tl, local_figure_title_text(params, plot_ctx), 'Interpreter', 'none');

out = struct();
out.run_info = plot_ctx.run_info;
out.observable_file = local_get_field(plot_ctx, 'observable_file', '');
out.roi_ts_file = plot_ctx.roi_ts_file;
out.selection_mode = params.selection_mode;
out.basis_indices = params.basis_indices(:).';
out.selection_labels = selection_labels;
out.raw_indices = raw_indices;
out.roi_labels_native = cellstr(string(plot_ctx.region_labels(:)));
out.roi_labels_plot = roi_labels_plot;
out.roi_labels_display = roi_labels_display;
out.plot_order = plot_order;
out.roi_values_native = roi_values_native;
out.roi_values_plot = roi_values_plot;
out.roi_value_mode = params.roi_value_mode;
out.roi_reduce = params.roi_reduce;
out.mode_normalization = params.mode_normalization;
out.flip_roi_order = params.flip_roi_order;
out.map_space = local_describe_map_space(plot_ctx);
out.source_space_kind = plot_ctx.source_space_kind;
out.mapping_summary = plot_ctx.mapping_summary;
out.map_infos = map_infos;
out.highlight_spec = highlight_spec;
out.separator_spec = separator_spec;
out.post_file = local_get_field(B, 'source_file', '');
end


function params = local_apply_defaults(params)
partition_defaults = build_bold_cortical_subcortical_partition_defaults();
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'roi_ts_file', '');
params = local_set_default(params, 'selection_mode', 'sorted');
params = local_set_default(params, 'basis_indices', 1:3);
params = local_set_default(params, 'feature_reduce', 'mean');
params = local_set_default(params, 'roi_reduce', 'mean');
params = local_set_default(params, 'roi_value_mode', 'mean_abs');
params = local_set_default(params, 'mode_normalization', 'range');
params = local_set_default(params, 'flip_roi_order', true);
params = local_set_default(params, 'layout_mode', 'row');
if ~isfield(params, 'show_roi_labels') || isempty(params.show_roi_labels)
    if strcmpi(char(string(params.layout_mode)), 'row')
        params.show_roi_labels = 'first';
    else
        params.show_roi_labels = 'all';
    end
end
params = local_set_default(params, 'highlight_regions', false);
params = local_set_default(params, 'highlight_exact_names', {});
params = local_set_default(params, 'highlight_contains_tokens', {});
params = local_set_default(params, 'highlight_row_color', [0.10 0.10 0.10]);
params = local_set_default(params, 'highlight_row_alpha', 0.08);
params = local_set_default(params, 'highlight_label_prefix', '* ');
params = local_set_default(params, 'show_cortical_subcortical_separator', false);
params = local_set_default(params, 'separator_after_roi_name', partition_defaults.separator_after_roi_name);
params = local_set_default(params, 'separator_line_color', [0.22 0.22 0.22]);
params = local_set_default(params, 'separator_line_width', 1.3);
params = local_set_default(params, 'separator_line_style', '-');
params = local_set_default(params, 'cortical_exact_names', partition_defaults.cortical_exact_names);
params = local_set_default(params, 'cortical_contains_tokens', partition_defaults.cortical_contains_tokens);
params = local_set_default(params, 'subcortical_exact_names', partition_defaults.subcortical_exact_names);
params = local_set_default(params, 'subcortical_contains_tokens', partition_defaults.subcortical_contains_tokens);
params = local_set_default(params, 'axis_font_size', 10);
params = local_set_default(params, 'roi_label_font_size', 7);
params = local_set_default(params, 'mode_colors', []);
params = local_set_default(params, 'mode_title_texts', {});
params = local_set_default(params, 'figure_title_text', '');
end


function [raw_indices, selection_labels] = local_resolve_raw_indices(evalues, selection_mode, basis_indices)
basis_indices = unique(round(double(basis_indices(:).')), 'stable');
basis_indices = basis_indices(isfinite(basis_indices) & basis_indices >= 1);
selection_labels = {};
if isempty(basis_indices)
    raw_indices = [];
    return;
end

switch lower(char(string(selection_mode)))
    case 'sorted'
        [~, order] = sort(abs(evalues), 'descend');
        valid = basis_indices(basis_indices <= numel(order));
        raw_indices = order(valid);
        selection_labels = arrayfun(@(x) sprintf('sorted%03d', x), valid, 'UniformOutput', false);
    case 'raw'
        valid = basis_indices(basis_indices <= numel(evalues));
        raw_indices = valid;
        selection_labels = arrayfun(@(x) sprintf('raw%03d', x), valid, 'UniformOutput', false);
    otherwise
        error('selection_mode must be ''sorted'' or ''raw''.');
end

raw_indices = raw_indices(:).';
end


function roi_values = local_reduce_voxels_to_regions(voxel_values, voxel_region_idx, n_roi, roi_value_mode, roi_reduce)
roi_values = nan(1, n_roi);
for i_roi = 1:n_roi
    vals = voxel_values(voxel_region_idx == i_roi);
    vals = vals(isfinite(real(vals)) & isfinite(imag(vals)));
    if isempty(vals)
        continue;
    end

    switch lower(char(string(roi_value_mode)))
        case 'abs_mean'
            roi_values(i_roi) = abs(local_reduce(vals, roi_reduce));
        case 'mean_abs'
            roi_values(i_roi) = local_reduce(abs(vals), roi_reduce);
        case 'real_mean'
            roi_values(i_roi) = local_reduce(real(vals), roi_reduce);
        case 'imag_mean'
            roi_values(i_roi) = local_reduce(imag(vals), roi_reduce);
        otherwise
            error(['roi_value_mode must be ''abs_mean'', ''mean_abs'', ' ...
                '''real_mean'', or ''imag_mean''.']);
    end
end
end


function value = local_reduce(x, reduce_name)
switch lower(char(string(reduce_name)))
    case 'mean'
        value = mean(x, 'omitnan');
    case 'median'
        value = median(x, 'omitnan');
    case 'sum'
        value = sum(x, 'omitnan');
    otherwise
        error('roi_reduce must be ''mean'', ''median'', or ''sum''.');
end
end


function X = local_normalize_mode_rows(X, normalization_name)
switch lower(char(string(normalization_name)))
    case 'none'
        return;
    case 'range'
        for i_row = 1:size(X, 1)
            x = X(i_row, :);
            mn = min(x, [], 'omitnan');
            mx = max(x, [], 'omitnan');
            if isfinite(mn) && isfinite(mx) && mx > mn
                X(i_row, :) = (x - mn) ./ (mx - mn);
            elseif isfinite(mx) && abs(mx) > 0
                X(i_row, :) = double(isfinite(x));
            else
                X(i_row, :) = zeros(size(x));
            end
        end
    case 'maxabs'
        for i_row = 1:size(X, 1)
            x = X(i_row, :);
            mx = max(abs(x), [], 'omitnan');
            if isfinite(mx) && mx > 0
                X(i_row, :) = x ./ mx;
            else
                X(i_row, :) = zeros(size(x));
            end
        end
    otherwise
        error('mode_normalization must be ''none'', ''range'', or ''maxabs''.');
end
end


function colors = local_resolve_mode_colors(n_mode, params)
if ~isempty(params.mode_colors)
    colors_in = double(params.mode_colors);
    if size(colors_in, 2) == 3 && size(colors_in, 1) >= n_mode
        colors = colors_in(1:n_mode, :);
        return;
    end
end

base = lines(max(n_mode, 3));
base = 0.15 + 0.85 * base;
colors = base(1:n_mode, :);
end


function tf = local_is_plot_context(input)
tf = isstruct(input) && ~isfield(input, 'EDMD_outputs') && ...
    isfield(input, 'run_info') && ...
    isfield(input, 'kpm_modes') && ...
    isfield(input, 'evalues') && ...
    isfield(input, 'region_labels') && ...
    isfield(input, 'voxel_region_idx');
end


function spec = local_build_highlight_spec(roi_labels_plot, params)
spec = struct();
spec.enabled = logical(local_get_field(params, 'highlight_regions', false));
spec.color = double(local_get_field(params, 'highlight_row_color', [0.10 0.10 0.10]));
spec.alpha = double(local_get_field(params, 'highlight_row_alpha', 0.08));
spec.exact_names = cellstr(string(local_get_field(params, 'highlight_exact_names', {})));
spec.contains_tokens = cellstr(string(local_get_field(params, 'highlight_contains_tokens', {})));
spec.mask = local_match_roi_names(roi_labels_plot, spec.exact_names, spec.contains_tokens);
spec.indices = find(spec.mask);
spec.labels = roi_labels_plot(spec.mask);
if ~spec.enabled || isempty(spec.indices)
    spec.enabled = false;
end
end


function spec = local_build_group_separator_spec(roi_labels_plot, params)
spec = struct();
spec.enabled = logical(local_get_field(params, 'show_cortical_subcortical_separator', false));
spec.position = [];
spec.after_roi_name = '';
spec.before_roi_name = '';
spec.line_color = double(local_get_field(params, 'separator_line_color', [0.22 0.22 0.22]));
spec.line_width = double(local_get_field(params, 'separator_line_width', 1.3));
spec.line_style = char(string(local_get_field(params, 'separator_line_style', '-')));
if ~spec.enabled
    return;
end

after_roi_name = char(string(local_get_field(params, 'separator_after_roi_name', '')));
if ~isempty(after_roi_name)
    idx = find(strcmp(roi_labels_plot, after_roi_name), 1, 'first');
else
    partition_defaults = build_bold_cortical_subcortical_partition_defaults();
    cortical_exact_names = cellstr(string(local_get_field(params, 'cortical_exact_names', partition_defaults.cortical_exact_names)));
    cortical_contains_tokens = cellstr(string(local_get_field(params, 'cortical_contains_tokens', {})));
    cortical_mask = local_match_roi_names(roi_labels_plot, cortical_exact_names, cortical_contains_tokens);
    idx = find(cortical_mask, 1, 'last');
end

if isempty(idx) || idx < 1 || idx >= numel(roi_labels_plot)
    spec.enabled = false;
    return;
end

spec.position = idx + 0.5;
spec.after_roi_name = roi_labels_plot{idx};
spec.before_roi_name = roi_labels_plot{idx + 1};
end


function labels_out = local_format_roi_labels(roi_labels_plot, highlight_spec, params)
labels_out = cellstr(string(roi_labels_plot));
if ~highlight_spec.enabled
    return;
end

prefix = char(string(local_get_field(params, 'highlight_label_prefix', '* ')));
if isempty(prefix)
    return;
end

for i_idx = highlight_spec.indices(:).'
    labels_out{i_idx} = [prefix labels_out{i_idx}];
end
end


function mask = local_match_roi_names(roi_labels, exact_names, contains_tokens)
roi_labels = cellstr(string(roi_labels(:)));
mask = false(numel(roi_labels), 1);

if nargin >= 2 && ~isempty(exact_names)
    exact_names = cellstr(string(exact_names(:)));
    mask = mask | ismember(roi_labels, exact_names);
end

if nargin >= 3 && ~isempty(contains_tokens)
    contains_tokens = cellstr(string(contains_tokens(:)));
    for i_token = 1:numel(contains_tokens)
        token = char(string(contains_tokens{i_token}));
        if isempty(token)
            continue;
        end
        mask = mask | contains(roi_labels, token, 'IgnoreCase', true);
    end
end
end


function local_draw_highlight_rows(ax, x_limits_i, highlight_spec)
if ~highlight_spec.enabled || isempty(highlight_spec.indices)
    return;
end

x1 = x_limits_i(1);
x2 = x_limits_i(2);
for i_idx = highlight_spec.indices(:).'
    patch(ax, [x1 x2 x2 x1], [i_idx - 0.48 i_idx - 0.48 i_idx + 0.48 i_idx + 0.48], ...
        highlight_spec.color, ...
        'FaceAlpha', highlight_spec.alpha, ...
        'EdgeColor', 'none', ...
        'HandleVisibility', 'off');
end
end


function local_draw_group_separator(ax, x_limits_i, separator_spec)
if ~separator_spec.enabled || isempty(separator_spec.position)
    return;
end

line(ax, x_limits_i, [separator_spec.position separator_spec.position], ...
    'Color', separator_spec.line_color, ...
    'LineWidth', separator_spec.line_width, ...
    'LineStyle', separator_spec.line_style, ...
    'HandleVisibility', 'off');
end


function tf = local_should_show_roi_labels(show_mode, i_mode, n_mode)
show_mode = lower(char(string(show_mode)));
switch show_mode
    case 'all'
        tf = true;
    case 'first'
        tf = (i_mode == 1);
    case 'last'
        tf = (i_mode == n_mode);
    case 'none'
        tf = false;
    otherwise
        error('show_roi_labels must be ''all'', ''first'', ''last'', or ''none''.');
end
end


function lim = local_x_limits(x, normalization_name, roi_value_mode)
x = double(x(:));
x = x(isfinite(x));
if isempty(x)
    lim = [0 1];
    return;
end

if strcmpi(normalization_name, 'range')
    lim = [0 1];
    return;
end

if strcmpi(normalization_name, 'maxabs') || any(strcmpi(roi_value_mode, {'real_mean', 'imag_mean'}))
    mx = max(abs(x));
    if ~isfinite(mx) || mx <= 0
        lim = [-1 1];
    else
        lim = [-1 1] * 1.05 * mx;
    end
    return;
end

mx = max(x);
if ~isfinite(mx) || mx <= 0
    lim = [0 1];
else
    lim = [0 1.05 * mx];
end
end


function label = local_value_label(params)
base = '';
switch lower(char(string(params.roi_value_mode)))
    case 'abs_mean'
        base = '|mean(mode weight)| within ROI';
    case 'mean_abs'
        base = 'mean(|mode weight|) within ROI';
    case 'real_mean'
        base = 'mean(real(mode weight)) within ROI';
    case 'imag_mean'
        base = 'mean(imag(mode weight)) within ROI';
end

switch lower(char(string(params.mode_normalization)))
    case 'none'
        label = base;
    case 'range'
        label = sprintf('%s | normalized per mode (range)', base);
    case 'maxabs'
        label = sprintf('%s | normalized per mode (maxabs)', base);
end
end


function title_text = local_mode_title_text(params, selection_label, raw_idx, eval_i, i_mode)
if isfield(params, 'mode_title_texts') && ~isempty(params.mode_title_texts) && ...
        numel(params.mode_title_texts) >= i_mode && ...
        ~isempty(params.mode_title_texts{i_mode})
    title_text = char(string(params.mode_title_texts{i_mode}));
    return;
end

title_text = sprintf('%s -> raw %d | \\lambda = %.4g%+.4gi | |\\lambda| = %.4g', ...
    selection_label, raw_idx, real(eval_i), imag(eval_i), abs(eval_i));
end


function title_text = local_figure_title_text(params, plot_ctx)
override = char(string(local_get_field(params, 'figure_title_text', '')));
if ~isempty(override)
    title_text = override;
    return;
end

title_text = sprintf('%s | ROI-average Koopman mode summary | %s | ROI order %s', ...
    local_compact_run_label(plot_ctx.run_info), local_describe_map_space(plot_ctx), ...
    ternary(params.flip_roi_order, 'flipped', 'native'));
end


function map_space = local_describe_map_space(plot_ctx)
kind = local_get_field(plot_ctx, 'source_space_kind', 'unknown');
reduce = local_get_field(plot_ctx, 'feature_reduce', 'mean');
used_backprojection = ~plot_ctx.direct_voxel_mode;

switch kind
    case 'voxel'
        if used_backprojection
            map_space = 'SVD back-projected to voxel observables';
        else
            map_space = 'direct voxel observables';
        end
    case 'region_mean'
        if used_backprojection
            map_space = 'SVD back-projected to region means, then expanded to voxels';
        else
            map_space = 'region means expanded to voxels';
        end
    case 'multi_feature_voxel'
        if used_backprojection
            map_space = sprintf('SVD back-projected to multi-feature voxels (%s by voxel)', reduce);
        else
            map_space = sprintf('multi-feature voxel observables (%s by voxel)', reduce);
        end
    otherwise
        if used_backprojection
            map_space = sprintf('SVD back-projected to %s source space', kind);
        else
            map_space = sprintf('%s source space', kind);
        end
end
end


function [n_rows, n_cols, show_mode] = local_resolve_layout(n_mode, layout_mode, show_mode)
layout_mode = lower(char(string(layout_mode)));
switch layout_mode
    case 'row'
        n_rows = 1;
        n_cols = n_mode;
        if strcmpi(show_mode, 'last')
            show_mode = 'first';
        end
    case 'column'
        n_rows = n_mode;
        n_cols = 1;
    otherwise
        error('layout_mode must be ''row'' or ''column''.');
end
end


function [fig_width, fig_height] = local_resolve_figure_size(n_mode, n_roi, layout_mode)
layout_mode = lower(char(string(layout_mode)));
switch layout_mode
    case 'row'
        fig_width = max(1400, 500 * n_mode);
        fig_height = max(760, 11 * n_roi + 140);
    case 'column'
        panel_height = max(320, 10 * n_roi);
        fig_width = 1320;
        fig_height = panel_height * n_mode + 100;
    otherwise
        error('layout_mode must be ''row'' or ''column''.');
end
end


function label = local_compact_run_label(run_info)
dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', '')));
observable_mode = char(string(local_get_field(run_info, 'observable_mode', '')));
residual_form = char(string(local_get_field(run_info, 'residual_form', '')));
run_name = char(string(local_get_field(run_info, 'run_name', '')));

parts = {};
if ~isempty(dataset_stem)
    parts{end + 1} = dataset_stem; %#ok<AGROW>
end
if ~isempty(observable_mode)
    parts{end + 1} = observable_mode; %#ok<AGROW>
end
if ~isempty(residual_form)
    parts{end + 1} = residual_form; %#ok<AGROW>
end

if isempty(parts)
    label = run_name;
else
    label = strjoin(parts, ' | ');
end
end


function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
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
