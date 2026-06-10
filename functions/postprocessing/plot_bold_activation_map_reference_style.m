function [fig, info] = plot_bold_activation_map_reference_style(plot_ctx, raw_index, params)
%PLOT_BOLD_ACTIVATION_MAP_REFERENCE_STYLE Plot one reference-style BOLD mode activation map.

if nargin < 1 || ~isstruct(plot_ctx)
    error('plot_ctx is required.');
end
if nargin < 2 || isempty(raw_index)
    error('raw_index is required.');
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);
highlight_spec = local_normalize_highlight_spec(params);

kpm_modes = plot_ctx.kpm_modes;
if raw_index < 1 || raw_index > size(kpm_modes, 1)
    error('raw_index=%d is outside [1, %d].', raw_index, size(kpm_modes, 1));
end

lambda = plot_ctx.evalues(raw_index);
mode_obs = double(kpm_modes(raw_index, :));

if ~plot_ctx.direct_voxel_mode
    source_values = mode_obs * plot_ctx.coeff_t;
    map_space = local_describe_map_space(plot_ctx, true);
else
    source_values = mode_obs;
    map_space = local_describe_map_space(plot_ctx, false);
end

[voxel_values, map_info] = map_bold_source_values_to_voxels(source_values(:), plot_ctx, ...
    struct('feature_reduce', params.feature_reduce));

vals = local_prepare_values(voxel_values(:), params.value_mode);
[cmap, clim_use, value_label] = local_colormap_and_limits(params.value_mode);
act_vol = local_values_to_volume(vals, plot_ctx.coords, size(plot_ctx.ana));

slice_list = params.slice_list;
if isempty(slice_list)
    slice_list = unique(plot_ctx.coords(:, 3))';
    slice_list = slice_list(slice_list >= 1 & slice_list <= size(plot_ctx.ana, 3));
end
if isempty(slice_list)
    error('No valid slices were resolved for activation plotting.');
end
region_outline_spec = build_bold_region_highlight_spec(plot_ctx, params);
annotation_spec = build_bold_region_annotation_spec(plot_ctx, slice_list, params);

n_slice = numel(slice_list);
n_col = params.tiles_per_row;
n_row = ceil(n_slice / n_col);
fig_width = 1450;
if region_outline_spec.enabled && region_outline_spec.show_legend
    fig_width = fig_width + 40;
end
fig = figure('Color', 'w', 'Name', 'Reference-style BOLD activation map');
fig.Position = [80 80 fig_width 170 * n_row + 170];
tl = tiledlayout(n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_tiles = gobjects(n_slice, 1);

for jj = 1:n_slice
    z = slice_list(jj);
    nexttile;
    ax_tiles(jj) = gca;
    if z < 1 || z > size(plot_ctx.ana, 3)
        axis off;
        title(sprintf('S%d', z));
        continue;
    end
    img = squeeze(plot_ctx.ana(:, :, z))';
    overlay = squeeze(act_vol(:, :, z))';
    alpha_data = params.overlay_alpha * double(~isnan(overlay));

    image(local_gray_rgb(img, plot_ctx.bg_limits));
    axis image xy off;
    set(gca, 'Color', 'k');
    hold on;
    h = imagesc(overlay);
    set(h, 'AlphaData', alpha_data);
    colormap(gca, cmap);
    clim(clim_use);
    if highlight_spec.enabled
        local_draw_highlight_outline(highlight_spec, z);
    end
    if region_outline_spec.enabled
        local_draw_region_outlines(region_outline_spec, z);
    end
    if annotation_spec.enabled
        local_draw_region_annotations(annotation_spec, z);
    end
    title(sprintf('S%d', z), 'FontSize', 9, 'FontWeight', 'bold');
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = value_label;
if region_outline_spec.enabled && region_outline_spec.show_legend
    local_add_region_outline_legend(ax_tiles, tl, region_outline_spec);
end

if isfield(params, 'title_text') && ~isempty(params.title_text)
    title_text = char(string(params.title_text));
else
    title_text = sprintf('%s | raw %d | |lambda|=%.4g | %s | ROI(s): %s', ...
        local_get_field(plot_ctx.run_info, 'run_name', 'BOLD'), ...
        raw_index, abs(lambda), map_space, local_region_summary(plot_ctx.region_labels));
end
if highlight_spec.enabled
    title_text = sprintf('%s | core ROI outlined', title_text);
end
if region_outline_spec.enabled
    title_text = sprintf('%s | ROI boundaries outlined', title_text);
end
if annotation_spec.enabled
    title_text = sprintf('%s | ROI labels annotated', title_text);
end
sgtitle(title_text, 'Interpreter', 'none');

info = struct();
info.raw_index = raw_index;
info.lambda = lambda;
info.value_mode = params.value_mode;
info.map_space = map_space;
info.slice_list = slice_list;
info.region_labels = plot_ctx.region_labels;
info.roi_ts_file = local_get_field(plot_ctx, 'roi_ts_file', '');
info.mapping_info = map_info;
info.core_roi_highlight = highlight_spec.info;
info.region_outline = region_outline_spec.info;
info.region_annotations = annotation_spec.info;
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'value_mode', 'abs');
params = local_set_default(params, 'slice_list', 1:20);
params = local_set_default(params, 'tiles_per_row', 10);
params = local_set_default(params, 'overlay_alpha', 0.86);
params = local_set_default(params, 'feature_reduce', 'mean');
params = local_set_default(params, 'title_text', '');
params = local_set_default(params, 'highlight_spec', struct());
params = local_set_default(params, 'annotate_regions', false);
params = local_set_default(params, 'annotation_observable_modes', {});
params = local_set_default(params, 'annotation_exact_names', {});
params = local_set_default(params, 'annotation_contains_tokens', {});
params = local_set_default(params, 'annotation_repeat_each_slice', false);
params = local_set_default(params, 'annotation_min_voxels_per_slice', 5);
params = local_set_default(params, 'annotation_text_color', [1 1 1]);
params = local_set_default(params, 'annotation_font_size', 8);
params = local_set_default(params, 'annotation_font_weight', 'bold');
params = local_set_default(params, 'annotation_background_color', [0 0 0]);
params = local_set_default(params, 'annotation_margin', 1);
params = local_set_default(params, 'annotation_text_offset', [0 -2]);
params = local_set_default(params, 'annotation_min_text_spacing', 5);
params = local_set_default(params, 'annotation_show_marker', false);
params = local_set_default(params, 'outline_regions', false);
params = local_set_default(params, 'outline_observable_modes', {});
params = local_set_default(params, 'outline_exact_names', {});
params = local_set_default(params, 'outline_contains_tokens', {});
params = local_set_default(params, 'outline_use_annotation_names', true);
params = local_set_default(params, 'outline_line_width', 1.6);
params = local_set_default(params, 'outline_colors', []);
params = local_set_default(params, 'outline_show_legend', true);
end


function vals = local_prepare_values(raw_vals, value_mode)
switch lower(value_mode)
    case 'abs'
        vals = local_minmax(abs(double(raw_vals)));
    case 'real'
        vals = real(double(raw_vals));
        max_abs = max(abs(vals), [], 'omitnan');
        if isfinite(max_abs) && max_abs > 0
            vals = vals ./ max_abs;
        end
    otherwise
        error('value_mode must be ''abs'' or ''real''.');
end
end


function vol = local_values_to_volume(vals, coords, vol_size)
coords = round(double(coords));
in_bounds = coords(:, 1) >= 1 & coords(:, 1) <= vol_size(1) & ...
    coords(:, 2) >= 1 & coords(:, 2) <= vol_size(2) & ...
    coords(:, 3) >= 1 & coords(:, 3) <= vol_size(3);
vol = nan(vol_size);
idx = sub2ind(vol_size, coords(in_bounds, 1), coords(in_bounds, 2), coords(in_bounds, 3));
vol(idx) = vals(in_bounds);
end


function [cmap, clim_use, value_label] = local_colormap_and_limits(value_mode)
switch lower(value_mode)
    case 'abs'
        cmap = turbo(256);
        clim_use = [0 1];
        value_label = 'normalized |mode weight|';
    case 'real'
        cmap = local_redblue(256);
        clim_use = [-1 1];
        value_label = 'normalized real(mode weight)';
end
end


function rgb = local_gray_rgb(img, lim)
x = (double(img) - lim(1)) ./ max(eps, lim(2) - lim(1));
x = min(max(x, 0), 1);
rgb = repmat(x, 1, 1, 3);
end


function y = local_minmax(x)
mn = min(x, [], 'omitnan');
mx = max(x, [], 'omitnan');
if ~isfinite(mn) || ~isfinite(mx) || mx <= mn
    y = zeros(size(x));
else
    y = (x - mn) ./ (mx - mn);
end
end


function cmap = local_redblue(n)
x = linspace(-1, 1, n)';
cmap = zeros(n, 3);
cmap(:, 1) = max(0, x);
cmap(:, 3) = max(0, -x);
cmap(:, 2) = 1 - abs(x);
cmap = min(max(0.15 + 0.85 * cmap, 0), 1);
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


function summary = local_region_summary(region_labels)
labels = cellstr(string(region_labels(:)));
n = numel(labels);
if n == 0
    summary = 'none';
elseif n <= 4
    summary = strjoin(labels(:).', ', ');
else
    preview = strjoin(labels(1:4).', ', ');
    summary = sprintf('%s, ... (%d total)', preview, n);
end
end


function map_space = local_describe_map_space(plot_ctx, used_backprojection)
kind = local_get_field(plot_ctx, 'source_space_kind', 'unknown');
reduce = local_get_field(plot_ctx, 'feature_reduce', 'mean');

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


function highlight_spec = local_normalize_highlight_spec(params)
highlight_spec = struct();
highlight_spec.enabled = false;
highlight_spec.mask_volume = [];
highlight_spec.line_color = [1 1 1];
highlight_spec.line_width = 1.6;
highlight_spec.info = struct('enabled', false, 'matched_region_labels', {cell(0, 1)});

if ~isfield(params, 'highlight_spec') || ~isstruct(params.highlight_spec) || ...
        isempty(params.highlight_spec)
    return;
end

highlight_spec = params.highlight_spec;
if ~isfield(highlight_spec, 'enabled') || isempty(highlight_spec.enabled)
    highlight_spec.enabled = false;
end
if ~isfield(highlight_spec, 'mask_volume')
    highlight_spec.mask_volume = [];
end
if ~isfield(highlight_spec, 'line_color') || isempty(highlight_spec.line_color)
    highlight_spec.line_color = [1 1 1];
end
if ~isfield(highlight_spec, 'line_width') || isempty(highlight_spec.line_width)
    highlight_spec.line_width = 1.6;
end
if ~isfield(highlight_spec, 'info') || ~isstruct(highlight_spec.info)
    highlight_spec.info = struct('enabled', logical(highlight_spec.enabled), ...
        'matched_region_labels', {cell(0, 1)});
end
end


function local_draw_highlight_outline(highlight_spec, z)
if isempty(highlight_spec.mask_volume) || z < 1 || z > size(highlight_spec.mask_volume, 3)
    return;
end
mask2d = squeeze(highlight_spec.mask_volume(:, :, z)).';
if ~any(mask2d(:))
    return;
end
contour(double(mask2d), [0.5 0.5], 'Color', highlight_spec.line_color, ...
    'LineWidth', highlight_spec.line_width, 'LineStyle', '-');
end


function local_draw_region_outlines(region_outline_spec, z)
if ~isfield(region_outline_spec, 'region_specs') || isempty(region_outline_spec.region_specs)
    return;
end

for i_region = 1:numel(region_outline_spec.region_specs)
    one_spec = region_outline_spec.region_specs(i_region);
    if isempty(one_spec.mask_volume) || z < 1 || z > size(one_spec.mask_volume, 3)
        continue;
    end
    mask2d = squeeze(one_spec.mask_volume(:, :, z)).';
    if ~any(mask2d(:))
        continue;
    end
    contour(double(mask2d), [0.5 0.5], 'Color', one_spec.line_color, ...
        'LineWidth', one_spec.line_width, 'LineStyle', '-');
end
end


function local_add_region_outline_legend(ax_tiles, tl, region_outline_spec)
if isempty(ax_tiles) || ~any(isgraphics(ax_tiles))
    return;
end
valid_ax = ax_tiles(find(isgraphics(ax_tiles), 1, 'first'));
if isempty(valid_ax)
    return;
end

n_region = numel(region_outline_spec.region_specs);
legend_handles = gobjects(n_region, 1);
legend_labels = cell(n_region, 1);
for i_region = 1:n_region
    one_spec = region_outline_spec.region_specs(i_region);
    legend_handles(i_region) = plot(valid_ax, nan, nan, '-', ...
        'Color', one_spec.line_color, ...
        'LineWidth', max(2.2, one_spec.line_width), ...
        'DisplayName', one_spec.label);
    legend_labels{i_region} = one_spec.label;
end

lgd = legend(valid_ax, legend_handles, legend_labels, ...
    'Interpreter', 'none', ...
    'Orientation', 'horizontal', ...
    'NumColumns', min(max(1, n_region), 4), ...
    'Box', 'off');
lgd.Layout.Tile = 'south';
lgd.FontSize = 9;
lgd.Title.String = 'ROI boundary colors';
lgd.Title.FontWeight = 'bold';
end


function local_draw_region_annotations(annotation_spec, z)
if isempty(annotation_spec.slice_annotations) || z < 1 || z > numel(annotation_spec.slice_annotations)
    return;
end
ann = annotation_spec.slice_annotations{z};
if isempty(ann)
    return;
end

for i_ann = 1:numel(ann)
    if annotation_spec.show_marker
        plot(ann(i_ann).x, ann(i_ann).y, '.', ...
            'Color', annotation_spec.text_color, 'MarkerSize', 10);
    end
    text(ann(i_ann).x_text, ann(i_ann).y_text, ann(i_ann).label, ...
        'Color', annotation_spec.text_color, ...
        'FontSize', annotation_spec.font_size, ...
        'FontWeight', annotation_spec.font_weight, ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'Interpreter', 'none', ...
        'BackgroundColor', annotation_spec.background_color, ...
        'Margin', annotation_spec.margin);
end
end
