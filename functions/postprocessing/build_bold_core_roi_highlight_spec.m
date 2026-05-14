function spec = build_bold_core_roi_highlight_spec(plot_ctx, params)
%BUILD_BOLD_CORE_ROI_HIGHLIGHT_SPEC Build core-ROI outline metadata for maps.

if nargin < 1 || ~isstruct(plot_ctx)
    error('plot_ctx is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

spec = struct();
spec.enabled = false;
spec.mask_volume = false(size(plot_ctx.ana));
spec.line_color = params.core_roi_outline_color;
spec.line_width = params.core_roi_line_width;
spec.info = struct( ...
    'enabled', false, ...
    'observable_mode', char(string(local_get_field(plot_ctx.run_info, 'observable_mode', ''))), ...
    'requested_exact_names', {cellstr(string(params.core_roi_exact_names(:)).')}, ...
    'requested_contains_tokens', {cellstr(string(params.core_roi_contains_tokens(:)).')}, ...
    'matched_region_labels', {cell(0, 1)}, ...
    'reason', '');

if ~params.highlight_core_rois
    spec.info.reason = 'highlight_core_rois_disabled';
    return;
end

observable_mode = char(string(local_get_field(plot_ctx.run_info, 'observable_mode', '')));
allowed_modes = cellstr(string(params.core_roi_observable_modes(:)).');
if ~isempty(allowed_modes) && ~any(strcmpi(observable_mode, allowed_modes))
    spec.info.reason = sprintf('observable_mode_%s_not_in_highlight_set', observable_mode);
    return;
end

region_labels = cellstr(string(plot_ctx.region_labels(:)));
match_region = false(numel(region_labels), 1);
for i_region = 1:numel(region_labels)
    match_region(i_region) = local_is_core_roi_label(region_labels{i_region}, params);
end
matched_region_labels = region_labels(match_region);
if isempty(matched_region_labels)
    spec.info.reason = 'no_core_roi_labels_present';
    return;
end

match_voxel = ismember(plot_ctx.voxel_region_idx(:), find(match_region));
coords = round(double(plot_ctx.coords(match_voxel, :)));
vol_size = size(plot_ctx.ana);
in_bounds = coords(:, 1) >= 1 & coords(:, 1) <= vol_size(1) & ...
    coords(:, 2) >= 1 & coords(:, 2) <= vol_size(2) & ...
    coords(:, 3) >= 1 & coords(:, 3) <= vol_size(3);
coords = coords(in_bounds, :);
if isempty(coords)
    spec.info.reason = 'matched_core_roi_voxels_out_of_bounds';
    return;
end

mask_volume = false(vol_size);
idx = sub2ind(vol_size, coords(:, 1), coords(:, 2), coords(:, 3));
mask_volume(idx) = true;

spec.enabled = true;
spec.mask_volume = mask_volume;
spec.info.enabled = true;
spec.info.matched_region_labels = matched_region_labels(:);
spec.info.reason = 'ok';
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'highlight_core_rois', true);
params = local_set_default(params, 'core_roi_observable_modes', ...
    {'svd', 'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100', ...
    'roi_mean', 'roi_mean_slow_band_power'});
params = local_set_default(params, 'core_roi_exact_names', ...
    {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'});
params = local_set_default(params, 'core_roi_contains_tokens', {'V1'});
params = local_set_default(params, 'core_roi_outline_color', [1 1 1]);
params = local_set_default(params, 'core_roi_line_width', 1.6);
end


function tf = local_is_core_roi_label(label, params)
label = char(string(label));

exact_names = cellstr(string(params.core_roi_exact_names(:)).');
tf = any(strcmpi(label, exact_names));
if tf
    return;
end

contains_tokens = cellstr(string(params.core_roi_contains_tokens(:)).');
for i_token = 1:numel(contains_tokens)
    if contains(label, contains_tokens{i_token}, 'IgnoreCase', true)
        tf = true;
        return;
    end
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
