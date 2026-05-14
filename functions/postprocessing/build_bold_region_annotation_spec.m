function spec = build_bold_region_annotation_spec(plot_ctx, slice_list, params)
%BUILD_BOLD_REGION_ANNOTATION_SPEC Build ROI text-annotation metadata.

if nargin < 1 || ~isstruct(plot_ctx)
    error('plot_ctx is required.');
end
if nargin < 2 || isempty(slice_list)
    slice_list = 1:size(plot_ctx.ana, 3);
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

slice_list = unique(round(double(slice_list(:).')), 'stable');
slice_list = slice_list(slice_list >= 1 & slice_list <= size(plot_ctx.ana, 3));

spec = struct();
spec.enabled = false;
spec.slice_annotations = cell(size(plot_ctx.ana, 3), 1);
spec.text_color = params.annotation_text_color;
spec.font_size = params.annotation_font_size;
spec.font_weight = params.annotation_font_weight;
spec.background_color = params.annotation_background_color;
spec.margin = params.annotation_margin;
spec.show_marker = params.annotation_show_marker;
spec.info = struct( ...
    'enabled', false, ...
    'observable_mode', char(string(local_get_field(plot_ctx.run_info, 'observable_mode', ''))), ...
    'requested_exact_names', {cellstr(string(params.annotation_exact_names(:)).')}, ...
    'requested_contains_tokens', {cellstr(string(params.annotation_contains_tokens(:)).')}, ...
    'matched_region_labels', {cell(0, 1)}, ...
    'slice_list', slice_list(:).', ...
    'repeat_each_slice', logical(params.annotation_repeat_each_slice), ...
    'reason', '');

if ~params.annotate_regions
    spec.info.reason = 'annotate_regions_disabled';
    return;
end

if isempty(slice_list)
    spec.info.reason = 'empty_slice_list';
    return;
end

observable_mode = char(string(local_get_field(plot_ctx.run_info, 'observable_mode', '')));
allowed_modes = cellstr(string(params.annotation_observable_modes(:)).');
if ~isempty(allowed_modes) && ~any(strcmpi(observable_mode, allowed_modes))
    spec.info.reason = sprintf('observable_mode_%s_not_in_annotation_set', observable_mode);
    return;
end

exact_names = cellstr(string(params.annotation_exact_names(:)).');
contains_tokens = cellstr(string(params.annotation_contains_tokens(:)).');
if isempty(exact_names) && isempty(contains_tokens)
    spec.info.reason = 'no_requested_region_names';
    return;
end

region_labels = cellstr(string(plot_ctx.region_labels(:)));
match_region = false(numel(region_labels), 1);
for i_region = 1:numel(region_labels)
    match_region(i_region) = local_is_requested_region_label(region_labels{i_region}, ...
        exact_names, contains_tokens);
end
matched_region_labels = region_labels(match_region);
if isempty(matched_region_labels)
    spec.info.reason = 'no_requested_region_labels_present';
    return;
end

coords_all = round(double(plot_ctx.coords));
vol_size = size(plot_ctx.ana);
in_bounds = coords_all(:, 1) >= 1 & coords_all(:, 1) <= vol_size(1) & ...
    coords_all(:, 2) >= 1 & coords_all(:, 2) <= vol_size(2) & ...
    coords_all(:, 3) >= 1 & coords_all(:, 3) <= vol_size(3);

for i_region = find(match_region(:)).'
    mask = plot_ctx.voxel_region_idx(:) == i_region & in_bounds(:);
    coords = coords_all(mask, :);
    if isempty(coords)
        continue;
    end

    slice_counts = zeros(1, numel(slice_list));
    for i_slice = 1:numel(slice_list)
        slice_counts(i_slice) = sum(coords(:, 3) == slice_list(i_slice));
    end

    if params.annotation_repeat_each_slice
        chosen_slices = slice_list(slice_counts >= params.annotation_min_voxels_per_slice);
    else
        [best_count, best_idx] = max(slice_counts);
        if isempty(best_count) || best_count < params.annotation_min_voxels_per_slice
            chosen_slices = [];
        else
            chosen_slices = slice_list(best_idx);
        end
    end

    for z = chosen_slices(:).'
        coords_z = coords(coords(:, 3) == z, 1:2);
        if isempty(coords_z)
            continue;
        end
        ann = struct();
        ann.label = region_labels{i_region};
        ann.region_index = i_region;
        ann.slice = z;
        ann.x = median(coords_z(:, 1), 'omitnan');
        ann.y = median(coords_z(:, 2), 'omitnan');
        ann.x_text = ann.x + params.annotation_text_offset(1);
        ann.y_text = ann.y + params.annotation_text_offset(2);
        ann.n_voxels_in_slice = size(coords_z, 1);
        spec.slice_annotations{z} = local_append_annotation(spec.slice_annotations{z}, ann);
    end
end

for z = slice_list(:).'
    spec.slice_annotations{z} = local_relax_slice_annotations( ...
        spec.slice_annotations{z}, vol_size, params.annotation_min_text_spacing);
end

spec.enabled = any(cellfun(@(x) ~isempty(x), spec.slice_annotations));
spec.info.enabled = spec.enabled;
spec.info.matched_region_labels = matched_region_labels(:);
if spec.enabled
    spec.info.reason = 'ok';
else
    spec.info.reason = 'no_slice_annotations_resolved';
end
end


function params = local_apply_defaults(params)
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
end


function tf = local_is_requested_region_label(label, exact_names, contains_tokens)
label_key = local_normalize_match_key(label);
tf = any(strcmp(label_key, local_normalize_match_key(exact_names)));
if tf
    return;
end

contains_keys = local_normalize_match_key(contains_tokens);
for i_token = 1:numel(contains_keys)
    if contains(label_key, contains_keys{i_token})
        tf = true;
        return;
    end
end
end


function key = local_normalize_match_key(value)
if iscell(value)
    key = cell(size(value));
    for i = 1:numel(value)
        key{i} = local_normalize_match_key(value{i});
    end
    return;
end
key = lower(strtrim(char(string(value))));
key = regexprep(key, '[^a-z0-9]+', '');
end


function out = local_append_annotation(in, ann)
if isempty(in)
    out = ann;
else
    out = [in; ann]; %#ok<AGROW>
end
end


function ann = local_relax_slice_annotations(ann, vol_size, min_spacing)
if isempty(ann) || numel(ann) < 2
    return;
end

[~, order] = sortrows([[ann.x_text].', [ann.y_text].'], [1 2]);
for i_order = 2:numel(order)
    prev_idx = order(i_order - 1);
    curr_idx = order(i_order);
    dx = abs(ann(curr_idx).x_text - ann(prev_idx).x_text);
    dy = abs(ann(curr_idx).y_text - ann(prev_idx).y_text);
    if dx <= min_spacing && dy < min_spacing
        ann(curr_idx).y_text = ann(prev_idx).y_text + min_spacing;
    end
end

y_max = vol_size(2);
x_max = vol_size(1);
for i = 1:numel(ann)
    ann(i).x_text = min(max(ann(i).x_text, 1), x_max);
    ann(i).y_text = min(max(ann(i).y_text, 1), y_max);
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
