function spec = build_bold_region_highlight_spec(plot_ctx, params)
%BUILD_BOLD_REGION_HIGHLIGHT_SPEC Build arbitrary ROI outline metadata.

if nargin < 1 || ~isstruct(plot_ctx)
    error('plot_ctx is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);
exact_names = local_resolve_outline_exact_names(params);
contains_tokens = local_resolve_outline_contains_tokens(params);

spec = struct();
spec.enabled = false;
spec.region_specs = repmat(local_empty_region_spec(), 0, 1);
spec.legend_labels = cell(0, 1);
spec.legend_colors = zeros(0, 3);
spec.show_legend = params.outline_show_legend;
spec.info = struct( ...
    'enabled', false, ...
    'observable_mode', char(string(local_get_field(plot_ctx.run_info, 'observable_mode', ''))), ...
    'requested_exact_names', {cellstr(string(exact_names(:)).')}, ...
    'requested_contains_tokens', {cellstr(string(contains_tokens(:)).')}, ...
    'matched_region_labels', {cell(0, 1)}, ...
    'legend_labels', {cell(0, 1)}, ...
    'reason', '');

if ~params.outline_regions
    spec.info.reason = 'outline_regions_disabled';
    return;
end

observable_mode = char(string(local_get_field(plot_ctx.run_info, 'observable_mode', '')));
allowed_modes = cellstr(string(params.outline_observable_modes(:)).');
if ~isempty(allowed_modes) && ~any(strcmpi(observable_mode, allowed_modes))
    spec.info.reason = sprintf('observable_mode_%s_not_in_outline_set', observable_mode);
    return;
end

exact_names = cellstr(string(exact_names(:)).');
contains_tokens = cellstr(string(contains_tokens(:)).');
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

region_indices = find(match_region);
vol_size = size(plot_ctx.ana);
colors = local_resolve_outline_colors(numel(region_indices), params);

for i_region = 1:numel(region_indices)
    region_idx = region_indices(i_region);
    coords = round(double(plot_ctx.coords(plot_ctx.voxel_region_idx(:) == region_idx, :)));
    in_bounds = coords(:, 1) >= 1 & coords(:, 1) <= vol_size(1) & ...
        coords(:, 2) >= 1 & coords(:, 2) <= vol_size(2) & ...
        coords(:, 3) >= 1 & coords(:, 3) <= vol_size(3);
    coords = coords(in_bounds, :);
    if isempty(coords)
        continue;
    end

    mask_volume = false(vol_size);
    idx = sub2ind(vol_size, coords(:, 1), coords(:, 2), coords(:, 3));
    mask_volume(idx) = true;

    one_spec = local_empty_region_spec();
    one_spec.label = region_labels{region_idx};
    one_spec.region_index = region_idx;
    one_spec.mask_volume = mask_volume;
    one_spec.line_color = colors(i_region, :);
    one_spec.line_width = params.outline_line_width;
    spec.region_specs(end + 1, 1) = one_spec; %#ok<AGROW>
end

if isempty(spec.region_specs)
    spec.info.reason = 'matched_region_voxels_out_of_bounds';
    return;
end

spec.legend_labels = {spec.region_specs.label}.';
spec.legend_colors = vertcat(spec.region_specs.line_color);
spec.enabled = true;
spec.info.enabled = true;
spec.info.matched_region_labels = matched_region_labels(:);
spec.info.legend_labels = spec.legend_labels;
spec.info.reason = 'ok';
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'outline_regions', false);
params = local_set_default(params, 'outline_observable_modes', {});
params = local_set_default(params, 'outline_exact_names', {});
params = local_set_default(params, 'outline_contains_tokens', {});
params = local_set_default(params, 'outline_use_annotation_names', true);
params = local_set_default(params, 'outline_line_width', 1.6);
params = local_set_default(params, 'outline_colors', []);
params = local_set_default(params, 'outline_show_legend', true);
params = local_set_default(params, 'annotation_exact_names', {});
params = local_set_default(params, 'annotation_contains_tokens', {});
end


function names = local_resolve_outline_exact_names(params)
names = params.outline_exact_names;
if isempty(names) && params.outline_use_annotation_names && isfield(params, 'annotation_exact_names')
    names = params.annotation_exact_names;
end
end


function tokens = local_resolve_outline_contains_tokens(params)
tokens = params.outline_contains_tokens;
if isempty(tokens) && params.outline_use_annotation_names && isfield(params, 'annotation_contains_tokens')
    tokens = params.annotation_contains_tokens;
end
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


function colors = local_resolve_outline_colors(n_region, params)
if n_region <= 0
    colors = zeros(0, 3);
    return;
end

if ~isempty(params.outline_colors)
    colors_in = double(params.outline_colors);
    if size(colors_in, 2) == 3 && size(colors_in, 1) >= n_region
        colors = colors_in(1:n_region, :);
        return;
    end
end

base = [ ...
    0.98 0.98 0.98; ...
    1.00 0.40 0.25; ...
    0.20 0.90 0.35; ...
    1.00 0.85 0.20; ...
    0.25 0.80 1.00; ...
    1.00 0.35 0.85; ...
    0.95 0.55 0.15; ...
    0.65 0.55 1.00];

if n_region <= size(base, 1)
    colors = base(1:n_region, :);
    return;
end

extra = lines(n_region);
extra = 0.2 + 0.8 * extra;
colors = extra;
colors(1:size(base, 1), :) = base;
end


function spec = local_empty_region_spec()
spec = struct( ...
    'label', '', ...
    'region_index', NaN, ...
    'mask_volume', [], ...
    'line_color', [1 1 1], ...
    'line_width', 1.6);
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
