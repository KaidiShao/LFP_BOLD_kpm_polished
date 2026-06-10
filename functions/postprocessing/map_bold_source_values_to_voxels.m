function [voxel_values, info] = map_bold_source_values_to_voxels(source_values, plot_ctx, params)
%MAP_BOLD_SOURCE_VALUES_TO_VOXELS Map observable-space weights into voxel space.

if nargin < 1 || isempty(source_values)
    error('source_values is required.');
end
if nargin < 2 || ~isstruct(plot_ctx)
    error('plot_ctx is required.');
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params, plot_ctx);

source_values = source_values(:);
if ~isfield(plot_ctx, 'source_to_voxel_indices') || isempty(plot_ctx.source_to_voxel_indices)
    error('plot_ctx.source_to_voxel_indices is required.');
end
if numel(source_values) ~= numel(plot_ctx.source_to_voxel_indices)
    error(['Observable-space mode has %d values but activation mapping expects %d ' ...
        'source variables.'], numel(source_values), numel(plot_ctx.source_to_voxel_indices));
end

n_voxel = size(plot_ctx.coords, 1);
accum = complex(zeros(n_voxel, 1));
counts = zeros(n_voxel, 1);
unresolved = {};

for i_source = 1:numel(source_values)
    target_idx = plot_ctx.source_to_voxel_indices{i_source};
    if isempty(target_idx)
        if isfield(plot_ctx, 'source_variable_labels') && ...
                numel(plot_ctx.source_variable_labels) >= i_source
            unresolved{end + 1, 1} = plot_ctx.source_variable_labels{i_source}; %#ok<AGROW>
        else
            unresolved{end + 1, 1} = sprintf('source_%04d', i_source); %#ok<AGROW>
        end
        continue;
    end
    accum(target_idx) = accum(target_idx) + source_values(i_source);
    counts(target_idx) = counts(target_idx) + 1;
end

if ~isempty(unresolved)
    preview = unresolved(1:min(6, numel(unresolved)));
    error(['Activation mapping still has %d unresolved source variable(s). ' ...
        'First labels: %s'], numel(unresolved), strjoin(preview(:).', ', '));
end

[accum, counts, duplicate_fill_info] = local_fill_uncovered_duplicate_coords( ...
    accum, counts, plot_ctx);
covered = counts > 0;

voxel_values = complex(nan(n_voxel, 1));
switch lower(params.feature_reduce)
    case 'mean'
        voxel_values(covered) = accum(covered) ./ counts(covered);
    case 'sum'
        voxel_values(covered) = accum(covered);
    otherwise
        error('feature_reduce must be ''mean'' or ''sum''.');
end

if any(~covered) && ~params.allow_uncovered_voxels
    error('Activation mapping left %d voxel(s) uncovered.', sum(~covered));
elseif any(~covered) && params.warn_uncovered_voxels
    warning('Activation mapping left %d voxel(s) uncovered; leaving them as NaN.', sum(~covered));
end

info = struct();
info.feature_reduce = params.feature_reduce;
info.n_voxels = n_voxel;
info.n_sources = numel(source_values);
info.n_covered_voxels = sum(covered);
info.n_uncovered_voxels = sum(~covered);
info.coverage_fraction = sum(covered) ./ max(1, n_voxel);
info.n_initial_covered_voxels = duplicate_fill_info.n_initial_covered_voxels;
info.n_initial_uncovered_voxels = duplicate_fill_info.n_initial_uncovered_voxels;
info.n_filled_by_duplicate_coords = duplicate_fill_info.n_filled_by_duplicate_coords;
info.min_features_per_voxel = min(counts);
info.max_features_per_voxel = max(counts);
info.mean_features_per_voxel = mean(counts);
end


function [accum, counts, info] = local_fill_uncovered_duplicate_coords(accum, counts, plot_ctx)
covered_initial = counts > 0;
info = struct();
info.n_initial_covered_voxels = sum(covered_initial);
info.n_initial_uncovered_voxels = sum(~covered_initial);
info.n_filled_by_duplicate_coords = 0;

if all(covered_initial) || ~isfield(plot_ctx, 'coords') || isempty(plot_ctx.coords)
    return;
end

coords = round(double(plot_ctx.coords));
if size(coords, 1) ~= numel(counts) || size(coords, 2) < 3
    return;
end

keys = string(coords(:, 1)) + "_" + string(coords(:, 2)) + "_" + string(coords(:, 3));
[~, ~, group_idx] = unique(keys, 'stable');
[sorted_group_idx, order] = sort(group_idx);
group_starts = [1; find(diff(sorted_group_idx) ~= 0) + 1];
group_ends = [group_starts(2:end) - 1; numel(sorted_group_idx)];

for i_group = 1:numel(group_starts)
    idx = order(group_starts(i_group):group_ends(i_group));
    covered_idx = idx(counts(idx) > 0);
    missing_idx = idx(counts(idx) == 0);
    if isempty(covered_idx) || isempty(missing_idx)
        continue;
    end

    values = accum(covered_idx) ./ counts(covered_idx);
    valid = isfinite(real(values)) & isfinite(imag(values));
    if ~any(valid)
        continue;
    end
    fill_value = mean(values(valid));
    accum(missing_idx) = fill_value;
    counts(missing_idx) = 1;
    info.n_filled_by_duplicate_coords = info.n_filled_by_duplicate_coords + numel(missing_idx);
end
end


function params = local_apply_defaults(params, plot_ctx)
default_reduce = 'mean';
if nargin >= 2 && isstruct(plot_ctx) && isfield(plot_ctx, 'feature_reduce') && ...
        ~isempty(plot_ctx.feature_reduce)
    default_reduce = plot_ctx.feature_reduce;
end
if ~isfield(params, 'feature_reduce') || isempty(params.feature_reduce)
    params.feature_reduce = default_reduce;
end
params.feature_reduce = char(string(params.feature_reduce));
if ~isfield(params, 'allow_uncovered_voxels') || isempty(params.allow_uncovered_voxels)
    params.allow_uncovered_voxels = false;
end
if ~isfield(params, 'warn_uncovered_voxels') || isempty(params.warn_uncovered_voxels)
    params.warn_uncovered_voxels = false;
end
end
