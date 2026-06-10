function stats = compute_bold_efun_density_cross_correlation_source_results(ctx)
%COMPUTE_BOLD_EFUN_DENSITY_CROSS_CORRELATION_SOURCE_RESULTS Compute pipeline 8 xcorr statistics.

all_rows = {};
source_results = {};

for i_src = 1:numel(ctx.density_sources)
    source = ctx.density_sources{i_src};
    candidates = local_load_density_candidates(source, ctx.params);
    if isempty(candidates)
        warning('No density candidates loaded for source %s.', source.name);
        continue;
    end

    aligned_candidates = cell(numel(candidates), 1);
    for i_cand = 1:numel(candidates)
        aligned_candidates{i_cand} = local_align_density_candidate( ...
            local_get_cell_or_array_item(candidates, i_cand), ...
            ctx.B, ctx.T, ctx.dt);
    end

    for i_feat = 1:numel(ctx.feature_specs)
        best_abs = -Inf;
        best_source = [];
        best_candidate_idx = NaN;
        best_maps = [];

        for i_cand = 1:numel(aligned_candidates)
            cand = aligned_candidates{i_cand};
            maps = local_lagged_corr_maps( ...
                cand.X, ctx.bold_feature_mats{i_feat}, ...
                ctx.lag_bins, ctx.session_id_by_sample, ctx.base_mask, ctx.params);
            current_best = local_nanmax(maps.peak_abs_corr(:));
            if current_best > best_abs
                best_abs = current_best;
                best_source = cand;
                best_candidate_idx = i_cand;
                best_maps = maps;
            end
        end

        if isempty(best_source)
            continue;
        end

        source_result = struct();
        source_result.source = source;
        source_result.best_candidate = best_source;
        source_result.best_candidate_idx = best_candidate_idx;
        source_result.best_feature = ctx.feature_specs{i_feat};
        source_result.best_bold_matrix = ctx.bold_feature_mats{i_feat};
        source_result.best_maps = best_maps;
        source_result.lag_bins = ctx.lag_bins;
        source_result.lag_sec = ctx.lag_sec;
        source_result.positive_lag_definition = ctx.positive_lag_definition;
        source_results{end + 1, 1} = source_result; %#ok<AGROW>

        rows_i = local_build_peak_rows(source_result, ctx.feature_specs{i_feat});
        all_rows = [all_rows; rows_i]; %#ok<AGROW>
    end
end

peak_table = local_rows_to_table(all_rows);
peak_table = local_sort_peak_table(peak_table);
if isempty(peak_table)
    top_table = peak_table;
else
    top_table = peak_table(1:min(height(peak_table), ctx.params.top_n), :);
end
[peak_table_by_density, top_table_by_density, density_group_names, ...
    density_group_fields] = local_build_density_group_tables( ...
    peak_table, ctx.params.top_n);
[peak_table_by_density_feature, top_table_by_density_feature, ...
    density_feature_group_names, density_feature_group_fields, ...
    density_feature_density_names, density_feature_feature_names] = ...
    local_build_density_feature_group_tables(peak_table, ctx.params.top_n);

stats = struct();
stats.source_results = source_results;
stats.peak_table = peak_table;
stats.top_table = top_table;
stats.peak_table_by_density = peak_table_by_density;
stats.top_table_by_density = top_table_by_density;
stats.density_group_names = density_group_names;
stats.density_group_fields = density_group_fields;
stats.peak_table_by_density_feature = peak_table_by_density_feature;
stats.top_table_by_density_feature = top_table_by_density_feature;
stats.density_feature_group_names = density_feature_group_names;
stats.density_feature_group_fields = density_feature_group_fields;
stats.density_feature_density_names = density_feature_density_names;
stats.density_feature_feature_names = density_feature_feature_names;
end


function candidates = local_load_density_candidates(source, params)
files = {};
if isfield(source, 'files') && ~isempty(source.files)
    files = local_as_cellstr(source.files);
elseif isfield(source, 'file') && ~isempty(source.file)
    files = local_as_cellstr(source.file);
elseif isfield(source, 'manifest_file') && ~isempty(source.manifest_file)
    files = local_density_files_from_manifest(source.manifest_file);
elseif isfield(source, 'file_pattern') && ~isempty(source.file_pattern)
    L = dir(source.file_pattern);
    files = arrayfun(@(d) fullfile(d.folder, d.name), L, 'UniformOutput', false);
end

candidates = {};
for i = 1:numel(files)
    file = files{i};
    if exist(file, 'file') ~= 2
        continue;
    end
    S = load(file);
    cand = local_extract_density_candidate(S, file, source, params);
    if isempty(cand)
        continue;
    end
    candidates{end + 1, 1} = cand; %#ok<AGROW>
end

if isempty(candidates) && isfield(source, 'D') && isstruct(source.D)
    cand = local_extract_density_candidate(source, '', source, params);
    if ~isempty(cand)
        candidates = {cand};
    end
end
end


function out = local_as_cellstr(value)
if isempty(value)
    out = {};
elseif ischar(value)
    out = {value};
elseif isstring(value)
    out = cellstr(value(:));
elseif iscell(value)
    out = cellstr(string(value(:)));
else
    error('Expected a char, string, or cell array of paths.');
end
out = out(:).';
end


function files = local_density_files_from_manifest(manifest_file)
files = {};
S = load(manifest_file);
if isfield(S, 'manifest') && isfield(S.manifest, 'table')
    T = S.manifest.table;
elseif isfield(S, 'T')
    T = S.T;
else
    return;
end

names = {'mat_file', 'density_mat_file', 'result_mat_file'};
for i = 1:numel(names)
    if ismember(names{i}, T.Properties.VariableNames)
        files = cellstr(string(T.(names{i})));
        files = files(~cellfun(@isempty, files));
        return;
    end
end
end


function cand = local_extract_density_candidate(S, file, source, params)
cand = [];
D = [];
if isfield(S, 'E')
    D = S.E;
elseif isfield(S, 'D')
    D = S.D;
elseif isfield(S, 'out')
    D = S.out;
elseif isfield(S, 'density')
    D = S;
end
if isempty(D) || ~isstruct(D)
    return;
end

[X, field_used] = local_pick_density_matrix(D, params.density_field_preference);
if isempty(X)
    return;
end

X = double(X);
if ~ismatrix(X)
    X = reshape(X, size(X, 1), []);
end
X = local_transform_matrix(X, params.density_transform);

t = [];
if isfield(D, 't_centers') && ~isempty(D.t_centers)
    t = double(D.t_centers(:));
elseif isfield(D, 'time') && ~isempty(D.time)
    t = double(D.time(:));
end

labels = local_density_labels(D, size(X, 2));
cand = struct();
cand.name = source.name;
cand.type = source.type;
cand.file = file;
cand.field_used = field_used;
cand.X = X;
cand.t = t;
cand.labels = labels;
cand.meta = local_density_meta(D);
end


function [X, field_used] = local_pick_density_matrix(D, preferences)
X = [];
field_used = '';
for i = 1:numel(preferences)
    name = preferences{i};
    if isfield(D, name) && ~isempty(D.(name))
        X = D.(name);
        field_used = name;
        return;
    end
end
if isfield(D, 'bands') && ~isempty(D.bands) && isfield(D.bands(1), 'smoothed_density_mean')
    X = cat(2, D.bands.smoothed_density_mean);
    field_used = 'bands.smoothed_density_mean';
end
end


function cand = local_align_density_candidate(cand, B, T, dt)
if size(cand.X, 1) == T
    return;
end
if isempty(cand.t)
    warning('Density %s has %d rows but BOLD has %d; truncating to common length.', ...
        cand.name, size(cand.X, 1), T);
    n = min(size(cand.X, 1), T);
    X = nan(T, size(cand.X, 2));
    X(1:n, :) = cand.X(1:n, :);
    cand.X = X;
    return;
end

if isfield(B, 'time_vec') && numel(B.time_vec) == T
    t_bold = double(B.time_vec(:));
else
    t_bold = (0:T-1)' * dt;
end
cand.X = interp1(cand.t, cand.X, t_bold, 'linear', NaN);
end


function maps = local_lagged_corr_maps(X_density, Y_bold, lag_bins, session_id, base_mask, params)
n_density = size(X_density, 2);
n_bold = size(Y_bold, 2);
n_lag = numel(lag_bins);
corr_cube = nan(n_density, n_bold, n_lag);
valid_count = zeros(n_density, n_bold, n_lag, 'uint32');
n = size(X_density, 1);

for i_lag = 1:n_lag
    lag = lag_bins(i_lag);
    if lag > 0
        ix = (1:(n - lag)).';
        iy = ((1 + lag):n).';
    elseif lag < 0
        a = -lag;
        ix = ((1 + a):n).';
        iy = (1:(n - a)).';
    else
        ix = (1:n).';
        iy = ix;
    end

    valid_pair = base_mask(ix) & base_mask(iy) & ...
        (session_id(ix) == session_id(iy));
    X = X_density(ix, :);
    Y = Y_bold(iy, :);
    X(~valid_pair, :) = NaN;
    Y(~valid_pair, :) = NaN;
    [C, N] = local_pairwise_corr_matrix(X, Y, params.min_valid_samples);
    corr_cube(:, :, i_lag) = C;
    valid_count(:, :, i_lag) = uint32(N);
end

[peak_corr, peak_abs_corr, peak_lag_bins, peak_lag_sec, zero_corr] = ...
    local_peak_maps(corr_cube, lag_bins, params);
maps = struct();
maps.corr_cube = corr_cube;
maps.valid_count = valid_count;
maps.peak_corr = peak_corr;
maps.peak_abs_corr = peak_abs_corr;
maps.peak_lag_bins = peak_lag_bins;
maps.peak_lag_sec = peak_lag_sec;
maps.zero_corr = zero_corr;
end


function [C, N] = local_pairwise_corr_matrix(X, Y, min_n)
X = double(X);
Y = double(Y);
nx = size(X, 2);
ny = size(Y, 2);
C = nan(nx, ny);
N = zeros(nx, ny);
for i = 1:nx
    xi = X(:, i);
    for j = 1:ny
        yj = Y(:, j);
        valid = isfinite(xi) & isfinite(yj);
        n = nnz(valid);
        N(i, j) = n;
        if n < min_n
            continue;
        end
        xv = xi(valid);
        yv = yj(valid);
        sx = std(xv);
        sy = std(yv);
        if sx == 0 || sy == 0
            continue;
        end
        xv = xv - mean(xv);
        yv = yv - mean(yv);
        C(i, j) = (xv.' * yv) ./ ((n - 1) * sx * sy);
    end
end
end


function [peak_corr, peak_abs_corr, peak_lag_bins, peak_lag_sec, zero_corr] = ...
        local_peak_maps(corr_cube, lag_bins, params)
[nd, nb, ~] = size(corr_cube);
peak_corr = nan(nd, nb);
peak_abs_corr = nan(nd, nb);
peak_lag_bins = nan(nd, nb);
peak_lag_sec = nan(nd, nb);
zero_idx = find(lag_bins == 0, 1, 'first');
if isempty(zero_idx)
    zero_corr = nan(nd, nb);
else
    zero_corr = corr_cube(:, :, zero_idx);
end
for i = 1:nd
    for j = 1:nb
        curve = squeeze(corr_cube(i, j, :));
        if all(~isfinite(curve))
            continue;
        end
        [~, idx] = max(abs(curve));
        peak_corr(i, j) = curve(idx);
        peak_abs_corr(i, j) = abs(curve(idx));
        peak_lag_bins(i, j) = lag_bins(idx);
        peak_lag_sec(i, j) = lag_bins(idx) * params.dt_current;
    end
end
end


function rows = local_build_peak_rows(source_result, feature_spec)
maps = source_result.best_maps;
cand = source_result.best_candidate;
[n_density, n_bold] = size(maps.peak_abs_corr);
rows = cell(n_density * n_bold, 13);
k = 0;
for i_den = 1:n_density
    for i_mode = 1:n_bold
        k = k + 1;
        rows(k, :) = { ...
            cand.name, cand.type, cand.file, cand.field_used, ...
            feature_spec.name, i_den, cand.labels{i_den}, i_mode, ...
            maps.zero_corr(i_den, i_mode), maps.peak_corr(i_den, i_mode), ...
            maps.peak_abs_corr(i_den, i_mode), ...
            maps.peak_lag_bins(i_den, i_mode), ...
            maps.peak_lag_sec(i_den, i_mode)};
    end
end
rows = rows(1:k, :);
end


function T = local_sort_peak_table(T)
if isempty(T)
    return;
end
finite_peak = isfinite(T.peak_abs_corr);
T = [ ...
    sortrows(T(finite_peak, :), 'peak_abs_corr', 'descend'); ...
    T(~finite_peak, :)];
end


function [peak_by_density, top_by_density, density_names, field_names] = ...
        local_build_density_group_tables(peak_table, top_n)
peak_by_density = struct();
top_by_density = struct();
density_names = {};
field_names = {};
if isempty(peak_table) || ~ismember('density_name', peak_table.Properties.VariableNames)
    return;
end

all_names = cellstr(string(peak_table.density_name));
[density_names, ~] = unique(all_names, 'stable');
if isempty(density_names)
    density_names = {};
    return;
end
density_names = density_names(:).';
field_names = cell(size(density_names));
for i_group = 1:numel(density_names)
    density_name = density_names{i_group};
    field_name = local_density_group_field_name(density_name);
    mask = strcmp(all_names, density_name);
    peak_table_i = local_sort_peak_table(peak_table(mask, :));
    if isempty(peak_table_i)
        top_table_i = peak_table_i;
    else
        top_table_i = peak_table_i(1:min(height(peak_table_i), top_n), :);
    end
    peak_by_density.(field_name) = peak_table_i;
    top_by_density.(field_name) = top_table_i;
    field_names{i_group} = field_name;
end
end


function [peak_by_group, top_by_group, group_names, group_fields, ...
        density_names, feature_names] = ...
        local_build_density_feature_group_tables(peak_table, top_n)
peak_by_group = struct();
top_by_group = struct();
group_names = {};
group_fields = {};
density_names = {};
feature_names = {};
if isempty(peak_table) || ...
        ~ismember('density_name', peak_table.Properties.VariableNames) || ...
        ~ismember('bold_feature', peak_table.Properties.VariableNames)
    return;
end

all_density = cellstr(string(peak_table.density_name));
all_feature = cellstr(string(peak_table.bold_feature));
group_keys = strcat(all_density, {'|'}, all_feature);
[unique_keys, first_idx] = unique(group_keys, 'stable');
if isempty(unique_keys)
    return;
end

n_groups = numel(unique_keys);
group_names = cell(1, n_groups);
group_fields = cell(1, n_groups);
density_names = cell(1, n_groups);
feature_names = cell(1, n_groups);
for i_group = 1:n_groups
    density_name = all_density{first_idx(i_group)};
    feature_name = all_feature{first_idx(i_group)};
    field_name = local_density_feature_group_field_name(density_name, feature_name);
    mask = strcmp(all_density, density_name) & strcmp(all_feature, feature_name);
    peak_table_i = local_sort_peak_table(peak_table(mask, :));
    if isempty(peak_table_i)
        top_table_i = peak_table_i;
    else
        top_table_i = peak_table_i(1:min(height(peak_table_i), top_n), :);
    end
    peak_by_group.(field_name) = peak_table_i;
    top_by_group.(field_name) = top_table_i;
    group_names{i_group} = sprintf('%s | %s', density_name, feature_name);
    group_fields{i_group} = field_name;
    density_names{i_group} = density_name;
    feature_names{i_group} = feature_name;
end
end


function field_name = local_density_group_field_name(name)
field_name = char(matlab.lang.makeValidName(char(string(name))));
if isempty(field_name)
    field_name = 'density_group';
end
end


function field_name = local_density_feature_group_field_name(density_name, feature_name)
field_name = char(matlab.lang.makeValidName(sprintf('%s__%s', ...
    char(string(density_name)), char(string(feature_name)))));
if isempty(field_name)
    field_name = 'density_feature_group';
end
end


function T = local_rows_to_table(rows)
if isempty(rows)
    T = table();
    return;
end
T = cell2table(rows, 'VariableNames', { ...
    'density_name', 'density_type', 'density_file', 'density_field', ...
    'bold_feature', 'density_index', 'density_label', 'bold_mode_index', ...
    'zero_corr', 'peak_corr', 'peak_abs_corr', ...
    'peak_lag_bins', 'peak_lag_sec'});
end


function X = local_transform_matrix(X, mode)
switch lower(char(mode))
    case 'none'
        return;
    case 'zscore'
        mu = mean(X, 1, 'omitnan');
        sig = std(X, 0, 1, 'omitnan');
        sig(sig == 0 | ~isfinite(sig)) = NaN;
        X = (X - mu) ./ sig;
    otherwise
        error('Unknown transform mode: %s', mode);
end
end


function labels = local_density_labels(D, K)
if isfield(D, 'band_labels') && numel(D.band_labels) >= K
    labels = cellstr(string(D.band_labels(1:K)));
elseif isfield(D, 'mode_index') && numel(D.mode_index) >= K
    labels = cellstr(compose('mode%d', double(D.mode_index(1:K))));
elseif isfield(D, 'component_index') && numel(D.component_index) >= K
    labels = cellstr(compose('comp%d', double(D.component_index(1:K))));
else
    labels = cellstr(compose('density%d', 1:K));
end
labels = labels(:);
end


function meta = local_density_meta(D)
meta = struct();
fields = {'threshold_by_mode', 'threshold_by_component', 'summary'};
for i = 1:numel(fields)
    if isfield(D, fields{i})
        meta.(fields{i}) = D.(fields{i});
    end
end
if isfield(D, 'params')
    meta.params = local_slim_metadata(D.params);
end
end


function value = local_slim_metadata(value)
max_elements = 10000;
if isnumeric(value) || islogical(value)
    if numel(value) > max_elements
        value = struct( ...
            'omitted_large_array', true, ...
            'class_name', class(value), ...
            'size', size(value));
    end
elseif isstring(value) || ischar(value)
    if numel(value) > max_elements
        value = struct( ...
            'omitted_large_text', true, ...
            'class_name', class(value), ...
            'size', size(value));
    end
elseif iscell(value)
    if numel(value) > max_elements
        value = struct( ...
            'omitted_large_cell', true, ...
            'class_name', class(value), ...
            'size', size(value));
    else
        for i = 1:numel(value)
            value{i} = local_slim_metadata(value{i});
        end
    end
elseif isstruct(value)
    for i_item = 1:numel(value)
        names = fieldnames(value(i_item));
        for i_name = 1:numel(names)
            name = names{i_name};
            value(i_item).(name) = local_slim_metadata(value(i_item).(name));
        end
    end
end
end


function v = local_nanmax(x)
x = x(isfinite(x));
if isempty(x)
    v = NaN;
else
    v = max(x);
end
end


function item = local_get_cell_or_array_item(value, idx)
if iscell(value)
    item = value{idx};
else
    item = value(idx);
end
end
