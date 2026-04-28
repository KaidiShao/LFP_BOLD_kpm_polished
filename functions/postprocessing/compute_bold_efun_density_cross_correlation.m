function out = compute_bold_efun_density_cross_correlation(bold_post_input, density_sources, params)
%COMPUTE_BOLD_EFUN_DENSITY_CROSS_CORRELATION
% Session-aware lagged correlation between BOLD eigenfunctions and densities.
%
% Positive lag means density leads BOLD:
%   corr_lag(k) = corr(density(t), bold_feature(t + k))
%
% The implementation uses explicit masks. For each lag, a pair is valid
% only when both samples are finite, both are away from session borders, and
% both samples come from the same session.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input must be a BOLD_POST struct or a saved BOLD post MAT file.');
end
if nargin < 2
    density_sources = struct([]);
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

B = local_load_bold_post(bold_post_input);
[session, dt, T] = local_resolve_bold_session(B);
params.dt_current = dt;
[feature_specs, bold_feature_mats] = local_build_bold_features(B, params);

max_lag_bins = round(params.max_lag_sec / dt);
border_pad_bins = round(params.border_mask_sec / dt);
lag_bins = (-max_lag_bins:max_lag_bins).';
lag_sec = lag_bins * dt;
session_id_by_sample = local_session_id_vector(T, session);
base_mask = local_border_mask(T, session, border_pad_bins);

if isempty(density_sources)
    density_sources = local_default_density_sources(B);
end
density_sources = local_normalize_density_sources(density_sources);

all_rows = {};
source_results = {};

for i_src = 1:numel(density_sources)
    source = density_sources{i_src};
    candidates = local_load_density_candidates(source, params);
    if isempty(candidates)
        warning('No density candidates loaded for source %s.', source.name);
        continue;
    end

    best_abs = -Inf;
    best_source = [];
    best_candidate_idx = NaN;
    best_feature_idx = NaN;
    best_maps = [];

    for i_cand = 1:numel(candidates)
        cand = local_align_density_candidate(local_get_cell_or_array_item(candidates, i_cand), B, T, dt);

        for i_feat = 1:numel(feature_specs)
            maps = local_lagged_corr_maps(cand.X, bold_feature_mats{i_feat}, ...
                lag_bins, session_id_by_sample, base_mask, params);
            current_best = local_nanmax(maps.peak_abs_corr(:));
            if current_best > best_abs
                best_abs = current_best;
                best_source = cand;
                best_candidate_idx = i_cand;
                best_feature_idx = i_feat;
                best_maps = maps;
            end
        end
    end

    if isempty(best_source)
        continue;
    end

    source_result = struct();
    source_result.source = source;
    source_result.best_candidate = best_source;
    source_result.best_candidate_idx = best_candidate_idx;
    source_result.best_feature = feature_specs{best_feature_idx};
    source_result.best_bold_matrix = bold_feature_mats{best_feature_idx};
    source_result.best_maps = best_maps;
    source_result.lag_bins = lag_bins;
    source_result.lag_sec = lag_sec;
    source_result.positive_lag_definition = 'density leads BOLD';
    source_results{end + 1, 1} = source_result; %#ok<AGROW>

    rows_i = local_build_peak_rows(source_result, feature_specs{best_feature_idx});
    all_rows = [all_rows; rows_i]; %#ok<AGROW>
end

peak_table = local_rows_to_table(all_rows);
if ~isempty(peak_table)
    finite_peak = isfinite(peak_table.peak_abs_corr);
    peak_table = [ ...
        sortrows(peak_table(finite_peak, :), 'peak_abs_corr', 'descend'); ...
        peak_table(~finite_peak, :)];
end
if isempty(peak_table)
    top_table = peak_table;
else
    top_table = peak_table(1:min(height(peak_table), params.top_n), :);
end

out = struct();
out.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
out.bold_post_file = local_get_field(B, 'source_file', '');
out.params = params;
out.dt = dt;
out.T = T;
out.session = session;
out.max_lag_bins = max_lag_bins;
out.border_pad_bins = border_pad_bins;
out.border_mask = base_mask;
out.lag_bins = lag_bins;
out.lag_sec = lag_sec;
out.positive_lag_definition = 'density leads BOLD';
out.feature_specs = feature_specs;
out.density_sources = density_sources;
out.source_results = source_results;
out.peak_table = peak_table;
out.top_table = top_table;

if params.save_results
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end
    save_file = fullfile(params.save_dir, [params.save_tag, '.mat']);
    csv_file = fullfile(params.save_dir, [params.save_tag, '_peaks.csv']);
    top_csv_file = fullfile(params.save_dir, [params.save_tag, '_top.csv']);
    save(save_file, 'out', '-v7.3');
    if ~isempty(peak_table), writetable(peak_table, csv_file); end
    if ~isempty(top_table), writetable(top_table, top_csv_file); end
    out.save_paths = struct('mat_file', save_file, ...
        'peak_csv', csv_file, 'top_csv', top_csv_file);
end

if params.make_figures
    out.figure_paths = plot_bold_efun_density_cross_correlation_summary(out, params.plot);
end
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'max_lag_sec', 120);
params = local_set_default(params, 'border_mask_sec', params.max_lag_sec);
params = local_set_default(params, 'min_valid_samples', 20);
params = local_set_default(params, 'top_n', 5);
params = local_set_default(params, 'feature_names', ...
    {'efun_abs', 'efun_real', 'deconv_abs', 'deconv_real'});
params = local_set_default(params, 'density_field_preference', ...
    {'smoothed_density_mean', 'density_mean', ...
    'density_time_by_mode', 'density_time_by_component', 'density'});
params = local_set_default(params, 'density_transform', 'zscore');
params = local_set_default(params, 'bold_transform', 'zscore');
params = local_set_default(params, 'save_results', true);
params = local_set_default(params, 'save_dir', pwd);
params = local_set_default(params, 'save_tag', 'bold_efun_density_xcorr');
params = local_set_default(params, 'make_figures', true);
if ~isfield(params, 'plot') || isempty(params.plot)
    params.plot = struct();
end
params.plot = local_set_default(params.plot, 'save_dir', params.save_dir);
params.plot = local_set_default(params.plot, 'save_tag', params.save_tag);
params.plot = local_set_default(params.plot, 'top_n', params.top_n);
params.plot = local_set_default(params.plot, 'resolution', 220);
end


function B = local_load_bold_post(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    available = whos('-file', file);
    available_names = {available.name};
    load_names = intersect({'BOLD_POST', 'EDMD_outputs'}, available_names, 'stable');
    S = load(file, load_names{:});
    if isfield(S, 'BOLD_POST')
        B = S.BOLD_POST;
    elseif isfield(S, 'EDMD_outputs')
        B = struct('EDMD_outputs', S.EDMD_outputs);
    else
        error('File %s must contain BOLD_POST or EDMD_outputs.', file);
    end
    B.source_file = file;
    B = local_slim_bold_post(B);
elseif isstruct(input)
    B = input;
    if ~isfield(B, 'source_file')
        B.source_file = '';
    end
    B = local_slim_bold_post(B);
else
    error('bold_post_input must be a path or struct.');
end
if ~isfield(B, 'EDMD_outputs')
    error('BOLD post input does not contain EDMD_outputs.');
end
end


function B = local_slim_bold_post(B)
keep = {'source_file', 'run_info', 'session', 'dt', 'time_vec', 'EDMD_outputs'};
out = struct();
for i = 1:numel(keep)
    name = keep{i};
    if isfield(B, name)
        out.(name) = B.(name);
    end
end
B = out;
end


function [session, dt, T] = local_resolve_bold_session(B)
E = B.EDMD_outputs;
T = size(E.efuns, 1);
if isfield(B, 'session') && ~isempty(B.session)
    session = B.session;
else
    session = struct();
end
fields = {'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
for i = 1:numel(fields)
    name = fields{i};
    if (~isfield(session, name) || isempty(session.(name))) && ...
            isfield(E, name) && ~isempty(E.(name))
        session.(name) = E.(name);
    end
end
if ~isfield(session, 'session_start_idx') || isempty(session.session_start_idx)
    session.session_start_idx = 1;
    session.session_end_idx = T;
    session.session_lengths = T;
    session.session_ids = 1;
    session.border_idx = [];
end
session.session_start_idx = double(session.session_start_idx(:));
session.session_end_idx = double(session.session_end_idx(:));
if ~isfield(session, 'session_ids') || isempty(session.session_ids)
    session.session_ids = (1:numel(session.session_start_idx)).';
end
if ~isfield(session, 'border_idx') || isempty(session.border_idx)
    session.border_idx = session.session_end_idx(1:end-1);
end

if isfield(B, 'dt') && ~isempty(B.dt)
    dt = double(B.dt);
elseif isfield(E, 'dt') && ~isempty(E.dt)
    dt = double(E.dt);
elseif isfield(E, 'dx') && ~isempty(E.dx)
    dt = double(E.dx);
else
    dt = 2;
end
dt = dt(1);
end


function [specs, mats] = local_build_bold_features(B, params)
E = B.EDMD_outputs;
names = cellstr(string(params.feature_names(:)).');
specs = {};
mats = {};
for i = 1:numel(names)
    name = names{i};
    switch lower(name)
        case 'efun_abs'
            X = abs(E.efuns);
            label = '|BOLD eigenfunction|';
        case 'efun_real'
            X = real(E.efuns);
            label = 'Re(BOLD eigenfunction)';
        case 'deconv_abs'
            X = local_deconv_matrix(E);
            if isempty(X), continue; end
            X = abs(X);
            label = '|deconvolved BOLD eigenfunction|';
        case 'deconv_real'
            X = local_deconv_matrix(E);
            if isempty(X), continue; end
            X = real(X);
            label = 'Re(deconvolved BOLD eigenfunction)';
        otherwise
            error('Unsupported BOLD feature name: %s', name);
    end
    X = local_transform_matrix(X, params.bold_transform);
    spec = struct('name', name, 'label', label, 'n_modes', size(X, 2));
    specs{end + 1, 1} = spec; %#ok<AGROW>
    mats{end + 1, 1} = X; %#ok<AGROW>
end
end


function U = local_deconv_matrix(E)
U = [];
if isfield(E, 'deconv_efuns') && isstruct(E.deconv_efuns)
    if isfield(E.deconv_efuns, 'u_sel') && ~isempty(E.deconv_efuns.u_sel)
        U = E.deconv_efuns.u_sel;
    elseif isfield(E.deconv_efuns, 'u_all') && ~isempty(E.deconv_efuns.u_all)
        U = E.deconv_efuns.u_all;
    end
end
end


function density_sources = local_default_density_sources(B)
density_sources = {};
if ~isfield(B, 'run_info')
    return;
end
dataset_stem = local_get_field(B.run_info, 'dataset_stem', '');
if isempty(dataset_stem)
    return;
end
processed_root = io_project.get_project_processed_root();
event_file = fullfile(processed_root, dataset_stem, 'event_density', ...
    sprintf('%s_event_density_2s.mat', dataset_stem));
if exist(event_file, 'file') == 2
    source = struct();
    source.name = 'blp_event_density';
    source.type = 'event_density';
    source.file = event_file;
    density_sources = {source};
end
end


function sources = local_normalize_density_sources(sources)
if isempty(sources)
    return;
end
if ~iscell(sources)
    sources = num2cell(sources(:));
end
for i = 1:numel(sources)
    source = sources{i};
    if ~isfield(source, 'name') || isempty(source.name)
        source.name = sprintf('density%d', i);
    end
    if ~isfield(source, 'type') || isempty(source.type)
        source.type = 'generic_density';
    end
    sources{i} = source;
end
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
if ndims(X) > 2
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


function ids = local_session_id_vector(T, session)
ids = nan(T, 1);
starts = double(session.session_start_idx(:));
ends = double(session.session_end_idx(:));
session_ids = double(session.session_ids(:));
for i = 1:numel(starts)
    a = max(1, starts(i));
    b = min(T, ends(i));
    if a <= b
        ids(a:b) = session_ids(i);
    end
end
end


function mask = local_border_mask(T, session, pad_bins)
mask = true(T, 1);
starts = double(session.session_start_idx(:));
ends = double(session.session_end_idx(:));
for i = 1:numel(starts)
    a = max(1, starts(i));
    b = min(T, ends(i));
    if a > b
        continue;
    end
    if pad_bins > 0
        mask(a:min(b, a + pad_bins - 1)) = false;
        mask(max(a, b - pad_bins + 1):b) = false;
    end
end
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
fields = {'threshold_by_mode', 'threshold_by_component', 'params', 'summary'};
for i = 1:numel(fields)
    if isfield(D, fields{i})
        meta.(fields{i}) = D.(fields{i});
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
