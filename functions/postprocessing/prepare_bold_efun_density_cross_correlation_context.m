function ctx = prepare_bold_efun_density_cross_correlation_context( ...
        bold_post_input, density_sources, params)
%PREPARE_BOLD_EFUN_DENSITY_CROSS_CORRELATION_CONTEXT Load and align pipeline 8 inputs.

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

ctx = struct();
ctx.B = B;
ctx.params = params;
ctx.session = session;
ctx.dt = dt;
ctx.T = T;
ctx.max_lag_bins = max_lag_bins;
ctx.border_pad_bins = border_pad_bins;
ctx.lag_bins = lag_bins;
ctx.lag_sec = lag_sec;
ctx.session_id_by_sample = session_id_by_sample;
ctx.base_mask = base_mask;
ctx.feature_specs = feature_specs;
ctx.bold_feature_mats = bold_feature_mats;
ctx.density_sources = density_sources;
ctx.positive_lag_definition = 'density leads BOLD';
end


function B = local_load_bold_post(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load_mat_file_with_short_path(file);
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
            if isempty(X)
                continue;
            end
            X = abs(X);
            label = '|deconvolved BOLD eigenfunction|';
        case 'deconv_real'
            X = local_deconv_matrix(E);
            if isempty(X)
                continue;
            end
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
event_file = fullfile( ...
    io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 2, 'event_density'), ...
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


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
