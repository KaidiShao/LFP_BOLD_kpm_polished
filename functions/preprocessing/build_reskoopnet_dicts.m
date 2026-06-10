function dict = build_reskoopnet_dicts(D, cfg, output_root, params)
% Build a ResKoopNet dictionary from loaded BLP data and saved spectrograms.
%
% Inputs
%   D           Output struct from load_blp_dataset
%   cfg         Dataset config struct
%   output_root Root folder for processed outputs
%   params      Struct with dictionary-building parameters
%
% Output
%   dict        Struct with saved file paths and metadata
%
% Notes
%   - BLP channels are kept as raw observables.
%   - Spectrogram observables are built from saved region-mean spectrograms.
%   - 0-low_full_max_hz: keep every frequency bin.
%   - (low_full_max_hz, high_max_hz]: average every high_group_size bins.
%   - For complex spectrograms, the recommended mode is 'complex_split',
%     which stores real and imaginary parts as separate real observables.

%% =========================
%  Defaults
%  =========================
if nargin < 3 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 4
    params = struct();
end

if ~isfield(params, 'spec_mode')
    params.spec_mode = 'complex_split';   % 'abs' or 'complex_split'
end
if ~isfield(params, 'low_full_max_hz')
    params.low_full_max_hz = 50;
end
if ~isfield(params, 'high_max_hz')
    params.high_max_hz = 250;
end
if ~isfield(params, 'high_group_size')
    params.high_group_size = 2;
end
if ~isfield(params, 'chunk_size')
    params.chunk_size = 200000;
end
if ~isfield(params, 'precision')
    params.precision = 'single';          % 'single' or 'double'
end
if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end

if ~ismember(params.spec_mode, {'abs', 'complex_split'})
    error('params.spec_mode must be ''abs'' or ''complex_split''.');
end

%% =========================
%  Locate spectrogram files
%  =========================
spec_files = resolve_regionmean_spectrogram_files(cfg, output_root);
abs_file = local_require_existing_file(spec_files.abs_file, spec_files.abs_name);
complex_file = local_require_existing_file(spec_files.complex_file, spec_files.complex_name);

%% =========================
%  Build feature plan
%  =========================
plan = build_feature_plan(D, cfg, abs_file, complex_file, params);
obs_info = plan.obs_info;
freqs = plan.freqs;
regions = plan.regions;
n_obs = plan.n_obs;
if isfield(D, 'n_time') && ~isempty(D.n_time)
    n_time = D.n_time;
else
    n_time = size(D.data, 1);
end

fprintf('Total time points   : %d\n', n_time);
fprintf('Total observables   : %d\n', n_obs);
fprintf('BLP observables     : %d\n', numel(D.selected_channels));
fprintf('Spectrogram regions : %d\n', numel(plan.regions));

%% =========================
%  Prepare save paths
%  =========================
save_dir = io_project.get_pipeline_stage_dir(output_root, cfg, 1, 'dictionary');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_tag = sprintf('%s_low%g_high%g_g%d_%s_%s', ...
    cfg.file_stem, ...
    params.low_full_max_hz, ...
    params.high_max_hz, ...
    params.high_group_size, ...
    params.spec_mode, ...
    params.precision);

save_tag = strrep(save_tag, '.', 'p');

save_file = fullfile(save_dir, [save_tag, '.mat']);
info_csv_file = fullfile(save_dir, [save_tag, '_obs_info.csv']);
work_save_file = fullfile(save_dir, [save_tag, '.building.mat']);
work_info_csv_file = fullfile(save_dir, [save_tag, '_obs_info.building.csv']);
work_progress_file = fullfile(save_dir, [save_tag, '.building_progress.txt']);

if ~params.force_recompute && exist(save_file, 'file') == 2 && exist(info_csv_file, 'file') == 2
    try
        local_validate_dictionary_output(save_file, D, plan, false);
        fprintf('Dictionary cache hit:\n  %s\n', save_file);
        dict = struct();
        dict.save_file = save_file;
        dict.info_csv_file = info_csv_file;
        dict.n_time = n_time;
        dict.n_obs = n_obs;
        dict.dx = [];
        dict.dt = [];
        dict.fs = [];
        dict.obs_info = obs_info;
        return;
    catch ME
        warning('Existing dictionary failed validation and will be rebuilt: %s', ME.message);
    end
end

resume_building = false;
if exist(work_save_file, 'file') == 2
    if params.force_recompute
        delete(work_save_file);
    else
        try
            local_validate_dictionary_shape(work_save_file, D, plan, params);
            resume_building = true;
            fprintf('Dictionary partial output found; attempting resumable build:\n  %s\n', work_save_file);
        catch ME
            warning('Existing partial dictionary failed shape validation and will be rebuilt: %s', ME.message);
            delete(work_save_file);
        end
    end
end
if exist(work_info_csv_file, 'file') == 2
    delete(work_info_csv_file);
end
if ~resume_building && exist(work_progress_file, 'file') == 2
    delete(work_progress_file);
end

%% =========================
%  Prepare source spectrogram matfile
%  =========================
Mspec = matfile(plan.source_file);
spec_size = size(Mspec, plan.source_var_name);
if spec_size(2) ~= n_time
    error('Spectrogram time dimension (%d) does not match raw data length (%d).', spec_size(2), n_time);
end

%% =========================
%  Preallocate output matfile
%  =========================
Mout = matfile(work_save_file, 'Writable', true);
chunk_size = params.chunk_size;
n_chunks = ceil(n_time / chunk_size);

if resume_building
    first_chunk = local_resume_chunk_from_progress( ...
        work_progress_file, params, n_chunks, n_time, n_obs, work_save_file);
    if isempty(first_chunk)
        first_chunk = local_find_first_incomplete_chunk(Mout, D, plan, n_time, n_chunks, chunk_size);
    end
    if first_chunk > 1
        first_chunk = first_chunk - 1;
        fprintf('Rewriting chunk %d as a resume guard.\n', first_chunk);
    end
    fprintf('Resuming dictionary chunks at %d/%d.\n', min(first_chunk, n_chunks), n_chunks);
elseif strcmp(params.precision, 'single')
    Mout.obs(n_time, n_obs) = single(0);
    first_chunk = 1;
    local_write_dictionary_progress(work_progress_file, params, 0, n_chunks, 0, n_time, n_obs, work_save_file);
elseif strcmp(params.precision, 'double')
    Mout.obs(n_time, n_obs) = double(0);
    first_chunk = 1;
    local_write_dictionary_progress(work_progress_file, params, 0, n_chunks, 0, n_time, n_obs, work_save_file);
else
    error('Unsupported params.precision.');
end

%% =========================
%  Chunk loop
%  =========================
for k = first_chunk:n_chunks
    idx1 = (k - 1) * chunk_size + 1;
    idx2 = min(k * chunk_size, n_time);
    idx = idx1:idx2;

    fprintf('[%d/%d] Building dictionary chunk %d:%d ...\n', k, n_chunks, idx1, idx2);

    obs_chunk = build_dictionary_chunk(D, Mspec, idx, plan);

    if strcmp(params.precision, 'single')
        obs_chunk = single(obs_chunk);
    else
        obs_chunk = double(obs_chunk);
    end

    Mout.obs(idx, :) = obs_chunk;
    local_write_dictionary_progress(work_progress_file, params, k, n_chunks, idx2, n_time, n_obs, work_save_file);
end

%% =========================
%  Save metadata
%  =========================
session_ids = D.session_ids;
session_lengths = D.session_lengths;
session_dx = D.session_dx;
dx = io_utils.resolve_uniform_dx(session_dx);
dt = dx;
fs = [];
if ~isempty(dx)
    fs = 1 / dx;
end
sampling_period = dx;
sample_period = dx;
sampling_frequency = fs;
session_start_idx = D.session_start_idx;
session_end_idx = D.session_end_idx;
border_idx = D.border_idx;
selected_channels = D.selected_channels;
channel_sites = cfg.channels.sites(D.selected_channels);

save(work_save_file, ...
    'obs_info', 'params', 'freqs', 'regions', ...
    'dx', 'dt', 'fs', 'sampling_period', 'sample_period', 'sampling_frequency', ...
    'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'channel_sites', ...
    '-append');

writetable(obs_info, work_info_csv_file);
local_validate_dictionary_output(work_save_file, D, plan, true);

if exist(save_file, 'file') == 2
    delete(save_file);
end
if exist(info_csv_file, 'file') == 2
    delete(info_csv_file);
end
movefile(work_save_file, save_file);
movefile(work_info_csv_file, info_csv_file);
if exist(work_progress_file, 'file') == 2
    delete(work_progress_file);
end

%% =========================
%  Pack output
%  =========================
dict = struct();
dict.save_file = save_file;
dict.info_csv_file = info_csv_file;
dict.n_time = n_time;
dict.n_obs = n_obs;
dict.dx = dx;
dict.dt = dt;
dict.fs = fs;
dict.obs_info = obs_info;
end

function plan = build_feature_plan(D, cfg, abs_file, complex_file, params)
% Build one canonical feature-layout plan used by both data writing and metadata.

meta = load(abs_file, 'freqs', 'regions');
freqs = double(meta.freqs(:));
regions = meta.regions;

low_idx = find(freqs >= 0 & freqs <= params.low_full_max_hz);
high_idx = find(freqs > params.low_full_max_hz & freqs <= params.high_max_hz);
high_groups = make_consecutive_groups(high_idx, params.high_group_size);

if isempty(low_idx) && isempty(high_groups)
    error('No spectrogram frequency bins were selected.');
end

spec_features = build_spec_feature_plan(regions, freqs, low_idx, high_groups, params.spec_mode);

plan = struct();
plan.freqs = freqs;
plan.regions = regions;
plan.n_blp = numel(D.selected_channels);
plan.spec_features = spec_features;
plan.obs_info = build_observable_info(D, cfg, spec_features);
plan.n_obs = height(plan.obs_info);

if strcmp(params.spec_mode, 'abs')
    plan.source_file = abs_file;
    plan.source_var_name = 'tmpall_mean_abs';
else
    plan.source_file = complex_file;
    plan.source_var_name = 'tmpall_mean_complex';
end
end

function obs_chunk = build_dictionary_chunk(D, Mspec, idx, plan)
% Build one output chunk using a precomputed feature plan.

raw_chunk = double(io_raw.read_blp_data_slice(D, idx));   % time x n_channels
n_time_chunk = size(raw_chunk, 1);
obs_chunk = zeros(n_time_chunk, plan.n_obs);
obs_chunk(:, 1:plan.n_blp) = raw_chunk;

spec_chunk = double(Mspec.(plan.source_var_name)(:, idx, :));   % freq x time x region

for i_feat = 1:numel(plan.spec_features)
    feat = plan.spec_features(i_feat);
    spec_region = spec_chunk(:, :, feat.region_idx);

    switch feat.part
        case 'abs'
            spec_part = spec_region;
        case 'real'
            spec_part = real(spec_region);
        case 'imag'
            spec_part = imag(spec_region);
        otherwise
            error('Unsupported spectrogram part: %s', feat.part);
    end

    if feat.freq_idx_start == feat.freq_idx_end
        obs_chunk(:, plan.n_blp + i_feat) = transpose(spec_part(feat.freq_idx_start, :));
    else
        obs_chunk(:, plan.n_blp + i_feat) = transpose(mean( ...
            spec_part(feat.freq_idx_start:feat.freq_idx_end, :), 1));
    end
end
end

function groups = make_consecutive_groups(idx_vec, group_size)
% Split a vector of indices into consecutive groups.

groups = {};
if isempty(idx_vec)
    return;
end

n = numel(idx_vec);
g = 1;

for i = 1:group_size:n
    j = min(i + group_size - 1, n);
    groups{g, 1} = idx_vec(i:j);
    g = g + 1;
end
end

function spec_features = build_spec_feature_plan(regions, freqs, low_idx, high_groups, spec_mode)
% Expand one spectrogram layout plan in the exact saved column order.

spec_features = struct( ...
    'part', {}, ...
    'region_idx', {}, ...
    'region_label', {}, ...
    'freq_idx_start', {}, ...
    'freq_idx_end', {}, ...
    'freq_start_hz', {}, ...
    'freq_end_hz', {}, ...
    'aggregation', {}, ...
    'name', {});

if strcmp(spec_mode, 'abs')
    parts_to_add = {'abs'};
else
    parts_to_add = {'real', 'imag'};
end

row_id = 0;
for r = 1:numel(regions)
    for p = 1:numel(parts_to_add)
        part_name = parts_to_add{p};

        for i = 1:numel(low_idx)
            fi = low_idx(i);
            row_id = row_id + 1;
            spec_features(row_id).part = part_name; %#ok<AGROW>
            spec_features(row_id).region_idx = r;
            spec_features(row_id).region_label = regions{r};
            spec_features(row_id).freq_idx_start = fi;
            spec_features(row_id).freq_idx_end = fi;
            spec_features(row_id).freq_start_hz = freqs(fi);
            spec_features(row_id).freq_end_hz = freqs(fi);
            spec_features(row_id).aggregation = 'single_bin';
            spec_features(row_id).name = sprintf('spec_%s_%s_f%.3f', ...
                part_name, regions{r}, freqs(fi));
        end

        for g = 1:numel(high_groups)
            idxg = high_groups{g};
            fi1 = idxg(1);
            fi2 = idxg(end);
            row_id = row_id + 1;
            spec_features(row_id).part = part_name; %#ok<AGROW>
            spec_features(row_id).region_idx = r;
            spec_features(row_id).region_label = regions{r};
            spec_features(row_id).freq_idx_start = fi1;
            spec_features(row_id).freq_idx_end = fi2;
            spec_features(row_id).freq_start_hz = freqs(fi1);
            spec_features(row_id).freq_end_hz = freqs(fi2);
            spec_features(row_id).aggregation = sprintf('mean_%d_bins', numel(idxg));
            spec_features(row_id).name = sprintf('spec_%s_%s_f%.3f_%.3f', ...
                part_name, regions{r}, freqs(fi1), freqs(fi2));
        end
    end
end
end

function obs_info = build_observable_info(D, cfg, spec_features)
% Build a table describing every observable dimension.

sites = cfg.channels.sites(D.selected_channels);
n_blp = numel(D.selected_channels);
n_spec = numel(spec_features);
n_rows = n_blp + n_spec;

dim_id = (1:n_rows).';
source = cell(n_rows, 1);
name = cell(n_rows, 1);
channel_idx = nan(n_rows, 1);
channel_site = cell(n_rows, 1);
region_idx = nan(n_rows, 1);
region_label = cell(n_rows, 1);
part = cell(n_rows, 1);
freq_idx_start = nan(n_rows, 1);
freq_idx_end = nan(n_rows, 1);
freq_start_hz = nan(n_rows, 1);
freq_end_hz = nan(n_rows, 1);
aggregation = cell(n_rows, 1);

for c = 1:n_blp
    source{c} = 'blp';
    name{c} = sprintf('blp_ch%02d_%s', D.selected_channels(c), sites{c});
    channel_idx(c) = D.selected_channels(c);
    channel_site{c} = sites{c};
    region_label{c} = '';
    part{c} = 'raw';
    aggregation{c} = 'raw';
end

for i_feat = 1:n_spec
    row = n_blp + i_feat;
    feat = spec_features(i_feat);
    source{row} = 'spectrogram';
    name{row} = feat.name;
    channel_site{row} = '';
    region_idx(row) = feat.region_idx;
    region_label{row} = feat.region_label;
    part{row} = feat.part;
    freq_idx_start(row) = feat.freq_idx_start;
    freq_idx_end(row) = feat.freq_idx_end;
    freq_start_hz(row) = feat.freq_start_hz;
    freq_end_hz(row) = feat.freq_end_hz;
    aggregation{row} = feat.aggregation;
end

obs_info = table( ...
    dim_id, source, name, ...
    channel_idx, channel_site, ...
    region_idx, region_label, ...
    part, ...
    freq_idx_start, freq_idx_end, ...
    freq_start_hz, freq_end_hz, ...
    aggregation);
end

function file = local_require_existing_file(file, target_name)
if exist(file, 'file') ~= 2
    error('Spectrogram file not found: %s', target_name);
end
end

function local_validate_dictionary_shape(save_file, D, plan, params)
% Validate that an existing MAT-file can be safely resumed in-place.

info = whos('-file', save_file, 'obs');
if isempty(info)
    error('obs variable was not saved in %s.', save_file);
end

n_time = D.n_time;
if isempty(n_time)
    n_time = size(D.data, 1);
end

expected_size = [n_time, plan.n_obs];
if ~isequal(info.size, expected_size)
    error('obs size [%s] does not match expected [%s].', ...
        num2str(info.size), num2str(expected_size));
end

if ~strcmp(info.class, params.precision)
    error('obs class is %s, expected %s.', info.class, params.precision);
end
end

function first_chunk = local_resume_chunk_from_progress(progress_file, params, n_chunks, ...
    n_time, n_obs, work_save_file)
% Return next chunk from a matching progress sidecar, or [] if it is stale.

first_chunk = [];
if exist(progress_file, 'file') ~= 2
    return;
end

progress = local_read_dictionary_progress(progress_file);
required_fields = {'spec_mode', 'precision', 'chunk_size', 'completed_chunk', ...
    'n_chunks', 'n_time', 'n_obs', 'work_save_file'};
for i = 1:numel(required_fields)
    if ~isfield(progress, required_fields{i})
        return;
    end
end

if ~strcmp(progress.spec_mode, params.spec_mode) || ...
        ~strcmp(progress.precision, params.precision) || ...
        progress.chunk_size ~= params.chunk_size || ...
        progress.n_chunks ~= n_chunks || ...
        progress.n_time ~= n_time || ...
        progress.n_obs ~= n_obs || ...
        ~strcmp(progress.work_save_file, work_save_file)
    return;
end

completed_chunk = max(0, min(progress.completed_chunk, n_chunks));
first_chunk = completed_chunk + 1;
fprintf('Dictionary progress sidecar found: completed chunk %d/%d.\n', ...
    completed_chunk, n_chunks);
end

function progress = local_read_dictionary_progress(progress_file)
progress = struct();
txt = fileread(progress_file);
lines = regexp(txt, '\r\n|\n|\r', 'split');

numeric_fields = {'chunk_size', 'completed_chunk', 'n_chunks', ...
    'completed_idx2', 'n_time', 'n_obs'};

for i = 1:numel(lines)
    line = strtrim(lines{i});
    if isempty(line)
        continue;
    end
    parts = regexp(line, '^([^=]+)=(.*)$', 'tokens', 'once');
    if isempty(parts)
        continue;
    end
    key = matlab.lang.makeValidName(strtrim(parts{1}));
    value = strtrim(parts{2});
    if ismember(key, numeric_fields)
        value_num = str2double(value);
        if ~isnan(value_num)
            progress.(key) = value_num;
        end
    else
        progress.(key) = value;
    end
end
end

function first_chunk = local_find_first_incomplete_chunk(Mout, D, plan, n_time, n_chunks, chunk_size)
% Infer the first chunk that still needs writing.
%
% Older interrupted builds did not have a progress sidecar. Probe each chunk
% sparsely and keep only the contiguous prefix whose raw and spectrogram
% columns are nonzero. The caller rewrites one preceding chunk as a guard.

first_chunk = n_chunks + 1;

for k = 1:n_chunks
    idx1 = (k - 1) * chunk_size + 1;
    idx2 = min(k * chunk_size, n_time);
    probe_idx = local_probe_indices(idx1, idx2);

    try
        raw_probe = local_read_matfile_probe(Mout, probe_idx, 1:plan.n_blp);
        spec_probe = [];
        if plan.n_obs > plan.n_blp
            spec_cols = local_probe_spec_columns(plan);
            spec_probe = local_read_matfile_probe(Mout, probe_idx, spec_cols);
        end
    catch
        first_chunk = k;
        return;
    end

    raw_written = any(raw_probe(:) ~= 0);
    spec_written = true;
    if plan.n_obs > plan.n_blp
        spec_written = any(spec_probe(:) ~= 0);
    end

    if ~(raw_written && spec_written)
        first_chunk = k;
        return;
    end
end
end

function probe_idx = local_probe_indices(idx1, idx2)
n_probe = min(5, idx2 - idx1 + 1);
probe_idx = unique(round(linspace(idx1, idx2, n_probe)));
end

function spec_cols = local_probe_spec_columns(plan)
n_spec = plan.n_obs - plan.n_blp;
relative_cols = unique([1, ceil(n_spec / 2), n_spec]);
spec_cols = plan.n_blp + relative_cols;
end

function local_write_dictionary_progress(progress_file, params, completed_chunk, n_chunks, ...
    completed_idx2, n_time, n_obs, work_save_file)
fid = fopen(progress_file, 'w');
if fid < 0
    warning('Could not write dictionary progress file: %s', progress_file);
    return;
end

cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, 'updated=%s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'spec_mode=%s\n', params.spec_mode);
fprintf(fid, 'precision=%s\n', params.precision);
fprintf(fid, 'chunk_size=%d\n', params.chunk_size);
fprintf(fid, 'completed_chunk=%d\n', completed_chunk);
fprintf(fid, 'n_chunks=%d\n', n_chunks);
fprintf(fid, 'completed_idx2=%d\n', completed_idx2);
fprintf(fid, 'n_time=%d\n', n_time);
fprintf(fid, 'n_obs=%d\n', n_obs);
fprintf(fid, 'work_save_file=%s\n', work_save_file);
clear cleanup_obj
end

function local_validate_dictionary_output(save_file, D, plan, check_content)
% Catch interrupted MAT-file writes before the temporary file is promoted.
if nargin < 4 || isempty(check_content)
    check_content = true;
end

params = struct('precision', '');
info = whos('-file', save_file, 'obs');
if ~isempty(info)
    params.precision = info.class;
end
local_validate_dictionary_shape(save_file, D, plan, params);

if ~check_content
    return;
end

M = matfile(save_file);
for k = 1:numel(D.session_ids)
    idx1 = D.session_start_idx(k);
    idx2 = D.session_end_idx(k);
    n_probe = min(1000, idx2 - idx1 + 1);
    probe_start = idx1 + floor((idx2 - idx1 + 1 - n_probe) / 2);
    probe_idx = probe_start:(probe_start + n_probe - 1);
    raw_probe = local_read_matfile_probe(M, probe_idx, 1:plan.n_blp);
    if all(raw_probe(:) == 0)
        error(['Dictionary validation failed: session %d has all-zero raw BLP ' ...
            'probe values in %s. This usually means the write was interrupted.'], ...
            D.session_ids(k), save_file);
    end
    if plan.n_obs > plan.n_blp
        spec_cols = local_probe_spec_columns(plan);
        spec_probe = local_read_matfile_probe(M, probe_idx, spec_cols);
        if all(spec_probe(:) == 0)
            error(['Dictionary validation failed: session %d has all-zero spectrogram ' ...
                'probe values in %s. This usually means the write was interrupted.'], ...
            D.session_ids(k), save_file);
        end
    end
end
end

function values = local_read_matfile_probe(M, row_idx, col_idx)
% matfile only supports regularly-spaced subscripts. Fall back to reading
% one row or column at a time for sparse validation probes.
row_idx = row_idx(:)';
col_idx = col_idx(:)';

if isempty(row_idx) || isempty(col_idx)
    values = [];
    return;
end

if local_is_regular_index(row_idx) && local_is_regular_index(col_idx)
    values = M.obs(row_idx, col_idx);
    return;
end

if local_is_regular_index(row_idx)
    first_col = M.obs(row_idx, col_idx(1));
    values = zeros(numel(row_idx), numel(col_idx), class(first_col));
    values(:, 1) = first_col(:);
    for j = 2:numel(col_idx)
        col_values = M.obs(row_idx, col_idx(j));
        values(:, j) = col_values(:);
    end
    return;
end

if local_is_regular_index(col_idx)
    first_row = M.obs(row_idx(1), col_idx);
    values = zeros(numel(row_idx), numel(col_idx), class(first_row));
    values(1, :) = first_row;
    for i = 2:numel(row_idx)
        values(i, :) = M.obs(row_idx(i), col_idx);
    end
    return;
end

sample = M.obs(row_idx(1), col_idx(1));
values = zeros(numel(row_idx), numel(col_idx), class(sample));
values(1, 1) = sample;
for i = 1:numel(row_idx)
    for j = 1:numel(col_idx)
        if i == 1 && j == 1
            continue;
        end
        values(i, j) = M.obs(row_idx(i), col_idx(j));
    end
end
end

function tf = local_is_regular_index(idx)
idx = idx(:)';
if numel(idx) <= 2
    tf = true;
    return;
end

d = diff(idx);
tf = all(d == d(1));
end
