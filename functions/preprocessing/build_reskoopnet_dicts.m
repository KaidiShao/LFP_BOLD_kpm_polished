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
save_dir = fullfile(output_root, cfg.file_stem, 'reskoopnet_dictionary');
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

if exist(save_file, 'file') == 2
    delete(save_file);
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
Mout = matfile(save_file, 'Writable', true);

if strcmp(params.precision, 'single')
    Mout.obs(n_time, n_obs) = single(0);
elseif strcmp(params.precision, 'double')
    Mout.obs(n_time, n_obs) = double(0);
else
    error('Unsupported params.precision.');
end

%% =========================
%  Chunk loop
%  =========================
chunk_size = params.chunk_size;
n_chunks = ceil(n_time / chunk_size);

for k = 1:n_chunks
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

save(save_file, ...
    'obs_info', 'params', 'freqs', 'regions', ...
    'dx', 'dt', 'fs', 'sampling_period', 'sample_period', 'sampling_frequency', ...
    'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'channel_sites', ...
    '-append');

writetable(obs_info, info_csv_file);

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
