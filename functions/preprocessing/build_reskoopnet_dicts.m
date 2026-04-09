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
    output_root = 'D:\DataPons_processed\';
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
pad_sec = 20;
pad_mode = 'mirror';

if isfield(cfg, 'spectrogram')
    if isfield(cfg.spectrogram, 'pad_sec')
        pad_sec = cfg.spectrogram.pad_sec;
    end
    if isfield(cfg.spectrogram, 'pad_mode')
        pad_mode = cfg.spectrogram.pad_mode;
    end
end

abs_file = find_regionmean_spectrogram_file_by_kind(cfg, output_root, pad_mode, pad_sec, 'abs');
complex_file = find_regionmean_spectrogram_file_by_kind(cfg, output_root, pad_mode, pad_sec, 'complex');

%% =========================
%  Read spectrogram metadata
%  =========================
meta = load(abs_file, 'freqs', 'regions');
freqs = double(meta.freqs(:));
regions = meta.regions;

low_idx = find(freqs >= 0 & freqs <= params.low_full_max_hz);
high_idx = find(freqs > params.low_full_max_hz & freqs <= params.high_max_hz);
high_groups = make_consecutive_groups(high_idx, params.high_group_size);

if isempty(low_idx) && isempty(high_groups)
    error('No spectrogram frequency bins were selected.');
end

%% =========================
%  Build observable metadata
%  =========================
obs_info = build_observable_info(D, cfg, regions, freqs, low_idx, high_groups, params.spec_mode);
n_obs = height(obs_info);
n_time = size(D.data, 1);

fprintf('Total time points   : %d\n', n_time);
fprintf('Total observables   : %d\n', n_obs);
fprintf('BLP observables     : %d\n', numel(D.selected_channels));
fprintf('Spectrogram regions : %d\n', numel(regions));

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
if strcmp(params.spec_mode, 'abs')
    Mspec = matfile(abs_file);
    spec_var_name = 'tmpall_mean_abs';
else
    Mspec = matfile(complex_file);
    spec_var_name = 'tmpall_mean_complex';
end

spec_size = size(Mspec, spec_var_name);
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

    obs_chunk = build_observable_chunk(D, Mspec, spec_var_name, idx, ...
        low_idx, high_groups, regions, params);

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
session_start_idx = D.session_start_idx;
session_end_idx = D.session_end_idx;
border_idx = D.border_idx;
selected_channels = D.selected_channels;
channel_sites = cfg.channels.sites(D.selected_channels);

save(save_file, ...
    'obs_info', 'params', 'freqs', 'regions', ...
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
dict.obs_info = obs_info;
end

function obs_chunk = build_observable_chunk(D, Mspec, spec_var_name, idx, ...
    low_idx, high_groups, regions, params)
% Build one time chunk of observables.

raw_chunk = double(D.data(idx, :));   % time x n_channels
n_time_chunk = size(raw_chunk, 1);
n_blp = size(raw_chunk, 2);
n_regions = numel(regions);
n_spec_per_region = numel(low_idx) + numel(high_groups);

if strcmp(params.spec_mode, 'abs')
    n_obs = n_blp + n_regions * n_spec_per_region;
else
    n_obs = n_blp + n_regions * n_spec_per_region * 2;
end

obs_chunk = zeros(n_time_chunk, n_obs);

% Keep all BLP channels
obs_chunk(:, 1:n_blp) = raw_chunk;

% Read spectrogram chunk
spec_chunk = Mspec.(spec_var_name)(:, idx, :);   % freq x time x region
spec_chunk = double(spec_chunk);

col0 = n_blp + 1;

for r = 1:n_regions
    spec_region = spec_chunk(:, :, r);   % freq x time

    if strcmp(params.spec_mode, 'abs')
        feat = aggregate_spectrogram_features(spec_region, low_idx, high_groups);
        n_feat = size(feat, 2);

        obs_chunk(:, col0:col0+n_feat-1) = feat;
        col0 = col0 + n_feat;

    else
        feat_real = aggregate_spectrogram_features(real(spec_region), low_idx, high_groups);
        feat_imag = aggregate_spectrogram_features(imag(spec_region), low_idx, high_groups);

        n_feat = size(feat_real, 2);

        obs_chunk(:, col0:col0+n_feat-1) = feat_real;
        col0 = col0 + n_feat;

        obs_chunk(:, col0:col0+n_feat-1) = feat_imag;
        col0 = col0 + n_feat;
    end
end
end

function feat = aggregate_spectrogram_features(spec_region, low_idx, high_groups)
% Aggregate a spectrogram matrix into ResKoopNet-friendly observables.
%
% Input
%   spec_region  freq x time
%
% Output
%   feat         time x n_features

n_time = size(spec_region, 2);
n_low = numel(low_idx);
n_high = numel(high_groups);

feat = zeros(n_time, n_low + n_high);

% Keep all low-frequency bins
if n_low > 0
    feat(:, 1:n_low) = transpose(spec_region(low_idx, :));
end

% Average high-frequency bins in small groups
col0 = n_low + 1;
for g = 1:n_high
    idxg = high_groups{g};
    feat(:, col0) = transpose(mean(spec_region(idxg, :), 1));
    col0 = col0 + 1;
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

function obs_info = build_observable_info(D, cfg, regions, freqs, low_idx, high_groups, spec_mode)
% Build a table describing every observable dimension.

dim_id = [];
source = {};
name = {};
channel_idx = [];
channel_site = {};
region_idx = [];
region_label = {};
part = {};
freq_idx_start = [];
freq_idx_end = [];
freq_start_hz = [];
freq_end_hz = [];
aggregation = {};

row_id = 0;

%% -------------------------
%  Raw BLP observables
%  -------------------------
sites = cfg.channels.sites(D.selected_channels);

for c = 1:numel(D.selected_channels)
    row_id = row_id + 1;

    dim_id(end+1, 1) = row_id;
    source{end+1, 1} = 'blp';
    name{end+1, 1} = sprintf('blp_ch%02d_%s', D.selected_channels(c), sites{c});

    channel_idx(end+1, 1) = D.selected_channels(c);
    channel_site{end+1, 1} = sites{c};

    region_idx(end+1, 1) = NaN;
    region_label{end+1, 1} = '';
    part{end+1, 1} = 'raw';

    freq_idx_start(end+1, 1) = NaN;
    freq_idx_end(end+1, 1) = NaN;
    freq_start_hz(end+1, 1) = NaN;
    freq_end_hz(end+1, 1) = NaN;

    aggregation{end+1, 1} = 'raw';
end

%% -------------------------
%  Spectrogram observables
%  -------------------------
if strcmp(spec_mode, 'abs')
    parts_to_add = {'abs'};
else
    % Match build_observable_chunk(): for complex_split, all real-valued
    % features are written first, followed by the corresponding imag block.
    parts_to_add = {'real', 'imag'};
end

for r = 1:numel(regions)
    for p = 1:numel(parts_to_add)
        % Low-frequency single-bin observables
        for i = 1:numel(low_idx)
            fi = low_idx(i);
            part_name = parts_to_add{p};

            row_id = row_id + 1;

            dim_id(end+1, 1) = row_id;
            source{end+1, 1} = 'spectrogram';
            name{end+1, 1} = sprintf('spec_%s_%s_f%.3f', part_name, regions{r}, freqs(fi));

            channel_idx(end+1, 1) = NaN;
            channel_site{end+1, 1} = '';

            region_idx(end+1, 1) = r;
            region_label{end+1, 1} = regions{r};
            part{end+1, 1} = part_name;

            freq_idx_start(end+1, 1) = fi;
            freq_idx_end(end+1, 1) = fi;
            freq_start_hz(end+1, 1) = freqs(fi);
            freq_end_hz(end+1, 1) = freqs(fi);

            aggregation{end+1, 1} = 'single_bin';
        end
        % High-frequency grouped observables
        for g = 1:numel(high_groups)
            idxg = high_groups{g};
            fi1 = idxg(1);
            fi2 = idxg(end);
            part_name = parts_to_add{p};

            row_id = row_id + 1;

            dim_id(end+1, 1) = row_id;
            source{end+1, 1} = 'spectrogram';
            name{end+1, 1} = sprintf('spec_%s_%s_f%.3f_%.3f', ...
                part_name, regions{r}, freqs(fi1), freqs(fi2));

            channel_idx(end+1, 1) = NaN;
            channel_site{end+1, 1} = '';

            region_idx(end+1, 1) = r;
            region_label{end+1, 1} = regions{r};
            part{end+1, 1} = part_name;

            freq_idx_start(end+1, 1) = fi1;
            freq_idx_end(end+1, 1) = fi2;
            freq_start_hz(end+1, 1) = freqs(fi1);
            freq_end_hz(end+1, 1) = freqs(fi2);

            aggregation{end+1, 1} = sprintf('mean_%d_bins', numel(idxg));
        end
    end
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

function spec_file = find_regionmean_spectrogram_file_by_kind(cfg, output_root, pad_mode, pad_sec, kind)
% Find a saved region-mean spectrogram file by kind.
%
% kind = 'abs' or 'complex'

if strcmpi(pad_mode, 'mirror')
    pad_tag = sprintf('_mirrorpad_%gs', pad_sec);
else
    pad_tag = '_nopad';
end
pad_tag = strrep(pad_tag, '.', 'p');

if strcmp(kind, 'abs')
    target_name = [cfg.file_stem, pad_tag, '_regionmean_spectrograms_abs.mat'];
elseif strcmp(kind, 'complex')
    target_name = [cfg.file_stem, pad_tag, '_regionmean_spectrograms_complex.mat'];
else
    error('Unsupported kind.');
end

search_dirs = { ...
    fullfile(output_root, cfg.file_stem, 'spectrograms'), ...
    fullfile(output_root, 'spectrograms', cfg.file_stem)};

for i = 1:numel(search_dirs)
    f = fullfile(search_dirs{i}, target_name);
    if exist(f, 'file') == 2
        spec_file = f;
        return;
    end
end

error('Spectrogram file not found: %s', target_name);
end
