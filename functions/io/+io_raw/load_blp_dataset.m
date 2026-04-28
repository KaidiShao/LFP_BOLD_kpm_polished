function out = load_blp_dataset(cfg, opts)
% Load selected BLP sessions and concatenate them in time.
%
% Optional opts fields:
%   .cache_to_disk   false (default) | true
%   .force_recompute false (default) | true
%   .cache_file      explicit MAT-file path for on-disk caching
%   .metadata_only   false (default) | true

if nargin < 2
    opts = struct();
end

% Resolve loader options
if ~isfield(opts, 'metadata_only') || isempty(opts.metadata_only)
    opts.metadata_only = false;
end

if ~isfield(opts, 'cache_to_disk') || isempty(opts.cache_to_disk)
    opts.cache_to_disk = false;
end

if ~isfield(opts, 'force_recompute') || isempty(opts.force_recompute)
    opts.force_recompute = false;
end

if ~isfield(opts, 'cache_file') || isempty(opts.cache_file)
    opts.cache_file = fullfile( ...
        io_project.get_project_processed_root(), cfg.file_stem, 'raw_cache', ...
        [cfg.file_stem, '_selected_blp.mat']);
end

data_dir = fullfile(cfg.raw_data_root, cfg.data_subfolder);

if opts.metadata_only
    opts.cache_to_disk = false;
end

if opts.metadata_only
    storage_mode = 'source';
elseif opts.cache_to_disk
    storage_mode = 'disk';
else
    storage_mode = 'memory';
end

% Prepare disk cache before loading any sessions
Mcache = [];
write_row0 = [];
if strcmp(storage_mode, 'disk')
    [cached_out, Mcache, write_row0] = local_prepare_disk_cache(opts);
    if ~isempty(cached_out)
        out = cached_out;
        return;
    end
end

session_ids = [];
session_lengths = [];
session_dx = [];
selected_channels_ref = [];

% Collect all included session IDs first
[all_session_ids, all_session_channels] = ...
    io_utils.collect_included_session_entries(cfg, cfg.channels.selected_all);

n_sessions = numel(all_session_ids);
if strcmp(storage_mode, 'memory')
    data_cells = cell(n_sessions, 1);
end

% Load sessions
for k = 1:n_sessions
    sid = all_session_ids(k);
    this_channels = all_session_channels{k};

    fprintf('[%d/%d] Loading session %04d ...\n', k, n_sessions, sid);

    if isempty(selected_channels_ref)
        selected_channels_ref = this_channels;
    elseif ~isequal(selected_channels_ref, this_channels)
        error('Selected channels must be the same across concatenated sessions.');
    end

    fname = fullfile(data_dir, sprintf('%s_%04d_blp.mat', cfg.file_stem, sid));

    S = load(fname, 'blp');

    % Use only the first page in the 3rd dimension
    x = S.blp.dat(:, this_channels, 1);
    x = reshape(x, size(S.blp.dat, 1), numel(this_channels));

    if strcmp(storage_mode, 'disk')
        idx_write = write_row0:(write_row0 + size(x, 1) - 1);
        Mcache.data(idx_write, 1:numel(this_channels)) = x;
        write_row0 = idx_write(end) + 1;
    elseif strcmp(storage_mode, 'memory')
        data_cells{k} = x;
    else
        % metadata_only keeps only session bookkeeping in memory
    end

    session_ids(end+1, 1) = sid;
    session_lengths(end+1, 1) = size(x, 1);
    session_dx(end+1, 1) = S.blp.dx;

    fprintf('         done, size = [%d, %d], dx = %.12g\n', ...
        size(x, 1), size(x, 2), S.blp.dx);

    clear S x
end

% Assemble output struct
out = struct();

if strcmp(storage_mode, 'disk')
    out.data = [];
    out.data_storage = 'disk';
    out.data_file = opts.cache_file;
    out.data_var = 'data';
    out.n_time = sum(session_lengths);
elseif strcmp(storage_mode, 'source')
    out.data = [];
    out.data_storage = 'source';
    out.n_time = sum(session_lengths);
else
    out.data = cat(1, data_cells{:});
    out.data_storage = 'memory';
    out.n_time = size(out.data, 1);
end

out.session_ids = session_ids;
out.session_lengths = session_lengths;
out.session_dx = session_dx;
[out.session_start_idx, out.session_end_idx, out.border_idx] = ...
    io_utils.build_session_index_metadata(session_lengths);
out.selected_channels = selected_channels_ref;
out.raw_data_root = cfg.raw_data_root;
out.data_subfolder = cfg.data_subfolder;
out.file_stem = cfg.file_stem;

% Append metadata once the cache file contains the full data matrix
if strcmp(storage_mode, 'disk')
    local_finalize_disk_cache(opts.cache_file, out);
end

[out.dx, out.fs] = io_utils.resolve_uniform_dx(session_dx, ...
    'Sampling period dx is inconsistent across sessions. Check out.session_dx.');
end

function [cached_out, Mcache, write_row0] = local_prepare_disk_cache(opts)
cached_out = [];
Mcache = [];
write_row0 = [];

cache_dir = fileparts(opts.cache_file);
if exist(cache_dir, 'dir') ~= 7
    mkdir(cache_dir);
end

if ~opts.force_recompute && exist(opts.cache_file, 'file') == 2 ...
        && local_has_complete_cached_result(opts.cache_file)
    cached_out = local_load_cached_result(opts.cache_file);
    return;
end

if exist(opts.cache_file, 'file') == 2
    delete(opts.cache_file);
end

Mcache = matfile(opts.cache_file, 'Writable', true);
write_row0 = 1;
end

function local_finalize_disk_cache(cache_file, out)
session_ids = out.session_ids;
session_lengths = out.session_lengths;
session_dx = out.session_dx;
session_start_idx = out.session_start_idx;
session_end_idx = out.session_end_idx;
border_idx = out.border_idx;
selected_channels = out.selected_channels;
n_time = out.n_time;

save(cache_file, ...
    'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'n_time', '-append');
end

function out = local_load_cached_result(cache_file)
meta_fields = { ...
    'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'n_time'};

S = load(cache_file, meta_fields{:});

out = S;
out.data = [];
out.data_storage = 'disk';
out.data_file = cache_file;
out.data_var = 'data';

if ~isfield(out, 'n_time') || isempty(out.n_time)
    info = whos('-file', cache_file, 'data');
    out.n_time = info.size(1);
end

[out.dx, out.fs] = io_utils.resolve_uniform_dx(out.session_dx, ...
    'Sampling period dx is inconsistent across sessions. Check out.session_dx.');
end

function tf = local_has_complete_cached_result(cache_file)
tf = false;

try
    vars = who('-file', cache_file);
catch
    return;
end

required_vars = { ...
    'data', 'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'n_time'};

tf = all(ismember(required_vars, vars));
end
