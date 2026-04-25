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
        get_project_processed_root(), cfg.file_stem, 'raw_cache', ...
        [cfg.file_stem, '_selected_blp.mat']);
end

data_dir = fullfile(cfg.raw_data_root, cfg.data_subfolder);

if opts.metadata_only
    opts.cache_to_disk = false;
end

if opts.cache_to_disk
    cache_dir = fileparts(opts.cache_file);
    if exist(cache_dir, 'dir') ~= 7
        mkdir(cache_dir);
    end

    if ~opts.force_recompute && exist(opts.cache_file, 'file') == 2 ...
            && local_has_complete_cached_result(opts.cache_file)
        out = local_load_cached_result(opts.cache_file);
        return;
    end

    if exist(opts.cache_file, 'file') == 2
        delete(opts.cache_file);
    end

    Mcache = matfile(opts.cache_file, 'Writable', true);
    write_row0 = 1;
else
    data_cells = {};
end

session_ids = [];
session_lengths = [];
session_dx = [];
selected_channels_ref = [];

% Collect all included session IDs first
all_session_ids = [];
all_session_channels = {};

for i = 1:numel(cfg.sessions)
    if ~cfg.sessions(i).include
        continue;
    end

    this_session_ids = cfg.sessions(i).session_id;

    if isempty(cfg.sessions(i).selected_channels)
        this_channels = cfg.channels.selected_all;
    else
        this_channels = cfg.sessions(i).selected_channels;
    end

    for j = 1:numel(this_session_ids)
        all_session_ids(end+1, 1) = this_session_ids(j);
        all_session_channels{end+1, 1} = this_channels;
    end
end

n_sessions = numel(all_session_ids);

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

    if opts.cache_to_disk
        idx_write = write_row0:(write_row0 + size(x, 1) - 1);
        Mcache.data(idx_write, 1:numel(this_channels)) = x;
        write_row0 = idx_write(end) + 1;
    elseif ~opts.metadata_only
        data_cells{end+1, 1} = x;
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

out = struct();

if opts.cache_to_disk
    out.data = [];
    out.data_storage = 'disk';
    out.data_file = opts.cache_file;
    out.data_var = 'data';
    out.n_time = sum(session_lengths);
elseif opts.metadata_only
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
out.session_end_idx = cumsum(session_lengths);
out.session_start_idx = [1; out.session_end_idx(1:end-1) + 1];
out.border_idx = out.session_end_idx(1:end-1);
out.selected_channels = selected_channels_ref;
out.raw_data_root = cfg.raw_data_root;
out.data_subfolder = cfg.data_subfolder;
out.file_stem = cfg.file_stem;

if opts.cache_to_disk
    session_end_idx = out.session_end_idx;
    session_start_idx = out.session_start_idx;
    border_idx = out.border_idx;
    selected_channels = out.selected_channels;
    n_time = out.n_time;

    save(opts.cache_file, ...
        'session_ids', 'session_lengths', 'session_dx', ...
        'session_start_idx', 'session_end_idx', 'border_idx', ...
        'selected_channels', 'n_time', '-append');
end

if ~isempty(session_dx)
    dx0 = session_dx(1);
    if all(abs(session_dx - dx0) < 1e-12)
        out.dx = dx0;
        out.fs = 1 / dx0;
    else
        out.dx = [];
        out.fs = [];
        warning('Sampling period dx is inconsistent across sessions. Check out.session_dx.');
    end
end
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

if ~isempty(out.session_dx)
    dx0 = out.session_dx(1);
    if all(abs(out.session_dx - dx0) < 1e-12)
        out.dx = dx0;
        out.fs = 1 / dx0;
    else
        out.dx = [];
        out.fs = [];
        warning('Sampling period dx is inconsistent across sessions. Check out.session_dx.');
    end
end
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
