function out = load_blp_dataset(cfg)
% Load selected BLP sessions and concatenate them in time.

data_dir = fullfile(cfg.raw_data_root, cfg.data_subfolder);

data_cells = {};
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

    data_cells{end+1, 1} = x;
    session_ids(end+1, 1) = sid;
    session_lengths(end+1, 1) = size(x, 1);
    session_dx(end+1, 1) = S.blp.dx;

    fprintf('         done, size = [%d, %d], dx = %.12g\n', ...
        size(x, 1), size(x, 2), S.blp.dx);
end

out = struct();
out.data = cat(1, data_cells{:});
out.session_ids = session_ids;
out.session_lengths = session_lengths;
out.session_dx = session_dx;
out.session_end_idx = cumsum(session_lengths);
out.session_start_idx = [1; out.session_end_idx(1:end-1) + 1];
out.border_idx = out.session_end_idx(1:end-1);
out.selected_channels = selected_channels_ref;

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