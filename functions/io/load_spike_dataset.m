function out = load_spike_dataset(cfg)
%LOAD_SPIKE_DATASET Load selected spike sessions and concatenate them in time.
%
% This mirrors load_blp_dataset.m, but reads session files of the form
%   <file_stem>_<session>_spkt.mat
% containing a top-level struct named "Spkt".

data_dir = fullfile(cfg.raw_data_root, 'spkt');

data_cells = {};
session_ids = [];
session_lengths = [];
session_dx = [];
selected_channels_ref = [];

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
        all_session_ids(end+1, 1) = this_session_ids(j); %#ok<AGROW>
        all_session_channels{end+1, 1} = this_channels; %#ok<AGROW>
    end
end

n_sessions = numel(all_session_ids);

for k = 1:n_sessions
    sid = all_session_ids(k);
    this_channels = all_session_channels{k};

    fprintf('[%d/%d] Loading spike session %04d ...\n', k, n_sessions, sid);

    if isempty(selected_channels_ref)
        selected_channels_ref = this_channels;
    elseif ~isequal(selected_channels_ref, this_channels)
        error('Selected spike channels must be the same across concatenated sessions.');
    end

    fname = fullfile(data_dir, sprintf('%s_%04d_spkt.mat', cfg.file_stem, sid));
    S = load(fname, 'Spkt');

    if ~isfield(S, 'Spkt')
        error('File %s does not contain variable Spkt.', fname);
    end

    x = double(S.Spkt.dat(:, this_channels));

    data_cells{end+1, 1} = x; %#ok<AGROW>
    session_ids(end+1, 1) = sid; %#ok<AGROW>
    session_lengths(end+1, 1) = size(x, 1); %#ok<AGROW>
    session_dx(end+1, 1) = S.Spkt.dx; %#ok<AGROW>

    fprintf('         done, size = [%d, %d], dx = %.12g\n', ...
        size(x, 1), size(x, 2), S.Spkt.dx);
end

out = struct();
out.data = cat(1, data_cells{:});
out.data_cells = data_cells;
out.session_ids = session_ids;
out.session_lengths = session_lengths;
out.session_dx = session_dx;
out.session_end_idx = cumsum(session_lengths);
out.session_start_idx = [1; out.session_end_idx(1:end-1) + 1];
out.border_idx = out.session_end_idx(1:end-1);
out.selected_channels = selected_channels_ref;

if isfield(cfg, 'channels') && isfield(cfg.channels, 'sites') && ~isempty(selected_channels_ref)
    out.channel_sites = cfg.channels.sites(selected_channels_ref);
else
    out.channel_sites = repmat({''}, numel(selected_channels_ref), 1);
end

if ~isempty(session_dx)
    dx0 = session_dx(1);
    if all(abs(session_dx - dx0) < 1e-12)
        out.dx = dx0;
        out.fs = 1 / dx0;
    else
        out.dx = [];
        out.fs = [];
        warning('Spike sampling period dx is inconsistent across sessions. Check out.session_dx.');
    end
end
end
