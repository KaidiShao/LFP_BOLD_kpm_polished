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

[all_session_ids, all_session_channels] = ...
    io_utils.collect_included_session_entries(cfg, cfg.channels.selected_all);

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
[out.session_start_idx, out.session_end_idx, out.border_idx] = ...
    io_utils.build_session_index_metadata(session_lengths);
out.selected_channels = selected_channels_ref;

if isfield(cfg, 'channels') && isfield(cfg.channels, 'sites') && ~isempty(selected_channels_ref)
    out.channel_sites = cfg.channels.sites(selected_channels_ref);
else
    out.channel_sites = repmat({''}, numel(selected_channels_ref), 1);
end

[out.dx, out.fs] = io_utils.resolve_uniform_dx(session_dx, ...
    'Spike sampling period dx is inconsistent across sessions. Check out.session_dx.');
end
