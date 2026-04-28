function x = read_blp_data_slice(D, idx)
% Read a time slice from an in-memory or on-disk BLP dataset struct.

if nargin < 2 || isempty(idx)
    x = [];
    return;
end

if isfield(D, 'data_storage') && strcmpi(D.data_storage, 'disk')
    if ~isfield(D, 'data_file') || exist(D.data_file, 'file') ~= 2
        error('BLP cache file was not found: %s', D.data_file);
    end

    if isfield(D, 'data_var') && ~isempty(D.data_var)
        data_var = D.data_var;
    else
        data_var = 'data';
    end

    M = matfile(D.data_file);
    x = M.(data_var)(idx, :);
elseif isfield(D, 'data_storage') && strcmpi(D.data_storage, 'source')
    idx = idx(:)';
    if any(diff(idx) ~= 1)
        error('Source-backed BLP reads currently require contiguous sample indices.');
    end

    if ~all(isfield(D, {'session_start_idx', 'session_end_idx', 'session_ids', ...
            'selected_channels', 'raw_data_root', 'data_subfolder', 'file_stem'}))
        error('Source-backed BLP dataset metadata is incomplete.');
    end

    idx1 = idx(1);
    idx2 = idx(end);
    session_mask = D.session_start_idx <= idx2 & D.session_end_idx >= idx1;
    session_ids = find(session_mask);

    x_parts = cell(numel(session_ids), 1);

    for i = 1:numel(session_ids)
        k = session_ids(i);
        sid = D.session_ids(k);
        [x_session, ~] = local_load_source_session(D, sid);

        take1 = max(idx1, D.session_start_idx(k));
        take2 = min(idx2, D.session_end_idx(k));

        local_idx1 = take1 - D.session_start_idx(k) + 1;
        local_idx2 = take2 - D.session_start_idx(k) + 1;

        x_parts{i} = x_session(local_idx1:local_idx2, :);
    end

    x = cat(1, x_parts{:});
else
    x = D.data(idx, :);
end
end

function [x, dx] = local_load_source_session(D, sid)
fname = fullfile(D.raw_data_root, D.data_subfolder, ...
    sprintf('%s_%04d_blp.mat', D.file_stem, sid));

S = load(fname, 'blp');
x = S.blp.dat(:, D.selected_channels, 1);
x = reshape(x, size(S.blp.dat, 1), numel(D.selected_channels));

if nargout > 1
    dx = S.blp.dx;
end
end
