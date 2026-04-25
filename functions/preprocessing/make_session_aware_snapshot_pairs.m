function S = make_session_aware_snapshot_pairs(X_or_struct, lag, params)
% Build snapshot pairs without crossing session borders.
%
% Input can be a numeric matrix plus params.session_* metadata, or a struct
% containing .data, .session_start_idx, .session_end_idx, and .session_ids.

if nargin < 2 || isempty(lag)
    lag = 1;
end
if nargin < 3
    params = struct();
end

if lag < 1 || lag ~= round(lag)
    error('lag must be a positive integer.');
end

if isstruct(X_or_struct)
    Xall = X_or_struct.data;
    session_start_idx = X_or_struct.session_start_idx;
    session_end_idx = X_or_struct.session_end_idx;
    session_ids = X_or_struct.session_ids;
else
    Xall = X_or_struct;
    required = {'session_start_idx', 'session_end_idx', 'session_ids'};
    for i = 1:numel(required)
        if ~isfield(params, required{i})
            error('params.%s is required when passing a numeric matrix.', required{i});
        end
    end
    session_start_idx = params.session_start_idx;
    session_end_idx = params.session_end_idx;
    session_ids = params.session_ids;
end

if isempty(Xall)
    error('Input data is empty.');
end

if ~isfield(params, 'drop_initial') || isempty(params.drop_initial)
    params.drop_initial = 0;
end
if ~isfield(params, 'train_ratio') || isempty(params.train_ratio)
    params.train_ratio = [];
end
if ~isfield(params, 'rng_seed') || isempty(params.rng_seed)
    params.rng_seed = [];
end

valid_idx_x = [];
valid_idx_y = [];
session_idx = [];
session_id = [];

for k = 1:numel(session_ids)
    idx1 = session_start_idx(k) + params.drop_initial;
    idx2 = session_end_idx(k);

    ix = (idx1:(idx2 - lag)).';
    iy = ix + lag;

    if isempty(ix)
        continue;
    end

    valid_idx_x = [valid_idx_x; ix]; %#ok<AGROW>
    valid_idx_y = [valid_idx_y; iy]; %#ok<AGROW>
    session_idx = [session_idx; repmat(k, numel(ix), 1)]; %#ok<AGROW>
    session_id = [session_id; repmat(session_ids(k), numel(ix), 1)]; %#ok<AGROW>
end

S = struct();
S.X = Xall(valid_idx_x, :);
S.Y = Xall(valid_idx_y, :);
S.valid_idx_x = valid_idx_x;
S.valid_idx_y = valid_idx_y;
S.session_idx = session_idx;
S.session_id = session_id;
S.lag = lag;
S.drop_initial = params.drop_initial;
S.n_pairs = numel(valid_idx_x);

if ~isempty(params.train_ratio)
    if params.train_ratio <= 0 || params.train_ratio >= 1
        error('params.train_ratio must be in (0, 1).');
    end

    if ~isempty(params.rng_seed)
        rng(params.rng_seed);
    end

    order = randperm(S.n_pairs).';
    n_train = floor(params.train_ratio * S.n_pairs);
    S.train_idx = order(1:n_train);
    S.valid_idx = order(n_train+1:end);
    S.X_train = S.X(S.train_idx, :);
    S.Y_train = S.Y(S.train_idx, :);
    S.X_valid = S.X(S.valid_idx, :);
    S.Y_valid = S.Y(S.valid_idx, :);
end
end
