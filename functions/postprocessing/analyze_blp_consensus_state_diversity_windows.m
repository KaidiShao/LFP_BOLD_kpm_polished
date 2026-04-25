function W = analyze_blp_consensus_state_diversity_windows(cfg, output_root, C, params, source_consensus_file)
% Analyze consensus-state diversity in non-overlapping sample windows.
%
% Inputs
%   cfg                   Dataset config struct
%   output_root           Root folder for processed outputs
%   C                     Loaded consensus-state result struct
%   params                Struct with window-analysis parameters
%   source_consensus_file Optional source consensus-state file path
%
% Output
%   W                     Struct containing per-window state-diversity
%                         metrics and a top-ranked table of diverse windows

%% =========================
%  Input defaults
%  =========================
if nargin < 2 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 3
    error('C must be provided. Load it first with load_consensus_state_results.');
end

if nargin < 4
    params = struct();
end

if nargin < 5
    source_consensus_file = '';
end

params = apply_default_params(params);

%% =========================
%  Validate consensus-state results
%  =========================
required_fields = {'state_catalog', 'state_windows', 'state_code_by_time', ...
    'session_ids', 'session_lengths', 'session_dx', 'session_start_idx'};
for i = 1:numel(required_fields)
    if ~isfield(C, required_fields{i})
        error('The loaded consensus-state result is missing field "%s".', required_fields{i});
    end
end

if isempty(source_consensus_file) && isfield(C, 'save_file') && ~isempty(C.save_file)
    source_consensus_file = C.save_file;
end

window_length_samples = double(params.window_length_samples);
top_k = double(params.top_k);
window_mode = validatestring(char(string(params.window_mode)), {'session', 'global'});

if ~isscalar(window_length_samples) || window_length_samples <= 0 || window_length_samples ~= round(window_length_samples)
    error('params.window_length_samples must be a positive integer.');
end

if ~isscalar(top_k) || top_k <= 0 || top_k ~= round(top_k)
    error('params.top_k must be a positive integer.');
end

state_catalog = C.state_catalog(:);
state_labels = string({state_catalog.label});
state_codes = double([state_catalog.code]).';
n_states = numel(state_codes);

session_ids = double(C.session_ids(:));
session_lengths = double(C.session_lengths(:));
session_dx = double(C.session_dx(:));
session_start_idx = double(C.session_start_idx(:));
state_code_by_time = double(C.state_code_by_time(:));
n_sessions = numel(session_ids);
source_consensus_file_signature = build_file_signature(source_consensus_file);

%% =========================
%  Prepare save path / cache
%  =========================
save_dir = fullfile(output_root, cfg.file_stem, 'consensus_state_diversity_windows');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_tag = build_save_tag(cfg.file_stem, window_length_samples, window_mode);
save_file = fullfile(save_dir, [save_tag, '.mat']);
csv_file = fullfile(save_dir, [save_tag, '.csv']);
top_csv_file = fullfile(save_dir, [save_tag, '_top.csv']);

if exist(save_file, 'file') == 2 && ~params.force_recompute
    L = load(save_file);
    if isfield(L, 'W') && can_reuse_saved_result(L.W, params, source_consensus_file_signature)
        W = L.W;
        maybe_write_csv_outputs(W, params);
        return;
    end
end

%% =========================
%  Collect state-window centers
%  =========================
state_window_centers = zeros(0, 1);
state_window_codes = zeros(0, 1);
state_window_session_idx = zeros(0, 1);

if ~isempty(C.state_windows)
    all_windows = C.state_windows(:);
    n_all = numel(all_windows);
    state_window_centers = zeros(n_all, 1);
    state_window_codes = zeros(n_all, 1);
    state_window_session_idx = zeros(n_all, 1);

    for i = 1:n_all
        state_window_centers(i) = round(mean(double(all_windows(i).win_global)));
        state_window_codes(i) = double(all_windows(i).state_code);
        state_window_session_idx(i) = double(all_windows(i).session_idx);
    end
end

%% =========================
%  Build per-window metrics
%  =========================
window_table_cells = cell(0, 1);

if strcmp(window_mode, 'session')
    for k = 1:n_sessions
        n_samples = session_lengths(k);
        n_full_windows = floor(n_samples / window_length_samples);
        has_partial = mod(n_samples, window_length_samples) > 0;

        if params.keep_partial_window
            n_windows = n_full_windows + double(has_partial);
        else
            n_windows = n_full_windows;
        end

        if n_windows == 0
            continue;
        end

        session_counts = zeros(n_windows, n_states);
        state_mask = state_window_session_idx == k;
        centers_global = state_window_centers(state_mask);
        codes_this = state_window_codes(state_mask);

        if ~isempty(centers_global)
            centers_local = centers_global - session_start_idx(k) + 1;

            if params.keep_partial_window
                valid_event = centers_local >= 1 & centers_local <= n_samples;
            else
                valid_event = centers_local >= 1 & centers_local <= n_full_windows * window_length_samples;
            end

            centers_local = centers_local(valid_event);
            codes_this = codes_this(valid_event);

            for s = 1:n_states
                code = state_codes(s);
                centers_state = centers_local(codes_this == code);
                if isempty(centers_state)
                    continue;
                end

                window_idx = ceil(centers_state / window_length_samples);
                session_counts(:, s) = session_counts(:, s) + accumarray(window_idx, 1, [n_windows, 1], @sum, 0);
            end
        end

        for w = 1:n_windows
            local_start = (w - 1) * window_length_samples + 1;
            local_end = min(w * window_length_samples, n_samples);
            global_start = session_start_idx(k) + local_start - 1;
            global_end = session_start_idx(k) + local_end - 1;

            row = local_build_state_diversity_row( ...
                session_ids(k), ...
                w, ...
                global_start, ...
                global_end, ...
                local_start, ...
                local_end, ...
                (global_end - global_start + 1) * session_dx(k), ...
                session_counts(w, :), ...
                state_labels, ...
                state_codes, ...
                n_states, ...
                state_code_by_time);
            window_table_cells{end+1, 1} = row; %#ok<AGROW>
        end
    end
else
    total_samples = numel(state_code_by_time);
    n_full_windows = floor(total_samples / window_length_samples);
    has_partial = mod(total_samples, window_length_samples) > 0;

    if params.keep_partial_window
        n_windows = n_full_windows + double(has_partial);
    else
        n_windows = n_full_windows;
    end

    global_counts = zeros(n_windows, n_states);

    if ~isempty(state_window_centers)
        if params.keep_partial_window
            valid_event = state_window_centers >= 1 & state_window_centers <= total_samples;
        else
            valid_event = state_window_centers >= 1 & state_window_centers <= n_full_windows * window_length_samples;
        end

        centers_global = state_window_centers(valid_event);
        codes_global = state_window_codes(valid_event);

        for s = 1:n_states
            code = state_codes(s);
            centers_state = centers_global(codes_global == code);
            if isempty(centers_state)
                continue;
            end

            window_idx = ceil(centers_state / window_length_samples);
            global_counts(:, s) = global_counts(:, s) + accumarray(window_idx, 1, [n_windows, 1], @sum, 0);
        end
    end

    for w = 1:n_windows
        global_start = (w - 1) * window_length_samples + 1;
        global_end = min(w * window_length_samples, total_samples);
        span_info = resolve_global_sample_span( ...
            global_start, global_end, session_ids, session_lengths, session_dx, session_start_idx);

        row = local_build_state_diversity_row( ...
            span_info.start_session_id, ...
            NaN, ...
            global_start, ...
            global_end, ...
            span_info.start_local_idx, ...
            span_info.end_local_idx, ...
            span_info.duration_sec, ...
            global_counts(w, :), ...
            state_labels, ...
            state_codes, ...
            n_states, ...
            state_code_by_time);
        row.global_window_idx = w;
        row.start_session_id = span_info.start_session_id;
        row.end_session_id = span_info.end_session_id;
        row.crosses_session_boundary = span_info.crosses_session_boundary;
        window_table_cells{end+1, 1} = row; %#ok<AGROW>
    end
end

if isempty(window_table_cells)
    window_table = table();
else
    window_table = vertcat(window_table_cells{:});
end

%% =========================
%  Rank state-diverse windows
%  =========================
if isempty(window_table)
    top_windows_table = table();
else
    positive_mask = window_table.total_state_window_count > 0;
    ranked_table = window_table(positive_mask, :);

    if ~isempty(ranked_table)
        ranked_table = sortrows(ranked_table, ...
            {'active_state_richness', 'normalized_state_entropy', 'total_state_window_count', 'labeled_fraction', 'state_entropy'}, ...
            {'descend', 'descend', 'descend', 'descend', 'descend'});

        top_windows_table = ranked_table(1:min(top_k, height(ranked_table)), :);
        top_windows_table.state_diversity_rank = transpose(1:height(top_windows_table));
        top_windows_table = movevars(top_windows_table, 'state_diversity_rank', 'Before', 1);
    else
        top_windows_table = table();
    end
end

W = struct();
W.save_file = save_file;
W.csv_file = csv_file;
W.top_csv_file = top_csv_file;
W.source_consensus_file = source_consensus_file;
W.source_consensus_file_signature = source_consensus_file_signature;
W.dataset_id = cfg.dataset_id;
W.file_stem = cfg.file_stem;

W.state_codes = state_codes;
W.state_labels = state_labels(:);
W.window_length_samples = window_length_samples;
W.window_mode = window_mode;
W.keep_partial_window = logical(params.keep_partial_window);
W.top_k = top_k;

W.window_table = window_table;
W.top_windows_table = top_windows_table;
W.params = params;

save(save_file, 'W', '-v7.3');
maybe_write_csv_outputs(W, params);
end


function params = apply_default_params(params)
if ~isfield(params, 'window_length_samples')
    params.window_length_samples = 5000;
end

if ~isfield(params, 'window_mode') || isempty(params.window_mode)
    params.window_mode = 'session';
end

if ~isfield(params, 'keep_partial_window')
    params.keep_partial_window = false;
end

if ~isfield(params, 'top_k')
    params.top_k = 10;
end

if ~isfield(params, 'save_csv')
    params.save_csv = true;
end

if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end
end


function [entropy_val, normalized_entropy] = compute_entropy(counts, n_states)
counts = double(counts(:));
total_count = sum(counts);

if total_count <= 0
    entropy_val = 0;
    normalized_entropy = 0;
    return;
end

p = counts / total_count;
p = p(p > 0);
entropy_val = -sum(p .* log(p));

if n_states <= 1
    normalized_entropy = 0;
else
    normalized_entropy = entropy_val / log(n_states);
end
end


function row = local_build_state_diversity_row( ...
    session_id, window_idx_in_session, global_start, global_end, local_start, local_end, ...
    duration_sec, state_window_counts, state_labels, state_codes, n_states, state_code_by_time)
state_window_counts = double(state_window_counts(:)).';
total_state_window_count = sum(state_window_counts);
active_state_richness = nnz(state_window_counts > 0);
[state_entropy, normalized_state_entropy] = compute_entropy(state_window_counts, n_states);

if total_state_window_count > 0
    dominant_idx = find(state_window_counts == max(state_window_counts), 1, 'first');
    dominant_state = state_labels(dominant_idx);
else
    dominant_state = "";
end

code_seg = state_code_by_time(global_start:global_end);
active_mask = code_seg > 0;
labeled_sample_count = nnz(active_mask);
unlabeled_sample_count = numel(code_seg) - labeled_sample_count;
labeled_fraction = labeled_sample_count / max(numel(code_seg), 1);

state_sample_counts = zeros(1, n_states);
for s = 1:n_states
    state_sample_counts(s) = nnz(code_seg == state_codes(s));
end

row = table( ...
    session_id, ...
    window_idx_in_session, ...
    global_start, ...
    global_end, ...
    local_start, ...
    local_end, ...
    global_end - global_start + 1, ...
    duration_sec, ...
    labeled_sample_count, ...
    unlabeled_sample_count, ...
    labeled_fraction, ...
    total_state_window_count, ...
    active_state_richness, ...
    state_entropy, ...
    normalized_state_entropy, ...
    string(dominant_state), ...
    'VariableNames', { ...
        'session_id', ...
        'window_idx_in_session', ...
        'global_start_idx', ...
        'global_end_idx', ...
        'session_local_start_idx', ...
        'session_local_end_idx', ...
        'window_samples', ...
        'window_duration_sec', ...
        'labeled_sample_count', ...
        'unlabeled_sample_count', ...
        'labeled_fraction', ...
        'total_state_window_count', ...
        'active_state_richness', ...
        'state_entropy', ...
        'normalized_state_entropy', ...
        'dominant_state'});

for s = 1:n_states
    count_var = matlab.lang.makeValidName(sprintf('%s_count', char(state_labels(s))));
    row.(count_var) = state_window_counts(s);

    sample_var = matlab.lang.makeValidName(sprintf('%s_sample_count', char(state_labels(s))));
    row.(sample_var) = state_sample_counts(s);
end
end


function tag = build_save_tag(file_stem, window_length_samples, window_mode)
tag = sprintf('%s_consensus_state_diversity_windows_%dsamp', file_stem, window_length_samples);
if strcmp(window_mode, 'global')
    tag = sprintf('%s_globalwin', tag);
end
end


function tf = can_reuse_saved_result(W, params, source_consensus_file_signature)
tf = isstruct(W);
if ~tf
    return;
end

required_fields = {'window_length_samples', 'window_mode', 'keep_partial_window', ...
    'top_k', 'source_consensus_file_signature'};
for i = 1:numel(required_fields)
    if ~isfield(W, required_fields{i}) || isempty(W.(required_fields{i}))
        tf = false;
        return;
    end
end

tf = double(W.window_length_samples) == double(params.window_length_samples) && ...
    strcmpi(char(string(W.window_mode)), char(string(params.window_mode))) && ...
    logical(W.keep_partial_window) == logical(params.keep_partial_window) && ...
    double(W.top_k) == double(params.top_k);
if ~tf
    return;
end

tf = file_signature_matches(W.source_consensus_file_signature, source_consensus_file_signature);
end


function maybe_write_csv_outputs(W, params)
if ~params.save_csv
    return;
end

writetable(W.window_table, W.csv_file);
writetable(W.top_windows_table, W.top_csv_file);
end
