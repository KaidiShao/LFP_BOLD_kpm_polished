function V = analyze_blp_consensus_event_variability_windows(cfg, output_root, C, params, source_consensus_file)
% Analyze consensus-state variability in non-overlapping sample windows.
%
% Inputs
%   cfg                   Dataset config struct
%   output_root           Root folder for processed outputs
%   C                     Loaded consensus-state result struct
%   params                Struct with window-analysis parameters
%   source_consensus_file Optional source consensus-state file path
%
% Output
%   V                     Struct containing per-window variability metrics
%                         and a top-ranked table of variable windows

%% =========================
%  Input defaults
%  =========================
if nargin < 2 || isempty(output_root)
    output_root = 'D:\DataPons_processed\';
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
required_fields = {'state_catalog', 'state_code_by_time', 'session_ids', ...
    'session_lengths', 'session_dx', 'session_start_idx'};
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

%% =========================
%  Build per-window metrics
%  =========================
window_table_cells = cell(0, 1);

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

    session_global_start = session_start_idx(k);

    for w = 1:n_windows
        local_start = (w - 1) * window_length_samples + 1;
        local_end = min(w * window_length_samples, n_samples);
        global_start = session_global_start + local_start - 1;
        global_end = session_global_start + local_end - 1;

        code_seg = state_code_by_time(global_start:global_end);
        active_mask = code_seg > 0;
        active_codes = code_seg(active_mask);

        labeled_sample_count = nnz(active_mask);
        unlabeled_sample_count = numel(code_seg) - labeled_sample_count;
        labeled_fraction = labeled_sample_count / max(numel(code_seg), 1);

        if numel(code_seg) >= 2
            transition_count = nnz(diff(code_seg) ~= 0);
        else
            transition_count = 0;
        end

        transition_rate_hz = transition_count / ((local_end - local_start + 1) * session_dx(k));

        if isempty(active_codes)
            active_state_richness = 0;
            state_entropy = 0;
            normalized_state_entropy = 0;
            dominant_state = "";
            state_sample_counts = zeros(1, n_states);
        else
            state_sample_counts = zeros(1, n_states);
            for s = 1:n_states
                state_sample_counts(s) = nnz(active_codes == state_codes(s));
            end

            active_state_richness = nnz(state_sample_counts > 0);
            [state_entropy, normalized_state_entropy] = compute_entropy(state_sample_counts, n_states);

            dominant_idx = find(state_sample_counts == max(state_sample_counts), 1, 'first');
            dominant_state = state_labels(dominant_idx);
        end

        row = table( ...
            session_ids(k), ...
            w, ...
            global_start, ...
            global_end, ...
            local_start, ...
            local_end, ...
            local_end - local_start + 1, ...
            (local_end - local_start + 1) * session_dx(k), ...
            labeled_sample_count, ...
            unlabeled_sample_count, ...
            labeled_fraction, ...
            transition_count, ...
            transition_rate_hz, ...
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
                'transition_count', ...
                'transition_rate_hz', ...
                'active_state_richness', ...
                'state_entropy', ...
                'normalized_state_entropy', ...
                'dominant_state'});

        for s = 1:n_states
            var_name = matlab.lang.makeValidName(sprintf('%s_sample_count', char(state_labels(s))));
            row.(var_name) = state_sample_counts(s);
        end

        window_table_cells{end+1, 1} = row; %#ok<AGROW>
    end
end

if isempty(window_table_cells)
    window_table = table();
else
    window_table = vertcat(window_table_cells{:});
end

%% =========================
%  Rank variable windows
%  =========================
if isempty(window_table)
    top_windows_table = table();
else
    positive_mask = window_table.labeled_sample_count > 0;
    ranked_table = window_table(positive_mask, :);

    if ~isempty(ranked_table)
        ranked_table = sortrows(ranked_table, ...
            {'transition_rate_hz', 'normalized_state_entropy', 'active_state_richness', 'labeled_fraction', 'transition_count'}, ...
            {'descend', 'descend', 'descend', 'descend', 'descend'});

        top_windows_table = ranked_table(1:min(top_k, height(ranked_table)), :);
        top_windows_table.variability_rank = transpose(1:height(top_windows_table));
        top_windows_table = movevars(top_windows_table, 'variability_rank', 'Before', 1);
    else
        top_windows_table = table();
    end
end

%% =========================
%  Prepare save path
%  =========================
save_dir = fullfile(output_root, cfg.file_stem, 'event_variability_windows');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_tag = build_save_tag(cfg.file_stem, window_length_samples);
save_file = fullfile(save_dir, [save_tag, '.mat']);
csv_file = fullfile(save_dir, [save_tag, '.csv']);
top_csv_file = fullfile(save_dir, [save_tag, '_top.csv']);

if exist(save_file, 'file') == 2 && ~params.force_recompute
    L = load(save_file);
    V = L.V;
    return;
end

%% =========================
%  Pack output
%  =========================
V = struct();
V.save_file = save_file;
V.csv_file = csv_file;
V.top_csv_file = top_csv_file;
V.source_consensus_file = source_consensus_file;
V.dataset_id = cfg.dataset_id;
V.file_stem = cfg.file_stem;

V.state_codes = state_codes;
V.state_labels = state_labels(:);
V.window_length_samples = window_length_samples;
V.keep_partial_window = logical(params.keep_partial_window);
V.top_k = top_k;

V.window_table = window_table;
V.top_windows_table = top_windows_table;
V.params = params;

save(save_file, 'V', '-v7.3');

if params.save_csv
    writetable(window_table, csv_file);
    writetable(top_windows_table, top_csv_file);
end
end


function params = apply_default_params(params)
if ~isfield(params, 'window_length_samples')
    params.window_length_samples = 5000;
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


function tag = build_save_tag(file_stem, window_length_samples)
tag = sprintf('%s_event_variability_windows_%dsamp', file_stem, window_length_samples);
end
