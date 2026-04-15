function S = summarize_blp_consensus_state_types(cfg, output_root, C, params, source_consensus_file)
% Summarize consensus-state types by count and duration.
%
% Inputs
%   cfg                   Dataset config struct
%   output_root           Root folder for processed outputs
%   C                     Loaded consensus-state result struct
%   params                Struct with summary parameters
%   source_consensus_file Optional source consensus-state file path
%
% Output
%   S                     Struct containing summary tables and metadata

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
required_fields = {'state_catalog', 'state_code_by_time', 'state_windows', ...
    'session_lengths', 'session_dx', 'session_start_idx', 'session_end_idx'};
for i = 1:numel(required_fields)
    if ~isfield(C, required_fields{i})
        error('The loaded consensus-state result is missing field "%s".', required_fields{i});
    end
end

if isempty(source_consensus_file) && isfield(C, 'save_file') && ~isempty(C.save_file)
    source_consensus_file = C.save_file;
end

state_catalog = C.state_catalog(:);
n_states = numel(state_catalog);

session_lengths = double(C.session_lengths(:));
session_dx = double(C.session_dx(:));
session_start_idx = double(C.session_start_idx(:));
session_end_idx = double(C.session_end_idx(:));

if numel(session_lengths) ~= numel(session_dx) || ...
        numel(session_lengths) ~= numel(session_start_idx) || ...
        numel(session_lengths) ~= numel(session_end_idx)
    error('C session-length, dx, and index metadata must all have the same length.');
end

recording_duration_sec = sum(session_lengths .* session_dx);
total_sample_count = sum(session_lengths);
labeled_sample_count = nnz(double(C.state_code_by_time) > 0);
unlabeled_sample_count = total_sample_count - labeled_sample_count;
labeled_duration_sec = compute_mask_duration_sec(double(C.state_code_by_time) > 0, session_dx, session_start_idx, session_end_idx);
unlabeled_duration_sec = max(0, recording_duration_sec - labeled_duration_sec);

if isempty(C.state_windows)
    window_codes = zeros(0, 1);
    window_durations = zeros(0, 1);
else
    window_codes = double([C.state_windows.state_code]).';
    window_durations = double([C.state_windows.duration_sec]).';
end

total_window_count = numel(window_codes);

%% =========================
%  Build summary table
%  =========================
state_code = zeros(n_states, 1);
state_label = strings(n_states, 1);
window_count = zeros(n_states, 1);
window_fraction = zeros(n_states, 1);
sample_count = zeros(n_states, 1);
duration_sec = zeros(n_states, 1);
occupancy_fraction_total = zeros(n_states, 1);
occupancy_fraction_labeled = zeros(n_states, 1);
mean_window_duration_sec = nan(n_states, 1);
median_window_duration_sec = nan(n_states, 1);
min_window_duration_sec = nan(n_states, 1);
max_window_duration_sec = nan(n_states, 1);

for s = 1:n_states
    code = double(state_catalog(s).code);
    label = string(state_catalog(s).label);
    idx = (window_codes == code);
    durations = window_durations(idx);

    state_code(s) = code;
    state_label(s) = label;
    window_count(s) = nnz(idx);

    if total_window_count > 0
        window_fraction(s) = window_count(s) / total_window_count;
    end

    if isfield(C, 'state_sample_count') && numel(C.state_sample_count) >= s
        sample_count(s) = double(C.state_sample_count(s));
    else
        sample_count(s) = nnz(double(C.state_code_by_time) == code);
    end

    if isfield(C, 'state_duration_sec') && numel(C.state_duration_sec) >= s
        duration_sec(s) = double(C.state_duration_sec(s));
    else
        duration_sec(s) = compute_mask_duration_sec(double(C.state_code_by_time) == code, ...
            session_dx, session_start_idx, session_end_idx);
    end

    if recording_duration_sec > 0
        occupancy_fraction_total(s) = duration_sec(s) / recording_duration_sec;
    end

    if labeled_duration_sec > 0
        occupancy_fraction_labeled(s) = duration_sec(s) / labeled_duration_sec;
    end

    if ~isempty(durations)
        mean_window_duration_sec(s) = mean(durations);
        median_window_duration_sec(s) = median(durations);
        min_window_duration_sec(s) = min(durations);
        max_window_duration_sec(s) = max(durations);
    end
end

summary_table = table( ...
    state_code, ...
    state_label, ...
    window_count, ...
    window_fraction, ...
    sample_count, ...
    duration_sec, ...
    occupancy_fraction_total, ...
    occupancy_fraction_labeled, ...
    mean_window_duration_sec, ...
    median_window_duration_sec, ...
    min_window_duration_sec, ...
    max_window_duration_sec);

%% =========================
%  Prepare save path
%  =========================
save_dir = fullfile(output_root, cfg.file_stem, 'consensus_state_summary');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_tag = build_save_tag(cfg.file_stem);
save_file = fullfile(save_dir, [save_tag, '.mat']);
csv_file = fullfile(save_dir, [save_tag, '.csv']);

if exist(save_file, 'file') == 2 && ~params.force_recompute
    L = load(save_file);
    S = L.S;
    return;
end

%% =========================
%  Pack output
%  =========================
S = struct();
S.save_file = save_file;
S.csv_file = csv_file;
S.source_consensus_file = source_consensus_file;
S.dataset_id = cfg.dataset_id;
S.file_stem = cfg.file_stem;

S.recording_duration_sec = recording_duration_sec;
S.total_sample_count = total_sample_count;
S.labeled_sample_count = labeled_sample_count;
S.unlabeled_sample_count = unlabeled_sample_count;
S.labeled_duration_sec = labeled_duration_sec;
S.unlabeled_duration_sec = unlabeled_duration_sec;
S.total_window_count = total_window_count;

S.summary_table = summary_table;
S.params = params;

save(save_file, 'S', '-v7.3');

if params.save_csv
    writetable(summary_table, csv_file);
end
end


function params = apply_default_params(params)
if ~isfield(params, 'save_csv')
    params.save_csv = true;
end

if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end
end


function duration_sec = compute_mask_duration_sec(mask, session_dx, session_start_idx, session_end_idx)
if isempty(mask)
    duration_sec = 0;
    return;
end

duration_sec = 0;

for k = 1:numel(session_dx)
    idx1 = session_start_idx(k);
    idx2 = session_end_idx(k);
    duration_sec = duration_sec + nnz(mask(idx1:idx2)) * double(session_dx(k));
end
end


function tag = build_save_tag(file_stem)
tag = sprintf('%s_consensus_state_type_summary', file_stem);
end
