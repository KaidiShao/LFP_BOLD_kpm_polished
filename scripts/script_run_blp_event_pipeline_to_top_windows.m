this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if exist('filterSignal_mirror', 'file') ~= 2
    error('filterSignal_mirror.m was not found on the MATLAB path.');
end

if exist('find_peak_loc', 'file') ~= 2
    error('find_peak_loc.m was not found on the MATLAB path.');
end

close all force;

%% -------------------------
%  Dataset selection
%  -------------------------
cfg = cfg_E10gb1();   % Replace with your target cfg_*.m
output_root = 'E:\DataPons_processed\';

fprintf('Running BLP event pipeline for %s\n', cfg.dataset_id);
fprintf('Raw data root:\n  %s\n', cfg.raw_data_root);
fprintf('Output root:\n  %s\n', output_root);

%% -------------------------
%  Step 1. Load raw BLP
%  -------------------------
D = load_blp_dataset(cfg);
fprintf('Loaded concatenated raw data: [%d, %d]\n', size(D.data, 1), size(D.data, 2));

%% -------------------------
%  Step 2. Event detection
%  -------------------------
event_params = struct();
event_params.passband = [2, 15; 30, 90; 90, 190];
event_params.band_labels = {'theta', 'gamma', 'ripple'};
event_params.L_start_range = [151, 101, 51];
event_params.L_extract_range = [301, 201, 101];
event_params.ThresRatio_range = [3.5, 4, 4];
event_params.force_recompute = false;

R = compute_blp_bandpass_events(D, cfg, output_root, event_params);
fprintf('Saved event detection result to:\n  %s\n', R.save_file);

%% -------------------------
%  Step 3. Event density
%  -------------------------
density_params = struct();
density_params.bin_sec = 2;
density_params.smooth_sigma_sec = 2;
density_params.force_recompute = false;

E = compute_blp_event_density( ...
    cfg, ...
    output_root, ...
    R, ...
    density_params, ...
    R.save_file);

fprintf('Saved event-density result to:\n  %s\n', E.save_file);

%% -------------------------
%  Step 4. Consensus states
%  -------------------------
consensus_params = struct();
consensus_params.min_channel_count = [];
consensus_params.require_region_presence = false;
consensus_params.required_regions = {'hp', 'pl'};
consensus_params.force_recompute = false;

C = compute_blp_consensus_states( ...
    cfg, ...
    output_root, ...
    R, ...
    consensus_params, ...
    R.save_file);

fprintf('Saved consensus-state result to:\n  %s\n', C.save_file);

%% -------------------------
%  Step 5. Consensus summary
%  -------------------------
summary_params = struct();
summary_params.save_csv = true;
summary_params.force_recompute = false;

S = summarize_blp_consensus_state_types( ...
    cfg, ...
    output_root, ...
    C, ...
    summary_params, ...
    C.save_file);

fprintf('Saved consensus-state summary to:\n  %s\n', S.save_file);
if ~isempty(S.csv_file)
    fprintf('Saved consensus-state summary CSV to:\n  %s\n', S.csv_file);
end

%% -------------------------
%  Step 6. Top event-diversity windows
%  -------------------------
window_params = struct();
window_params.window_length_samples = 5000;
window_params.keep_partial_window = false;
window_params.top_k = 10;
window_params.save_csv = true;
window_params.force_recompute = false;

W = analyze_blp_consensus_event_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    window_params, ...
    C.save_file);

fprintf('Saved per-window event-diversity result to:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('Saved all-window CSV to:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('Saved top-window CSV to:\n  %s\n', W.top_csv_file);
end

disp('Top event-diversity windows:');
disp(W.top_windows_table);
