% Example single-dataset BLP consensus-state mainline pipeline.
%
% Role:
%   - reference script for one concrete dataset
%   - easy to read top-to-bottom
%   - convenient when you want all intermediate variables left in workspace
%
% This is not the reusable entry. For running the same pipeline on one
% arbitrary cfg_*.m definition, use:
%   scripts/script_run_one_cfg_to_consensus_state_top_windows.m
%
% For running multiple cfg_*.m definitions in one batch, use:
%   scripts/script_run_cfgs_to_consensus_state_top_windows.m

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

%% -------------------------
%  Toolboxes
%  -------------------------
if exist('filterSignal_mirror', 'file') ~= 2
    error('filterSignal_mirror.m was not found on the MATLAB path.');
end

if exist('find_peak_loc', 'file') ~= 2
    error('find_peak_loc.m was not found on the MATLAB path.');
end

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

%% -------------------------
%  Dataset and output root
%  -------------------------
cfg = cfg_E10gb1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

output_root = io_project.get_project_processed_root();
params = build_blp_consensus_state_pipeline_params();

fprintf('Running BLP consensus-state pipeline for %s\n', cfg.dataset_id);
fprintf('Raw data root:\n  %s\n', cfg.raw_data_root);
fprintf('Output root:\n  %s\n', output_root);

%% -------------------------
%  Step 1. Load raw BLP
%  -------------------------
D = io_raw.load_blp_dataset(cfg);
fprintf('Loaded concatenated raw data: [%d, %d]\n', size(D.data, 1), size(D.data, 2));

%% -------------------------
%  Step 2. Event detection
%  -------------------------
R = compute_blp_bandpass_events(D, cfg, output_root, params.event_params);
fprintf('Saved event detection result to:\n  %s\n', R.save_file);

%% -------------------------
%  Step 3. Event density
%  -------------------------
E = compute_blp_event_density(cfg, output_root, R, params.density_params, R.save_file);
fprintf('Saved event-density result to:\n  %s\n', E.save_file);

%% -------------------------
%  Step 4. Consensus states
%  -------------------------
C = compute_blp_consensus_states(cfg, output_root, R, params.consensus_params, R.save_file);
fprintf('Saved consensus-state result to:\n  %s\n', C.save_file);

%% -------------------------
%  Step 5. Consensus summary
%  -------------------------
S = summarize_blp_consensus_state_types(cfg, output_root, C, params.summary_params, C.save_file);
fprintf('Saved consensus-state summary to:\n  %s\n', S.save_file);
if ~isempty(S.csv_file)
    fprintf('Saved consensus-state summary CSV to:\n  %s\n', S.csv_file);
end

%% -------------------------
%  Step 6. Consensus-state diversity windows
%  -------------------------
W = analyze_blp_consensus_state_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    params.window_params, ...
    C.save_file);

fprintf('Saved consensus-state-diversity result to:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('Saved all-window CSV to:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('Saved top-window CSV to:\n  %s\n', W.top_csv_file);
end

%% -------------------------
%  Step 7. Require saved spectrogram support, then draw top-window plots
%  -------------------------
Spec = require_saved_blp_spectrogram_support(cfg, output_root);
fprintf('Using saved spectrogram plotting support:\n  %s\n', Spec.abs_file);

P = export_top_consensus_state_diversity_window_plots(cfg, output_root, W, params.plot_params);
fprintf('Saved top-window plots to:\n  %s\n', P.save_dir);

disp('Top consensus-state-diversity windows:');
disp(W.top_windows_table);
