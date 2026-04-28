% Canonical interactive BLP consensus-state pipeline entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_to_consensus_state_top_windows.m');
%
% Optional controls, set before run(...) if needed:
%   force_event_recompute = false;
%   force_density_recompute = false;
%   force_consensus_recompute = false;
%   force_summary_recompute = false;
%   force_window_recompute = false;
%   make_top_window_plots = true;
%   top_k = 30;
%   window_length_samples = 6000;
%
% Role:
%   - canonical single-dataset entry for the long-term BLP event/state line
%   - primary branch ends at consensus-state diversity top windows
%   - top-window plots are treated as the standard readout layer
%   - plotting requires a previously saved region-mean spectrogram file
%   - leaves D / R / E / C / S / W / Spec / P in the workspace
%
% This is the mainline script for the second pipeline.
% Event-diversity scripts are now treated as quick-check branches rather than
% the primary downstream path.
%
% For a fixed reference example, use:
%   scripts/script_run_e10gb1_to_consensus_state_top_windows.m
%
% For running multiple cfg_*.m definitions in one batch, use:
%   scripts/script_run_cfgs_to_consensus_state_top_windows.m

%% -------------------------
%  Required user input
%  -------------------------
if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

%% -------------------------
%  Repo and toolboxes
%  -------------------------
this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

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
%  User options
%  -------------------------
cfg_name = char(string(cfg_name));

params = build_blp_consensus_state_pipeline_params();
if exist('force_event_recompute', 'var') && ~isempty(force_event_recompute)
    params.event_params.force_recompute = logical(force_event_recompute);
end
if exist('force_density_recompute', 'var') && ~isempty(force_density_recompute)
    params.density_params.force_recompute = logical(force_density_recompute);
end
if exist('force_consensus_recompute', 'var') && ~isempty(force_consensus_recompute)
    params.consensus_params.force_recompute = logical(force_consensus_recompute);
end
if exist('force_summary_recompute', 'var') && ~isempty(force_summary_recompute)
    params.summary_params.force_recompute = logical(force_summary_recompute);
end
if exist('force_window_recompute', 'var') && ~isempty(force_window_recompute)
    params.window_params.force_recompute = logical(force_window_recompute);
end
if exist('top_k', 'var') && ~isempty(top_k)
    params.window_params.top_k = top_k;
end
if exist('window_length_samples', 'var') && ~isempty(window_length_samples)
    params.window_params.window_length_samples = window_length_samples;
end
if exist('make_top_window_plots', 'var') && ~isempty(make_top_window_plots)
    params.plot_params.enable = logical(make_top_window_plots);
end

%% -------------------------
%  Config and output root
%  -------------------------
cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

output_root = io_project.get_project_processed_root();

fprintf('Running BLP consensus-state pipeline for %s (%s)\n', cfg_name, cfg.dataset_id);
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
Spec = struct();
P = struct();
if params.plot_params.enable
    Spec = require_saved_blp_spectrogram_support(cfg, output_root);
    fprintf('Using saved spectrogram plotting support:\n  %s\n', Spec.abs_file);

    P = export_top_consensus_state_diversity_window_plots(cfg, output_root, W, params.plot_params);
    fprintf('Saved top-window plots to:\n  %s\n', P.save_dir);
end

disp('Top consensus-state-diversity windows:');
disp(W.top_windows_table);
