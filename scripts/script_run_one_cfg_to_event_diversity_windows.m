% Canonical interactive event-diversity quick-check entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_to_event_diversity_windows.m');
%
% Optional controls, set before run(...) if needed:
%   force_event_recompute = false;
%   force_density_recompute = false;
%   force_consensus_recompute = false;
%   force_summary_recompute = false;
%   force_window_recompute = false;
%   top_k = 10;
%   window_length_samples = 6000;
%   window_mode = 'global';
%
% Role:
%   - canonical single-dataset quick-check branch
%   - useful for early event-line diagnostics
%   - leaves D / R / E / C / S / W in the workspace
%
% For the long-term mainline, use:
%   scripts/script_run_one_cfg_to_consensus_state_top_windows.m
%
% For a fixed e10gb1 quick-check example, use:
%   scripts/script_run_e10gb1_to_event_diversity_windows.m
%
% For batch execution across cfg_*.m entries, use:
%   scripts/script_run_cfgs_to_event_diversity_windows.m

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

close all force;

%% -------------------------
%  User options
%  -------------------------
cfg_name = char(string(cfg_name));

params = build_blp_event_diversity_pipeline_params();
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
if exist('window_mode', 'var') && ~isempty(window_mode)
    params.window_params.window_mode = char(string(window_mode));
end

%% -------------------------
%  Config and output root
%  -------------------------
cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
output_root = io_project.get_project_processed_root();

fprintf('Running BLP event-diversity quick-check for %s (%s)\n', cfg_name, cfg.dataset_id);
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
%  Step 6. Quick-check event-diversity windows
%  -------------------------
W = analyze_blp_consensus_event_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    params.window_params, ...
    C.save_file);

fprintf('Saved event-diversity result to:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('Saved all-window CSV to:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('Saved top-window CSV to:\n  %s\n', W.top_csv_file);
end

disp('Top event-diversity windows:');
disp(W.top_windows_table);
