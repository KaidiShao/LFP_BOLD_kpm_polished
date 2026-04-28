% Canonical single-dataset event-detection stage entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_compute_one_cfg_blp_bandpass_events.m');
%
% Optional controls, set before run(...) if needed:
%   force_event_recompute = false;

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

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

cfg_fun = ['cfg_' char(string(cfg_name))];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
output_root = io_project.get_project_processed_root();

params = build_blp_consensus_state_pipeline_params();
if exist('force_event_recompute', 'var') && ~isempty(force_event_recompute)
    params.event_params.force_recompute = logical(force_event_recompute);
end

D = io_raw.load_blp_dataset(cfg);
R = compute_blp_bandpass_events(D, cfg, output_root, params.event_params);

fprintf('Saved event detection result to:\n  %s\n', R.save_file);
