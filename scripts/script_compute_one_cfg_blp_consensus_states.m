% Canonical single-dataset consensus-state stage entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_compute_one_cfg_blp_consensus_states.m');
%
% Optional controls, set before run(...) if needed:
%   force_consensus_recompute = false;

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_fun = ['cfg_' char(string(cfg_name))];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
output_root = io_project.get_project_processed_root();

[R, source_event_file] = io_results.load_event_results(cfg, output_root, []);

params = build_blp_consensus_state_pipeline_params();
if exist('force_consensus_recompute', 'var') && ~isempty(force_consensus_recompute)
    params.consensus_params.force_recompute = logical(force_consensus_recompute);
end

C = compute_blp_consensus_states( ...
    cfg, ...
    output_root, ...
    R, ...
    params.consensus_params, ...
    source_event_file);

fprintf('Saved consensus-state result to:\n  %s\n', C.save_file);
