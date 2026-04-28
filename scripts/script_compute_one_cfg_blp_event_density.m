% Canonical single-dataset event-density stage entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_compute_one_cfg_blp_event_density.m');
%
% Optional controls, set before run(...) if needed:
%   force_density_recompute = false;
%   density_bin_sec = 2;
%   density_smooth_sigma_sec = 2;

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
if exist('force_density_recompute', 'var') && ~isempty(force_density_recompute)
    params.density_params.force_recompute = logical(force_density_recompute);
end
if exist('density_bin_sec', 'var') && ~isempty(density_bin_sec)
    params.density_params.bin_sec = double(density_bin_sec);
end
if exist('density_smooth_sigma_sec', 'var') && ~isempty(density_smooth_sigma_sec)
    params.density_params.smooth_sigma_sec = double(density_smooth_sigma_sec);
end

E = compute_blp_event_density( ...
    cfg, ...
    output_root, ...
    R, ...
    params.density_params, ...
    source_event_file);

fprintf('Saved event-density result to:\n  %s\n', E.save_file);
