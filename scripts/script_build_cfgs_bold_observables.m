% Batch driver for the BOLD observable pipeline.
%
% Usage from MATLAB:
%   cfg_names = {'E10gb1', 'E10gH1'};
%   run('scripts/script_build_cfgs_bold_observables.m');
%
% Optional controls, set before run(...) if needed:
%   observable_modes = {'eleHP', 'HP', 'roi_mean'};
%   n_svd_components = 50;
%   n_svd100_components = 100;
%   snapshot_lag = 1;
%   save_precision = 'single';
%   force_recompute = false;
%   pre = struct(...);
%   bold_output_root = '';

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ' ...
        'cfg_names = {''E10gb1'', ''E10gH1''};']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_names = cellstr(string(cfg_names(:)).');

single_cfg_script = 'scripts/script_build_one_cfg_bold_observables.m';
batch_outputs = struct([]);

fprintf('Running BOLD observable batch for %d dataset(s)\n', numel(cfg_names));

for i = 1:numel(cfg_names)
    cfg_name = cfg_names{i};
    fprintf('\n============================================================\n');
    fprintf('Building BOLD observables for %s\n', cfg_name);
    run(sprintf('cfg_name = ''%s''; run(''%s'');', ...
        strrep(cfg_name, '''', ''''''), single_cfg_script));

    batch_outputs(i).cfg_name = cfg_name; %#ok<SAGROW>
    if exist('cfg', 'var') && isstruct(cfg) && isfield(cfg, 'dataset_id')
        batch_outputs(i).dataset_id = cfg.dataset_id; %#ok<SAGROW>
    else
        batch_outputs(i).dataset_id = ''; %#ok<SAGROW>
    end
    if exist('bold_outputs', 'var')
        batch_outputs(i).outputs = bold_outputs; %#ok<SAGROW>
    else
        batch_outputs(i).outputs = struct([]); %#ok<SAGROW>
    end
end

fprintf('\nFinished BOLD observable batch.\n');
