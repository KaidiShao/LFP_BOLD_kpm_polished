% Batch driver for canonical BLP eigenfunction postprocessing runs.
%
% Usage from MATLAB:
%   cfg_names = {'E10fV1', 'E10gH1'};
%   run('scripts/script_run_cfgs_blp_eigenfunction_postprocessing.m');

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ' ...
        'cfg_names = {''E10fV1'', ''E10gH1''};']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg_names = cellstr(string(cfg_names(:)).');
batch_manifests = struct([]);

for i_cfg = 1:numel(cfg_names)
    cfg_name = cfg_names{i_cfg}; %#ok<NASGU>
    run(fullfile(this_script_dir, 'script_run_one_cfg_blp_eigenfunction_postprocessing.m'));

    batch_manifests(i_cfg).cfg_name = cfg_name; %#ok<SAGROW>
    if exist('manifest', 'var')
        batch_manifests(i_cfg).manifest = manifest; %#ok<SAGROW>
    else
        batch_manifests(i_cfg).manifest = struct(); %#ok<SAGROW>
    end
end
