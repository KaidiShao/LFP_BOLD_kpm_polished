% Batch driver for the canonical pipeline 5 single-dataset script.
%
% Usage from MATLAB:
%   cfg_names = {'E10fV1','E10gH1'};
%   run('scripts/script_run_cfgs_blp_eigenfunction_reduction.m');

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ', ...
        'cfg_names = {''E10fV1'',''E10gH1''};']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg_names = cellstr(string(cfg_names(:)).');
batch_outputs = cell(numel(cfg_names), 1);

for i = 1:numel(cfg_names)
    cfg_name = cfg_names{i}; %#ok<NASGU>
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s\n', i, numel(cfg_names), cfg_name);
    run(fullfile(this_script_dir, 'script_run_one_cfg_blp_eigenfunction_reduction.m'));
    batch_outputs{i} = manifest; %#ok<NASGU>
end
