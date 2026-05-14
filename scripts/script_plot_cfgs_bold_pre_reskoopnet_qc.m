% Batch driver for BOLD pre-ResKoopNet QC.
%
% Usage from MATLAB:
%   cfg_names = {'E10gb1', 'E10gH1'};
%   run('scripts/script_plot_cfgs_bold_pre_reskoopnet_qc.m');

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ' ...
        'cfg_names = {''E10gb1'', ''E10gH1''};']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

clear info
clc;

cfg_names = cellstr(string(cfg_names(:)).');
single_cfg_script = 'scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m';
qc_batch_outputs = struct([]);

for i = 1:numel(cfg_names)
    cfg_name = cfg_names{i};
    fprintf('\n============================================================\n');
    fprintf('Plotting BOLD QC for %s\n', cfg_name);
    run(sprintf('cfg_name = ''%s''; run(''%s'');', ...
        strrep(cfg_name, '''', ''''''), single_cfg_script));

    qc_batch_outputs(i).cfg_name = cfg_name; %#ok<SAGROW>
    if exist('cfg', 'var') && isstruct(cfg) && isfield(cfg, 'dataset_id')
        qc_batch_outputs(i).dataset_id = cfg.dataset_id; %#ok<SAGROW>
    else
        qc_batch_outputs(i).dataset_id = ''; %#ok<SAGROW>
    end
    if exist('qc_outputs', 'var')
        qc_batch_outputs(i).outputs = qc_outputs; %#ok<SAGROW>
    else
        qc_batch_outputs(i).outputs = struct([]); %#ok<SAGROW>
    end
end
