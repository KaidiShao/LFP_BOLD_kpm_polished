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
%   run_log_file = '';

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ' ...
        'cfg_names = {''E10gb1'', ''E10gH1''};']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_names = cellstr(string(cfg_names(:)).');

if ~exist('run_log_file', 'var') || isempty(run_log_file)
    run_log_file = '';
end

if ~isempty(run_log_file)
    fid = fopen(run_log_file, 'w');
    if fid < 0
        error('Could not open run_log_file for writing: %s', run_log_file);
    end
    fclose(fid);
end

single_cfg_script = 'scripts/script_build_one_cfg_bold_observables.m';
batch_outputs = struct([]);

local_emit(run_log_file, 'Running BOLD observable batch for %d dataset(s)\n', numel(cfg_names));

for i = 1:numel(cfg_names)
    cfg_name = cfg_names{i};
    local_emit(run_log_file, '\n============================================================\n');
    local_emit(run_log_file, 'Building BOLD observables for %s\n', cfg_name);

    run_cmd = sprintf('cfg_name = ''%s''; run(''%s'');', ...
        strrep(cfg_name, '''', ''''''), single_cfg_script);
    run_text = evalc(run_cmd);
    if isempty(run_log_file)
        fprintf('%s', run_text);
    else
        local_append_text(run_log_file, run_text);
    end

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

local_emit(run_log_file, '\nFinished BOLD observable batch.\n');


function local_emit(run_log_file, varargin)
txt = sprintf(varargin{:});
local_append_text(run_log_file, txt);
if isempty(run_log_file)
    fprintf('%s', txt);
end
end


function local_append_text(run_log_file, txt)
if isempty(run_log_file) || isempty(txt)
    return;
end

fid = fopen(run_log_file, 'a');
if fid < 0
    error('Could not open run_log_file for append: %s', run_log_file);
end

cleanup_obj = onCleanup(@() fclose(fid));
fprintf(fid, '%s', txt);
clear cleanup_obj
end
