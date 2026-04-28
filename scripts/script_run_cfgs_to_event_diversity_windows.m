% Batch driver for the BLP event-diversity quick-check branch.
%
% Usage from MATLAB:
%   cfg_names = {'E10fV1', 'E10gH1'};
%   run('scripts/script_run_cfgs_to_event_diversity_windows.m');
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
%   run_log_file = '';
%
% Role:
%   - batch driver for the canonical event-diversity quick-check branch
%   - runs the canonical single-dataset quick-check script once per cfg_name

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ' ...
        'cfg_names = {''E10fV1'', ''E10gH1''};']);
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

cfg_names = cellstr(string(cfg_names));

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

batch_outputs = struct([]);
single_cfg_script = 'scripts/script_run_one_cfg_to_event_diversity_windows.m';

local_emit(run_log_file, 'Running BLP event-diversity quick-check for %d dataset(s)\n', numel(cfg_names));

for i = 1:numel(cfg_names)
    cfg_name = cfg_names{i};
    cfg_fun = ['cfg_' cfg_name];
    if exist(cfg_fun, 'file') ~= 2
        error('Config function not found on path: %s', cfg_fun);
    end

    cfg = feval(cfg_fun);
    local_emit(run_log_file, '\n============================================================\n');
    local_emit(run_log_file, 'Running BLP event-diversity quick-check for %s (%s)\n', cfg_name, cfg.dataset_id);
    local_emit(run_log_file, 'Raw data root:\n  %s\n', cfg.raw_data_root);

    run_cmd = sprintf('cfg_name = ''%s''; run(''%s'');', ...
        strrep(cfg_name, '''', ''''''), single_cfg_script);
    run_text = evalc(run_cmd);
    if isempty(run_log_file)
        fprintf('%s', run_text);
    else
        local_append_text(run_log_file, run_text);
    end

    batch_outputs(i).cfg_name = cfg_name; %#ok<SAGROW>
    batch_outputs(i).dataset_id = cfg.dataset_id; %#ok<SAGROW>
    batch_outputs(i).event_diversity_file = W.save_file; %#ok<SAGROW>
    batch_outputs(i).top_csv_file = W.top_csv_file; %#ok<SAGROW>
end

local_emit(run_log_file, '\nFinished BLP event-diversity quick-check batch.\n');


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
