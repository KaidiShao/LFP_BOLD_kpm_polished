% Batch driver for the BLP consensus-state mainline pipeline.
%
% Usage from MATLAB:
%   cfg_names = {'E10fV1', 'E10gH1'};
%   run('scripts/script_run_cfgs_to_consensus_state_top_windows.m');
%
% Optional controls, set before run(...) if needed:
%   force_event_recompute = false;
%   force_density_recompute = false;
%   force_consensus_recompute = false;
%   force_summary_recompute = false;
%   force_window_recompute = false;
%   make_top_window_plots = true;
%   top_k = 30;
%   window_length_samples = 6000;
%   run_log_file = '';
%
% Role:
%   - batch driver for the canonical BLP consensus-state mainline
%   - accepts multiple manually defined cfg_*.m entries
%   - intended for running several datasets in one pass
%   - runs the canonical single-dataset script once per cfg_name
%
% For a lighter single-dataset interactive entry, use:
%   scripts/script_run_one_cfg_to_consensus_state_top_windows.m

%% -------------------------
%  Required user input
%  -------------------------
if ~exist('cfg_names', 'var') || isempty(cfg_names)
    error(['Set cfg_names before running this script, for example: ' ...
        'cfg_names = {''E10fV1'', ''E10gH1''};']);
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

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

%% -------------------------
%  User options
%  -------------------------
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

%% -------------------------
%  Output root and batch loop
%  -------------------------
output_root = io_project.get_project_processed_root();
batch_outputs = struct([]);
single_cfg_script = 'scripts/script_run_one_cfg_to_consensus_state_top_windows.m';

local_emit(run_log_file, 'Running BLP consensus-state batch pipeline for %d dataset(s)\n', numel(cfg_names));
local_emit(run_log_file, 'Output root:\n  %s\n', output_root);

for i = 1:numel(cfg_names)
    cfg_name = cfg_names{i};
    cfg_fun = ['cfg_' cfg_name];
    if exist(cfg_fun, 'file') ~= 2
        error('Config function not found on path: %s', cfg_fun);
    end

    cfg = feval(cfg_fun);
    cfg.spectrogram.pad_sec = 20;
    cfg.spectrogram.pad_mode = 'mirror';

    local_emit(run_log_file, '\n============================================================\n');
    local_emit(run_log_file, 'Running BLP consensus-state pipeline for %s (%s)\n', cfg_name, cfg.dataset_id);
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
    batch_outputs(i).state_diversity_file = W.save_file; %#ok<SAGROW>
    batch_outputs(i).top_csv_file = W.top_csv_file; %#ok<SAGROW>
    batch_outputs(i).plot_dir = ''; %#ok<SAGROW>
    if exist('P', 'var') && isstruct(P) && isfield(P, 'save_dir') && ~isempty(P.save_dir)
        batch_outputs(i).plot_dir = P.save_dir; %#ok<SAGROW>
    end
end

local_emit(run_log_file, '\nFinished BLP consensus-state batch pipeline.\n');


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
