% Batch BLP dictionary pipeline driver for multiple cfg_*.m definitions.
%
% Usage from MATLAB:
%   cfg_names = {'E10fV1', 'E10gH1'};
%   run('scripts/script_preprocess_cfgs_to_observables_streamed.m');
%
% Optional controls, set before run(...) if needed:
%   force_recompute = false;
%   save_precision = 'single';
%   chunk_size = 200000;
%   low_full_max_hz = 50;
%   high_max_hz = 250;
%   high_group_size = 2;
%   dict_modes = {'abs', 'complex_split'};
%   cache_raw_to_disk = false;
%   load_metadata_only = true;
%   run_log_file = '';
%
% Role:
%   - batch driver for the BLP dictionary pipeline
%   - accepts any manually defined cfg_*.m that matches the current BLP
%     loader contract
%   - intended for running multiple datasets in one pass
%
% For a lighter single-dataset interactive entry, use:
%   scripts/script_preprocess_one_cfg_to_observables_streamed.m

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

eeglab_root = 'D:\Onedrive\Toolbox\eeglab10_2_5_8b\';
if exist(eeglab_root, 'dir') == 7 && exist('finputcheck', 'file') ~= 2
    addpath(genpath(eeglab_root));
end

if exist('timefreqMB', 'file') ~= 2
    error('timefreqMB was not found on the MATLAB path.');
end

if exist('finputcheck', 'file') ~= 2
    error(['finputcheck was not found on the MATLAB path, and the expected ' ...
        'EEGLAB directory does not provide it: %s'], eeglab_root);
end

%% -------------------------
%  User options
%  -------------------------
cfg_names = cellstr(string(cfg_names));

if ~exist('force_recompute', 'var') || isempty(force_recompute)
    force_recompute = false;
end

if ~exist('save_precision', 'var') || isempty(save_precision)
    save_precision = 'single';
end

if ~exist('chunk_size', 'var') || isempty(chunk_size)
    chunk_size = 200000;
end

if ~exist('low_full_max_hz', 'var') || isempty(low_full_max_hz)
    low_full_max_hz = 50;
end

if ~exist('high_max_hz', 'var') || isempty(high_max_hz)
    high_max_hz = 250;
end

if ~exist('high_group_size', 'var') || isempty(high_group_size)
    high_group_size = 2;
end

if ~exist('dict_modes', 'var') || isempty(dict_modes)
    dict_modes = {'abs', 'complex_split'};
end
dict_modes = cellstr(string(dict_modes));

if ~exist('cache_raw_to_disk', 'var') || isempty(cache_raw_to_disk)
    cache_raw_to_disk = false;
end

if ~exist('load_metadata_only', 'var') || isempty(load_metadata_only)
    load_metadata_only = true;
end

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
%  Output root and run summary
%  -------------------------
output_root = io_project.get_project_processed_root();

local_emit(run_log_file, 'Processed-data root:\n  %s\n', output_root);
local_emit(run_log_file, 'Datasets to preprocess:\n');
for i = 1:numel(cfg_names)
    local_emit(run_log_file, '  - %s\n', cfg_names{i});
end

%% -------------------------
%  Dataset loop
%  -------------------------
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
    local_emit(run_log_file, 'Running streamed preprocessing for %s (%s)\n', cfg_name, cfg.dataset_id);
    local_emit(run_log_file, 'Raw data root:\n  %s\n', cfg.raw_data_root);

    % Step 1. Load raw BLP
    load_opts = struct();
    load_opts.metadata_only = load_metadata_only;
    load_opts.cache_to_disk = cache_raw_to_disk;
    load_opts.force_recompute = force_recompute;

    load_log = evalc('D = io_raw.load_blp_dataset(cfg, load_opts);');
    local_append_text(run_log_file, load_log);
    local_emit(run_log_file, 'Loaded concatenated raw data: [%d, %d]\n', D.n_time, numel(D.selected_channels));

    % Step 2. Streamed spectrograms
    spec_opts = struct();
    spec_opts.save_precision = save_precision;
    spec_opts.return_data = false;
    spec_opts.force_recompute = force_recompute;

    spec_log = evalc('S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, spec_opts);');
    local_append_text(run_log_file, spec_log);
    local_emit(run_log_file, 'ABS spectrogram file:\n  %s\n', S.abs_file);
    local_emit(run_log_file, 'COMPLEX spectrogram file:\n  %s\n', S.complex_file);
    local_emit(run_log_file, 'ABS spectrogram size:\n  [%s]\n', num2str(S.tmpall_mean_abs_size));

    % Step 3. Dictionary modes
    params = struct();
    params.low_full_max_hz = low_full_max_hz;
    params.high_max_hz = high_max_hz;
    params.high_group_size = high_group_size;
    params.chunk_size = chunk_size;
    params.precision = save_precision;

    for j = 1:numel(dict_modes)
        params.spec_mode = dict_modes{j};
        dict_log = evalc('dict = build_reskoopnet_dicts(D, cfg, output_root, params);');
        local_append_text(run_log_file, dict_log);

        local_emit(run_log_file, '%s observable file:\n  %s\n', upper(params.spec_mode), dict.save_file);
        local_emit(run_log_file, '%s observable metadata CSV:\n  %s\n', upper(params.spec_mode), dict.info_csv_file);
    end
end

local_emit(run_log_file, '\nFinished preprocessing %d dataset(s).\n', numel(cfg_names));

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
