% Canonical interactive single-dataset BLP dictionary pipeline entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_preprocess_one_cfg_to_observables_streamed.m');
%
% Optional controls, set before run(...) if needed:
%   force_recompute = false;
%   save_precision = 'single';
%   chunk_size = 200000;
%   low_full_max_hz = 50;
%   high_max_hz = 250;
%   high_group_size = 2;
%   dict_modes = {'abs', 'complex_split'};
%   load_metadata_only = true;
%   cache_raw_to_disk = false;
%
% Role:
%   - canonical interactive entry for one manually defined cfg_*.m
%   - minimal top-to-bottom script that keeps intermediate variables in
%     workspace
%   - lighter than the batch cfg_names driver

%% -------------------------
%  Required user input
%  -------------------------
if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
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
cfg_name = char(string(cfg_name));

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

if ~exist('load_metadata_only', 'var') || isempty(load_metadata_only)
    load_metadata_only = true;
end

if ~exist('cache_raw_to_disk', 'var') || isempty(cache_raw_to_disk)
    cache_raw_to_disk = false;
end

%% -------------------------
%  Config and output root
%  -------------------------
cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

output_root = io_project.get_project_processed_root();

fprintf('Running streamed preprocessing for %s (%s)\n', cfg_name, cfg.dataset_id);
fprintf('Raw data root:\n  %s\n', cfg.raw_data_root);
fprintf('Output root:\n  %s\n', output_root);

%% -------------------------
%  Step 1. Load raw BLP
%  -------------------------
load_opts = struct();
load_opts.metadata_only = load_metadata_only;
load_opts.cache_to_disk = cache_raw_to_disk;
load_opts.force_recompute = force_recompute;

D = io_raw.load_blp_dataset(cfg, load_opts);
fprintf('Loaded concatenated raw data: [%d, %d]\n', D.n_time, numel(D.selected_channels));

%% -------------------------
%  Step 2. Streamed spectrograms
%  -------------------------
spec_opts = struct();
spec_opts.save_precision = save_precision;
spec_opts.return_data = false;
spec_opts.force_recompute = force_recompute;

S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, spec_opts);
fprintf('ABS spectrogram file:\n  %s\n', S.abs_file);
fprintf('COMPLEX spectrogram file:\n  %s\n', S.complex_file);
fprintf('ABS spectrogram size:\n  [%s]\n', num2str(S.tmpall_mean_abs_size));

%% -------------------------
%  Step 3. Dictionary modes
%  -------------------------
params = struct();
params.low_full_max_hz = low_full_max_hz;
params.high_max_hz = high_max_hz;
params.high_group_size = high_group_size;
params.chunk_size = chunk_size;
params.precision = save_precision;

dict_outputs = struct([]);
for i_mode = 1:numel(dict_modes)
    params.spec_mode = dict_modes{i_mode};
    dict = build_reskoopnet_dicts(D, cfg, output_root, params);

    dict_outputs(i_mode).spec_mode = params.spec_mode; %#ok<SAGROW>
    dict_outputs(i_mode).save_file = dict.save_file; %#ok<SAGROW>
    dict_outputs(i_mode).info_csv_file = dict.info_csv_file; %#ok<SAGROW>
    dict_outputs(i_mode).n_time = dict.n_time; %#ok<SAGROW>
    dict_outputs(i_mode).n_obs = dict.n_obs; %#ok<SAGROW>

    fprintf('%s observable file:\n  %s\n', upper(params.spec_mode), dict.save_file);
    fprintf('%s observable metadata CSV:\n  %s\n', upper(params.spec_mode), dict.info_csv_file);
end
