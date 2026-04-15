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

cfg = cfg_E10gb1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

output_root = 'E:\DataPons_processed\';

fprintf('Running preprocessing for %s\n', cfg.dataset_id);
fprintf('Raw data root:\n  %s\n', cfg.raw_data_root);
fprintf('Output root:\n  %s\n', output_root);

D = load_blp_dataset(cfg);
fprintf('Loaded concatenated raw data: [%d, %d]\n', size(D.data, 1), size(D.data, 2));

S = compute_blp_region_spectrograms(D, cfg, output_root);
fprintf('Saved/loaded region-mean spectrogram with abs size [%d, %d, %d]\n', ...
    size(S.tmpall_mean_abs, 1), size(S.tmpall_mean_abs, 2), size(S.tmpall_mean_abs, 3));

params = struct();
params.spec_mode = 'complex_split';
params.low_full_max_hz = 50;
params.high_max_hz = 250;
params.high_group_size = 2;
params.chunk_size = 200000;
params.precision = 'single';

dict = build_reskoopnet_dicts(D, cfg, output_root, params);

fprintf('Observable file:\n  %s\n', dict.save_file);
fprintf('Observable metadata CSV:\n  %s\n', dict.info_csv_file);
