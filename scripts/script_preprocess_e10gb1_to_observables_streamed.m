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

output_root = get_project_processed_root();

fprintf('Running streamed preprocessing for %s\n', cfg.dataset_id);
fprintf('Raw data root:\n  %s\n', cfg.raw_data_root);
fprintf('Output root:\n  %s\n', output_root);

D = load_blp_dataset(cfg);
fprintf('Loaded concatenated raw data: [%d, %d]\n', size(D.data, 1), size(D.data, 2));

spec_opts = struct();
spec_opts.save_precision = 'single';
spec_opts.return_data = false;
spec_opts.force_recompute = false;

S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, spec_opts);
fprintf('Streamed abs spectrogram file:\n  %s\n', S.abs_file);
fprintf('Streamed complex spectrogram file:\n  %s\n', S.complex_file);
fprintf('Saved/loaded region-mean spectrogram with abs size [%s]\n', ...
    num2str(S.tmpall_mean_abs_size));

params = struct();
params.low_full_max_hz = 50;
params.high_max_hz = 250;
params.high_group_size = 2;
params.chunk_size = 200000;
params.precision = 'single';

params.spec_mode = 'abs';
dict_abs = build_reskoopnet_dicts(D, cfg, output_root, params);

params.spec_mode = 'complex_split';
dict_complex = build_reskoopnet_dicts(D, cfg, output_root, params);

fprintf('ABS observable file:\n  %s\n', dict_abs.save_file);
fprintf('ABS observable metadata CSV:\n  %s\n', dict_abs.info_csv_file);
fprintf('COMPLEX_SPLIT observable file:\n  %s\n', dict_complex.save_file);
fprintf('COMPLEX_SPLIT observable metadata CSV:\n  %s\n', dict_complex.info_csv_file);
