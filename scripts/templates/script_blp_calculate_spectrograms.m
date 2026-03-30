addpath(genpath('D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished\'))
addpath(genpath('D:\Onedrive\Toolbox\eeglab10_2_5_8b\'));

cfg = cfg_F12m01();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

D = load_blp_dataset(cfg);
output_root =  'D:\DataPons_processed\';
S = compute_blp_region_spectrograms(D, cfg, output_root);