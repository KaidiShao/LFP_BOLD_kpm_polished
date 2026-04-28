addpath(genpath('D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished\'))

cfg = cfg_F12m01();
D = io_raw.load_blp_dataset(cfg);

X = D.data;
session_ids = D.session_ids;
session_lengths = D.session_lengths;
border_idx = D.border_idx;
selected_channels = D.selected_channels;