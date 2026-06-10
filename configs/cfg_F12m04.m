function cfg = cfg_F12m04()
% Configuration for dataset F12m04
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\F12m04.m
% Kaidi Shao, 20260605

%% basic info
cfg.dataset_id = 'F12.m04';
cfg.raw_data_root = 'E:\DataPons\F12.m04\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'f12m04';

%% channel selection
% AFTER CHECKFFT sites from the sesmonkeys source file.
% This source file is LGN/ST-PONS, not HP/PL. The current mainline keeps only
% HP/PL channels, so no session is marked runnable here.
cfg.channels.sites = {'lgn', 'st', 'lgn', 'lgn', 'NOS', 'st', 'st', 'NOS', ...
                      'NOS', 'po', 'po', 'po', 'po', 'pbn', 'po'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [];
cfg.channels.selected_site2 = [];
cfg.channels.selected_all = [];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 7:13;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'spont; excluded because no HP/PL channels in AFTER CHECKFFT sites';

cfg.sessions(2).session_id = 1:4;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'vstim; LED flicker, not spontaneous mainline';

cfg.sessions(3).session_id = 21:30;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'ephys'; % spontaneous without MRI

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'eleHP';
cfg.bold.role_map.hp = 'HP';

end
