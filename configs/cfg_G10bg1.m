function cfg = cfg_G10bg1()
% Configuration for dataset G10.bg1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\G10bg1.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'G10.bg1';
cfg.raw_data_root = 'G:\DataPons\G10.bg1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'g10bg1';

%% channel selection
cfg.channels.sites = {'hp', 'pl', 'pl', 'pl', 'pl', 'hp', 'hp', 'hp', 'hp', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 6:10];
cfg.channels.selected_site2 = 2:5;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 11:20;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous continuation

cfg.sessions(3).session_id = 21:25;
cfg.sessions(3).include = true;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'spont2'; % spontaneous continuation

cfg.sessions(4).session_id = 26:35;
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = [];
cfg.sessions(4).notes = 'estim'; % stimulation; matching BLP files absent under current G root

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
