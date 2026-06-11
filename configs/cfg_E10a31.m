function cfg = cfg_E10a31()
% Configuration for dataset E10.a31
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10a31.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'E10.a31';
cfg.raw_data_root = 'G:\DataPons\E10.a31\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10a31';

%% channel selection
cfg.channels.sites = {'hp', 'pl', 'pl', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 4:7];
cfg.channels.selected_site2 = [2, 3, 8:10];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = [11, 13:14, 20, 21, 25];
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % exp12/22 BOLD transient in source notes

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
