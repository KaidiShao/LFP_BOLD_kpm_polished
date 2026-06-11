function cfg = cfg_E10ed1()
% Configuration for dataset E10.ed1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10ed1.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'E10.ed1';
cfg.raw_data_root = 'G:\DataPons\E10.ed1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10ed1';

%% channel selection
cfg.channels.sites = {'pl', 'hp', 'pl', 'pl', 'pl', 'hp', 'hp', 'hp', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [2, 6:9];
cfg.channels.selected_site2 = [1, 3:5];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:8;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 11:25;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'estim'; % DES-Hippocampus, 200 us, 100 Hz, 100 ms burst, 8 uA

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
