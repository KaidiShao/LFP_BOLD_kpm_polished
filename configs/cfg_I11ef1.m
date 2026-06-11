function cfg = cfg_I11ef1()
% Configuration for dataset I11.ef1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\I11ef1.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'I11.ef1';
cfg.raw_data_root = 'G:\DataPons\I11.ef1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'i11ef1';

%% channel selection
cfg.channels.sites = {'hp', 'hp', 'pl', 'pl', 'pl', 'pl', 'pl', 'hp', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 2, 8, 9];
cfg.channels.selected_site2 = 3:7;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 11:30;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'estim'; % DES-Hippocampus, 200 us, 100 Hz, 100 ms burst, 8 uA

cfg.sessions(3).session_id = 31:40;
cfg.sessions(3).include = true;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'spont2'; % spontaneous after DES

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
