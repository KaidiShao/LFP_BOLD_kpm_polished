function cfg = cfg_K13m11()
% Configuration for dataset K13.m11
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m11.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'K13.m11';
cfg.raw_data_root = 'G:\DataPons\K13.m11\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m11';

%% channel selection
cfg.channels.sites = {'hp', 'pl', 'hp', 'hp', 'hp', 'pl', 'pbn', ...
                      'pbn', 'po', 'NOS', 'NOS', 'NOS', 'NOS'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 3:5];
cfg.channels.selected_site2 = [2, 6];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = [2:3, 5:12];
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 13:20;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous with fMRI 10 min

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
