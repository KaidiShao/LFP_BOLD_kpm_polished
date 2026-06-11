function cfg = cfg_K13m04()
% Configuration for dataset K13.m04
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m04.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'K13.m04';
cfg.raw_data_root = 'G:\DataPons\K13.m04\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m04';

%% channel selection
cfg.channels.sites = {'NOS', 'hp', 'hp', 'hp', 'pl', 'pl', 'hp', ...
                      'pbn', 'pbn', 'po', 'po', 'pbn', 'pbn', 'pbn', 'NOS', ...
                      'pbn', 'po'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [2:4, 7];
cfg.channels.selected_site2 = 5:6;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 2:9;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 10:15;
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
