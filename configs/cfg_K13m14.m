function cfg = cfg_K13m14()
% Configuration for dataset K13.m14
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m14.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'K13.m14';
cfg.raw_data_root = 'G:\DataPons\K13.m14\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m14';

%% channel selection
cfg.channels.sites = {'hp', 'pl', 'pl', 'hp', 'hp', 'pl', 'hp', ...
                      'po', 'pbn', 'pbn', 'pbn', 'pbn', 'NOS', 'NOS', 'NOS', 'NOS'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 4:5, 7];
cfg.channels.selected_site2 = [2:3, 6];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 2:14;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 17:18;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'fspont'; % fast spontaneous, 4 slices, Tvol=0.5s

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
