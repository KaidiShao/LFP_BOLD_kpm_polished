function cfg = cfg_K13m16()
% Configuration for dataset K13.m16
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m16.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'K13.m16';
cfg.raw_data_root = 'G:\DataPons\K13.m16\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m16';

%% channel selection
cfg.channels.sites = {'pl', 'NOS', 'pl', 'hp', 'hp', 'hp', 'hp', 'NOS', ...
                      'NOS', 'NOS', 'pbn', 'pbn', 'po', 'po', 'NOS', 'po'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = 4:7;
cfg.channels.selected_site2 = [1, 3];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = [2, 4:8];
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 9:16;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous with fMRI 10 min

cfg.sessions(3).session_id = 18:20;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'fspont'; % fast spontaneous, 4 slices, Tvol=0.5s

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
