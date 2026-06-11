function cfg = cfg_K13m09()
% Configuration for dataset K13.m09
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m09.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'K13.m09';
cfg.raw_data_root = 'G:\DataPons\K13.m09\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m09';

%% channel selection
cfg.channels.sites = {'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', ...
                      'po', 'pbn', 'pbn', 'NOS', 'po', 'pbn', 'po'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = 1:4;
cfg.channels.selected_site2 = 5:7;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 2:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 11:20;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous with fMRI 10 min

cfg.sessions(3).session_id = 23;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'ephys'; % spontaneous without fMRI 5 min

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
