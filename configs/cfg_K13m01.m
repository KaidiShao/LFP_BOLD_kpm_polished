function cfg = cfg_K13m01()
% Configuration for dataset K13.m01
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m01.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'K13.m01';
cfg.raw_data_root = 'G:\DataPons\K13.m01\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m01';

%% channel selection
cfg.channels.sites = {'hp', 'hp', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', ...
                      'hp', 'pl', 'NOS', 'po', 'pbn', 'NOS', 'NOS'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1:6, 9];
cfg.channels.selected_site2 = [7:8, 10];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 3:8;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
