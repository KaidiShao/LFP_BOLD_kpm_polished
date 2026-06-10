function cfg = cfg_A13pF1()
% Configuration for dataset A13.pF1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\A13pF1.m
% Kaidi Shao, 20260608

%% basic info
cfg.dataset_id = 'A13.pF1';
cfg.raw_data_root = 'G:\DataPons\A13.pF1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'a13pf1';

%% channel selection
% The source removes one NOS channel before BLP export; sites below match the
% 19-channel BLP files under G:\DataPons\A13.pF1\blp.
cfg.channels.sites = {'hp', 'hp', 'hp', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', ...
                      'lgn', 'lgn', 'lgn', 'st', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1:7];
cfg.channels.selected_site2 = [8:9];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:6;
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
