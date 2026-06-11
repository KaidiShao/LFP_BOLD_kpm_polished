function cfg = cfg_I11bu1()
% Configuration for dataset I11.bu1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\I11bu1.m
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'I11.bu1';
cfg.raw_data_root = 'G:\DataPons\I11.bu1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'i11bu1';

%% channel selection
cfg.channels.sites = {'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1:4, 8];
cfg.channels.selected_site2 = 5:7;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 8:15;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 16:20;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous continuation

cfg.sessions(3).session_id = 21:30;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'estim'; % stimulation; matching BLP files absent under current G root

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
