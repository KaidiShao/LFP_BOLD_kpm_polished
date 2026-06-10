function cfg = cfg_F12m03()
% Configuration for dataset F12m03
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\F12m03.m
% Kaidi Shao, 20260605

%% basic info
cfg.dataset_id = 'F12.m03';
cfg.raw_data_root = 'E:\DataPons\F12.m03\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'f12m03';

%% channel selection
% AFTER CHECKFFT sites from the sesmonkeys source file.
cfg.channels.sites = {'hp', 'pl', 'pl', 'pl', 'hp', 'hp', 'hp', 'NOS', 'hp', 'hp', ...
                      'NOS', 'pbn', 'po', 'po', 'NOS', 'po', 'po', 'po', 'po'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 5:7, 9:10];
cfg.channels.selected_site2 = [2:4];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 2:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min; exp01 shifted trigger

cfg.sessions(2).session_id = [12, 14];
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'fspont'; % fast spontaneous, 4 slices, Tvol=0.5s; exp11/13 shifted trigger

cfg.sessions(3).session_id = 15:19;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'ephys'; % spontaneous without MRI

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'eleHP';
cfg.bold.role_map.hp = 'HP';

end
