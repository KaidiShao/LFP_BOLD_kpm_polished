function cfg = cfg_F12m02()
% Configuration for dataset F12m02
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\F12m02.m
% Kaidi Shao, 20260605

%% basic info
cfg.dataset_id = 'F12.m02';
cfg.raw_data_root = 'E:\DataPons\F12.m02\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'f12m02';

%% channel selection
% AFTER CHECKFFT sites from the sesmonkeys source file.
cfg.channels.sites = {'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', 'NOS', 'hp', 'hp', ...
                      'NOS', 'po', 'po', 'po', 'pbn', 'po', 'po', 'po', 'pbn'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1:4, 9:10];
cfg.channels.selected_site2 = [5:7];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 11:14;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'fspont'; % fast spontaneous, 4 slices, Tvol=0.5s

cfg.sessions(3).session_id = 15:18;
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
