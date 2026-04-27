function cfg = cfg_E10fV1()
% Configuration for dataset E10.fV1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10fV1.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260414

%% basic info
cfg.dataset_id = 'E10.fV1';
cfg.raw_data_root = 'E:\DataPons\E10.fV1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10fV1';

%% channel selection
cfg.channels.sites =  {'lgn' 'lgn' 'lgn' 'st' 'lgn' 'st' 'lgn' 'hp' 'hp' 'hp' 'hp' 'pl' 'hp' 'hp' 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [8:11, 13, 14];
cfg.channels.selected_site2 = [12, 15];
cfg.channels.selected_all = [8:15];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:5;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'polar'; % polar visual stimulation

cfg.sessions(2).session_id = 6:25;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont'; % spontaneous

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'eleHP';
cfg.bold.role_map.hp = 'HP';
end
