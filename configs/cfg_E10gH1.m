function cfg = cfg_E10gH1()
% Configuration for dataset E10.gH1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10gH1.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260414

%% basic info
cfg.dataset_id = 'E10.gH1';
cfg.raw_data_root = 'E:\DataPons\E10.gH1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10gh1';

%% channel selection
cfg.channels.sites =  {'lgn' 'lgn' 'lgn' 'lgn' 'lgn' 'st' 'st' 'st' 'st'  'pl' 'pl' 'pl' 'hp' 'hp'  ...
                 'hp' 'hp' 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [13:17];
cfg.channels.selected_site2 = [10:12];
cfg.channels.selected_all = [10:17];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:5;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'polar'; % polar visual stimulation

cfg.sessions(2).session_id = [6:16, 21:25];
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont'; % spontaneous
end
