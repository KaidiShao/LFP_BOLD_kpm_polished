function cfg = cfg_E10bv1()
% Configuration for dataset E10.bv1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10bv1.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260414

%% basic info
cfg.dataset_id = 'E10.bv1';
cfg.raw_data_root = 'E:\DataPons\E10.bv1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10bv1';

%% channel selection
cfg.channels.sites = {'pl', 'hp', 'pl', 'pl', 'hp', 'hp', 'hp','pl','pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [2, 5:7];
cfg.channels.selected_site2 = [1, 3, 4, 8, 9];
cfg.channels.selected_all = [1:9];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:25;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

end
