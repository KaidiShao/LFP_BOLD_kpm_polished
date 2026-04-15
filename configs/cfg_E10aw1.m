function cfg = cfg_E10aw1()
% Configuration for dataset E10.aW1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10aW1.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: channel labels follow the updated manual labeling in the primary reference.
% Kaidi Shao, 20260414

%% basic info
cfg.dataset_id = 'E10.aW1';
cfg.raw_data_root = 'E:\DataPons\E10.aW1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10aw1';

%% channel selection
cfg.channels.sites = {'pl', 'hp', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [2:6];
cfg.channels.selected_site2 = [1, 7:10];
cfg.channels.selected_all = [1:10];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:25;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

end
