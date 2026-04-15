function cfg = cfg_E10eA1()
% Configuration for dataset E10.eA1
% Primary reference: D:\Onedrive\experimental_event_data\ripple_monkeys\E10eA1.m
% Secondary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10eA1.m
% NOTE: the two reference files are not identical; this cfg uses the more
%       conservative Michel channel map and keeps the 1:10 / 11:30 session split.
% Kaidi Shao, 20260414

%% basic info
cfg.dataset_id = 'E10.eA1';
cfg.raw_data_root = 'E:\DataPons\E10.eA1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10ea1';

%% channel selection (SECONDARY)
% cfg.channels.sites = {'pl', 'hp', 'hp', 'pl', 'pl', 'hp', 'hp','hp','hp'};
% cfg.channels.selected_labels = {'hp', 'pl'};
% cfg.channels.selected_site1 = [2, 3, 6:8];
% cfg.channels.selected_site2 = [1, 4, 5];
% cfg.channels.selected_all = [1:9];

%% channel selection (MORE CONSERVATIVE, from Michel)
cfg.channels.sites = {'pl', 'nos', 'nos', 'sr', 'sr', 'sr', 'nos','nos','pl'};
cfg.channels.selected_labels = {'sr', 'pl'};
cfg.channels.selected_site1 = [4:6];
cfg.channels.selected_site2 = [1, 9];
cfg.channels.selected_all = [1, 4:6, 9];


%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:10;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

cfg.sessions(2).session_id = 11:30; % MORE CONSERVATIVE FROM DESNET
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'vstim'; % visual stimulation
end
