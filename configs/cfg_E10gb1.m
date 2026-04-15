function cfg = cfg_E10gb1()
% Configuration for dataset E10.gb1
% Role: validation baseline for reproducing the previously used pipeline input
% Primary reference: D:\Onedrive\experimental_event_data\ripple_monkeys\E10.gb1\E10gb1.m
% Secondary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\E10gb1.m
% NOTE: this cfg intentionally keeps the broader manually annotated session
%       subset used in prior runs, even though the reference files disagree
%       on the more conservative spont-session exclusions.
% Kaidi Shao, 20260414

%% basic info
cfg.dataset_id = 'E10.gb1';
cfg.raw_data_root = 'E:\DataPons\E10.gb1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10gb1';

%% channel selection
cfg.channels.sites = {'lgn' 'lgn' 'st' 'lgn' 'lgn' 'lgn' 'st' 'st' 'hp' 'hp' 'hp' 'hp' 'pl' 'pl' 'pl' 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [9:12];
cfg.channels.selected_site2 = [13:16];
cfg.channels.selected_all = [9:16];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:5;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'polar'; % polar visual stimulation

cfg.sessions(2).session_id = [6:7 9:13 20:25]; % validation-baseline spont subset; 
% more conservative ripple_monkeys subset = [6:7 10:11 20:24]
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont';

cfg.sessions(3).session_id = [31:35]; % available as vspont in the primary reference, but excluded here
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'vspont'; % spontaneous block collected during the visual-stim protocol

end
