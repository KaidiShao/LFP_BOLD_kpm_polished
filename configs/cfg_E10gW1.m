function cfg = cfg_E10gW1()
% Configuration for dataset E10.gW1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET.OLD\sesmonkeys\E10gW1.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: session inclusion currently follows the local pipeline choice for
%       this dataset; comments here are aligned to the E10.gW1 reference.
% Kaidi Shao, 20260428

%% basic info
cfg.dataset_id = 'E10.gW1';
cfg.raw_data_root = 'E:\DataPons\E10.gW1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'e10gw1';

%% channel selection
cfg.channels.sites = {'lgn' 'lgn' 'lgn' 'lgn' 'lgn' 'st' 'st' 'st' 'st'  'pl' 'pl' 'pl' 'hp' 'hp' 'hp' 'hp' 'hp'};
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

cfg.sessions(2).session_id = [6:16]; % validation-baseline spont subset; 
% more conservative ripple_monkeys subset = [6:7 10:11 20:24]
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont';

cfg.sessions(3).session_id = [21:25]; % validation-baseline spont subset; 
% more conservative ripple_monkeys subset = [6:7 10:11 20:24]
cfg.sessions(3).include = false; %% different voxels (eleHp)
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'spont1'; 


cfg.sessions(4).session_id = [31:35]; % available as vspont in the primary reference, but excluded here
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = cfg.channels.selected_all;
cfg.sessions(4).notes = 'vspont'; % spontaneous block collected during the visual-stim protocol

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
