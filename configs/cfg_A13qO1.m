function cfg = cfg_A13qO1()
% Configuration for dataset A13.qO1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\A13qO1.m
% Kaidi Shao, 20260608

%% basic info
cfg.dataset_id = 'A13.qO1';
cfg.raw_data_root = 'G:\DataPons\A13.qO1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'a13qo1';

%% channel selection
% The source lists 27 ADF channels, but the fMRI spontaneous BLP exports
% (sessions 3:14) already omit the NOS/unusable channels. Sites below match
% the 22-column BLP matrices for those sessions.
cfg.channels.sites = {'cx', 'cx', 'cx', 'vcx', 'cx', 'cx', 'cx', 'vcx', 'cx', ...
                      'lgn', 'st', 'lgn', 'lgn', 'st', 'st', ...
                      'pl', 'pl', 'hp', 'hp', 'pl', 'pl', 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [18, 19];
cfg.channels.selected_site2 = [16, 17, 20:22];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:2;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'vstim'; % flicker 32 Hz; 27-column BLP layout, not mixed with 22-column fMRI block

cfg.sessions(2).session_id = 3:14;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(3).session_id = 15:16;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'ephys'; % spontaneous without EPI/BOLD

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
