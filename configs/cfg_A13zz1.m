function cfg = cfg_A13zz1()
% Configuration for dataset A13.zz1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\A13zz1.m
% NOTE: G:\DataPons\A13.zz1\roits currently contains no *_roits.mat files.
% Kaidi Shao, 20260608

%% basic info
cfg.dataset_id = 'A13.zz1';
cfg.raw_data_root = 'G:\DataPons\A13.zz1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'a13zz1';

%% channel selection
cfg.channels.sites = {'vcx', 'vcx', 'vcx', 'cx', 'vcx', 'vcx', 'vcx', 'cx', 'cx', ...
                      'st', 'lgn', 'lgn', 'lgn', 'st', 'lgn', 'lgn', 'lgn', 'st', ...
                      'hp', 'pl', 'pl', 'pl', 'pl', 'pl', 'pl', 'hp', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [19, 26, 27];
cfg.channels.selected_site2 = [20:25];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = [1, 3];
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'vstim'; % flicker 32 Hz; exp02 collapsed adfx

cfg.sessions(2).session_id = 4:10;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont'; % spontaneous with fMRI 10 min; BLP only in current G root

cfg.sessions(3).session_id = 21:30;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'spont2'; % source group, but matching BLP/roits files are absent under current G root

cfg.sessions(4).session_id = 31:40;
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = [];
cfg.sessions(4).notes = 'spont3'; % source group, but matching BLP/roits files are absent under current G root

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
