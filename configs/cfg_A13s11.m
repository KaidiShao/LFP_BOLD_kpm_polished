function cfg = cfg_A13s11()
% Configuration for dataset A13.s11
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\A13s11.m
% Kaidi Shao, 20260608

%% basic info
cfg.dataset_id = 'A13.s11';
cfg.raw_data_root = 'G:\DataPons\A13.s11\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'a13s11';

%% channel selection
% CHECKING STATE DEFINED BY XCOR OF WAVELETS sites from the source file.
cfg.channels.sites = {'vcx', 'vcx', 'vcx', 'vcx', 'vcx', 'cx', 'cx', 'cx', 'cx', ...
                      'lgn', 'lgn', 'st', 'st', 'lgn', 'lgn', 'st', ...
                      'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl', 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [17:20];
cfg.channels.selected_site2 = [21:24];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:3;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'vstim'; % flicker 32 Hz

cfg.sessions(2).session_id = 4:10;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(3).session_id = 11:18;
cfg.sessions(3).include = true;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'spont1'; % spontaneous continuation

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
