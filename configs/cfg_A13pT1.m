function cfg = cfg_A13pT1()
% Configuration for dataset A13.pT1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\A13pT1.m
% Kaidi Shao, 20260608

%% basic info
cfg.dataset_id = 'A13.pT1';
cfg.raw_data_root = 'G:\DataPons\A13.pT1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'a13pt1';

%% channel selection
% CHECKING STATE DEFINED BY XCOR OF WAVELETS sites from the source file.
cfg.channels.sites = {'vcx', 'vcx', 'vcx', 'cx', 'cx', 'cx', 'cx', 'cx', ...
                      'st', 'st', 'st', 'lgn', 'st', 'lgn', 'st', 'lgn', 'st', ...
                      'hp', 'pl', 'hp', 'pl', 'pl', 'hp', 'pl', 'pl', 'pl', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [18, 20, 23, 27];
cfg.channels.selected_site2 = [19, 21, 22, 24:26];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = [1, 3, 5, 6];
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 16:18;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'ephys'; % spontaneous without EPI/BOLD

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
