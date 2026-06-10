function cfg = cfg_A13qk1()
% Configuration for dataset A13.qk1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\A13qk1.m
% Kaidi Shao, 20260608

%% basic info
cfg.dataset_id = 'A13.qk1';
cfg.raw_data_root = 'G:\DataPons\A13.qk1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'a13qk1';

%% channel selection
% CHECKING STATE DEFINED BY XCOR OF WAVELETS sites from the source file.
cfg.channels.sites = {'cx', 'vcx', 'vcx', 'vcx', 'vcx', 'cx', 'cx', 'cx', 'cx', 'cx', ...
                      'st', 'lgn', 'lgn', 'lgn', 'st', 'st', 'st', 'st', 'st', ...
                      'hp', 'hp', 'pl', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [20, 21, 23:26];
cfg.channels.selected_site2 = [22, 27:29];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:8;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous with fMRI 10 min

cfg.sessions(2).session_id = 9:17;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous continuation

cfg.sessions(3).session_id = [18, 20];
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
