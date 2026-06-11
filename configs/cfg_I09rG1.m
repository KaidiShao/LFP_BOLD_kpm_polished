function cfg = cfg_I09rG1()
% Configuration for dataset I09.rG1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\I09rG1.m
% NOTE: current G:\DataPons\I09.rG1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'I09.rG1';
cfg.raw_data_root = 'G:\DataPons\I09.rG1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'i09rg1';

%% channel selection
% Post-NOS-removal source site order.
cfg.channels.sites = {'cx', 'cx', 'cx', 'cx', 'vcx', 'cx', 'vcx', 'vcx', ...
                      'lgn', 'st', 'lgn', 'lgn', 'st', 'st', ...
                      'pl', 'pl', 'pl', 'pl', 'pl', 'pl', 'hp', 'hp', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = 21:23;
cfg.channels.selected_site2 = 15:20;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:3;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'vstim; flicker 32 Hz; excluded because no data files in current G root';

cfg.sessions(2).session_id = 4:9;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'spont; excluded because no data files in current G root';

cfg.sessions(3).session_id = 10:16;
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'spont1; excluded because no data files in current G root';

cfg.sessions(4).session_id = 17:18;
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = [];
cfg.sessions(4).notes = 'ephys'; % spontaneous without EPI/BOLD

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
