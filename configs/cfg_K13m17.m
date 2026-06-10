function cfg = cfg_K13m17()
% Configuration for dataset K13.m17
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m17.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260519

%% basic info
cfg.dataset_id = 'K13.m17';
cfg.raw_data_root = 'E:\DataPons\K13.m17\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m17';

%% channel selection
cfg.channels.sites = { 'hp'  'hp'  'pl'  'hp'        'hp'  'pl'  'pl'  'hp' ...
                 'po'  'po'  'pbn'  'po'   'pbn'  'po'  'po'  'po'  'po'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 2, 4, 5, 8];
cfg.channels.selected_site2 = [3,6,7];
cfg.channels.selected_all = [1:8];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 2:18;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

cfg.sessions(2).session_id = 20:22;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'fspont'; % spontaneous

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
