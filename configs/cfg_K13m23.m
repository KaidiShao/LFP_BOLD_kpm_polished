function cfg = cfg_K13m23()
% Configuration for dataset K13.m23
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m23.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260519

%% basic info
cfg.dataset_id = 'K13.m23';
cfg.raw_data_root = 'E:\DataPons\K13.m23\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m23';

%% channel selection
cfg.channels.sites = { 'NOS'  'pbn'  'pbn'  'pbn'        'NOS'  'po'  'po'  'po'  'NOS' ...
                       'st'        'st'  'lgn'  'NOS'    'lgn' 'lgn'  'st'  'NOS' ...
                 'hp'  'pl'  'pl'  'pl'  'hp'      'pl'  'NOS'     'hp' };
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [18, 22,25];
cfg.channels.selected_site2 = [19,20,21,23];
cfg.channels.selected_all = [18,19,20,21,22,23,25];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 27:36;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

cfg.sessions(2).session_id = 43:56;
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
