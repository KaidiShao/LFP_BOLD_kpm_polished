function cfg = cfg_K13m21()
% Configuration for dataset K13.m21
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m21.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260531

%% basic info
cfg.dataset_id = 'K13.m21';
cfg.raw_data_root = 'E:\DataPons\K13.m21\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m21';

%% channel selection
cfg.channels.sites = {'pl'        'pl'  'pl'  'pl'  'hp'        'NOS'  'hp'  'hp' ...
                 'NOS'  'pbn'  'pbn'  'po'        'po'  'po'  'po'  'po'  'po' };
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [5, 7, 8];
cfg.channels.selected_site2 = [1:4];
cfg.channels.selected_all = [1:5, 7:8];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = [5, 8:12];
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

cfg.sessions(2).session_id = 13:15;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'fspont'; % fast spontaneous, 4 slices, Tvol=0.5s

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
