function cfg = cfg_K13m20()
% Configuration for dataset K13.m20
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\K13m20.m
% Secondary reference: none found under D:\Onedrive\experimental_event_data\ripple_monkeys
% NOTE: this dataset currently only has a matching reference under forNikos.
% Kaidi Shao, 20260531

%% basic info
cfg.dataset_id = 'K13.m20';
cfg.raw_data_root = 'E:\DataPons\K13.m20\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'k13m20';

%% channel selection
cfg.channels.sites = {'hp'        'pl'  'pl'  'hp'  'hp'        'NOS'  'pl'  'pl' ...
                 'NOS'  'pbn'  'po'  'po'  'pbn'        'po'  'po'  'po'  'po' };
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [1, 4, 5];
cfg.channels.selected_site2 = [2, 3, 7, 8];
cfg.channels.selected_all = [1:5, 7:8];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 2:8;
cfg.sessions(1).include = true;
cfg.sessions(1).selected_channels = cfg.channels.selected_all;
cfg.sessions(1).notes = 'spont'; % spontaneous

cfg.sessions(2).session_id = 9:15;
cfg.sessions(2).include = true;
cfg.sessions(2).selected_channels = cfg.channels.selected_all;
cfg.sessions(2).notes = 'spont1'; % spontaneous continuation

cfg.sessions(3).session_id = [17, 19:20];
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = cfg.channels.selected_all;
cfg.sessions(3).notes = 'fspont'; % fast spontaneous, 4 slices, Tvol=0.5s

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
