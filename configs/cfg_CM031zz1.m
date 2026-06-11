function cfg = cfg_CM031zz1()
% Configuration for dataset CM031.zz1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\CM031zz1.m
% NOTE: current G:\DataPons\CM031.zz1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'CM031.zz1';
cfg.raw_data_root = 'G:\DataPons\CM031.zz1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'cm031zz1';

%% channel selection
cfg.channels.sites = {'po', 'pbn', 'pbn', 'pbn', 'pbn', 'pbn', 'pbn', 'pbn', 'po', 'po', ...
                      'pl', 'pl', 'pl', 'hp', 'hp', 'hp', 'hp', 'hp', 'hp', 'hp'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = 14:20;
cfg.channels.selected_site2 = 11:13;
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:12;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'spont; excluded because no data files in current G root';

cfg.sessions(2).session_id = [14:16, 18:22];
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'spont1; exp13/17 MR-trigger shifted; excluded because no data files in current G root';

cfg.sessions(3).session_id = [23:28, 30:32];
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'spont2; exp29 MR-trigger shifted; excluded because no data files in current G root';

cfg.sessions(4).session_id = [33:35, 37:42];
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = [];
cfg.sessions(4).notes = 'spont3; exp36 MR-trigger shifted; excluded because no data files in current G root';

cfg.sessions(5).session_id = [43:48, 50, 52];
cfg.sessions(5).include = false;
cfg.sessions(5).selected_channels = [];
cfg.sessions(5).notes = 'spont4; exp49/51 MR-trigger shifted; excluded because no data files in current G root';

cfg.sessions(6).session_id = 54:62;
cfg.sessions(6).include = false;
cfg.sessions(6).selected_channels = [];
cfg.sessions(6).notes = 'spont5; exp53 MR-trigger shifted; excluded because no data files in current G root';

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
