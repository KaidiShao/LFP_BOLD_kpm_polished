function cfg = cfg_CM029zz1()
% Configuration for dataset CM029.zz1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\CM029zz1.m
% NOTE: current G:\DataPons\CM029.zz1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'CM029.zz1';
cfg.raw_data_root = 'G:\DataPons\CM029.zz1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'cm029zz1';

%% channel selection
cfg.channels.sites = {'pl', 'pl', 'hp', 'hp', 'pl', 'pl', 'pl', 'pl', 'pl', ...
                      'pbl', 'pbl', 'pbl', 'pbl', 'po', 'po', 'po', 'po', 'NOS', 'NOS'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [3, 4];
cfg.channels.selected_site2 = [1, 2, 5:9];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = [2:4, 6:8, 10:12];
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'spont; exp01/05/09 MR-trigger shifted; excluded because no data files in current G root';

cfg.sessions(2).session_id = 13:22;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'spont1; excluded because no data files in current G root';

cfg.sessions(3).session_id = [23:30, 32, 33];
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'spont2; excluded because no data files in current G root';

cfg.sessions(4).session_id = 34:42;
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = [];
cfg.sessions(4).notes = 'spont3; excluded because no data files in current G root';

cfg.sessions(5).session_id = 43:51;
cfg.sessions(5).include = false;
cfg.sessions(5).selected_channels = [];
cfg.sessions(5).notes = 'spont4; excluded because no data files in current G root';

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
