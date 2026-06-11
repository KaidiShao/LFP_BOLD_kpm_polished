function cfg = cfg_CM029pQ1()
% Configuration for dataset CM029.pQ1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\CM029pQ1.m
% NOTE: current G:\DataPons\CM029.pQ1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'CM029.pQ1';
cfg.raw_data_root = 'G:\DataPons\CM029.pQ1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'cm029pq1';

%% channel selection
% Source sites are LGN/ST-PONS, not HP/PL. The current mainline keeps HP/PL
% channels, so no session is marked runnable here.
cfg.channels.sites = {'lgn', 'NOS', 'st', 'lgn', 'st', 'lgn', 'NOS', 'lgn', 'lgn', ...
                      'po', 'pbn', 'NOS', 'po', 'pbn', 'pbn', 'NOS'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [];
cfg.channels.selected_site2 = [];
cfg.channels.selected_all = [];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:6;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'spont; excluded because no HP/PL channels and no data files in current G root';

cfg.sessions(2).session_id = 7:9;
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
