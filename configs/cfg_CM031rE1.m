function cfg = cfg_CM031rE1()
% Configuration for dataset CM031.rE1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\CM031rE1.m
% NOTE: current G:\DataPons\CM031.rE1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'CM031.rE1';
cfg.raw_data_root = 'G:\DataPons\CM031.rE1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'cm031re1';

%% channel selection
% Source sites are LGN/ST-PONS, not HP/PL. The current mainline keeps HP/PL
% channels, so no session is marked runnable here.
cfg.channels.sites = {'st', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'NOS', ...
                      'po', 'pbn', 'po', 'NOS', 'po', 'NOS', 'NOS', 'NOS'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [];
cfg.channels.selected_site2 = [];
cfg.channels.selected_all = [];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 5:10;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'spont; excluded because no HP/PL channels and no data files in current G root';

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
