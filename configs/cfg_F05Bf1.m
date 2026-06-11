function cfg = cfg_F05Bf1()
% Configuration for dataset F05.Bf1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\f05bf1.m
% NOTE: current G:\DataPons\F05.Bf1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'F05.Bf1';
cfg.raw_data_root = 'G:\DataPons\F05.Bf1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'f05bf1';

%% channel selection
% This is an LGN microstimulation dataset. The source uses one physiology
% channel for visual/electrical stimulation timing rather than an HP/PL map.
cfg.channels.sites = {'lgn'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = [];
cfg.channels.selected_site2 = [];
cfg.channels.selected_all = [];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:3;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'estim; biphasic microstimulation aggregate; excluded because no data files in current G root';

cfg.sessions(2).session_id = 1;
cfg.sessions(2).include = false;
cfg.sessions(2).selected_channels = [];
cfg.sessions(2).notes = 'estim1; biphasic 100 Hz / 1000 uA; excluded because no data files in current G root';

cfg.sessions(3).session_id = [2, 3];
cfg.sessions(3).include = false;
cfg.sessions(3).selected_channels = [];
cfg.sessions(3).notes = 'estim2; biphasic 100 Hz / 750 uA; excluded because no data files in current G root';

cfg.sessions(4).session_id = 4:8;
cfg.sessions(4).include = false;
cfg.sessions(4).selected_channels = [];
cfg.sessions(4).notes = 'estim3; sinc 500 uA, 100/18/14/10 Hz; excluded because no data files in current G root';

cfg.sessions(5).session_id = 4:8;
cfg.sessions(5).include = false;
cfg.sessions(5).selected_channels = [];
cfg.sessions(5).notes = 'estim3o; sinc 500 uA ROI-only variant; excluded because no data files in current G root';

cfg.sessions(6).session_id = 9:13;
cfg.sessions(6).include = false;
cfg.sessions(6).selected_channels = [];
cfg.sessions(6).notes = 'estim4; spike 500 uA, 100/18/14/10 Hz; excluded because no data files in current G root';

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
