function cfg = cfg_D11oE1()
% Configuration for dataset D11.oE1
% Primary reference: D:\Onedrive\ICPBR\forNikos\MatLab\DESNET\sesmonkeys\D11oE1.m
% NOTE: current G:\DataPons\D11.oE1 contains no data files.
% Kaidi Shao, 20260611

%% basic info
cfg.dataset_id = 'D11.oE1';
cfg.raw_data_root = 'G:\DataPons\D11.oE1\';
cfg.data_subfolder = 'blp';
cfg.file_stem = 'd11oe1';

%% channel selection
cfg.channels.sites = {'vcx', 'vcx', 'vcx', 'vcx', 'vcx', 'vcx', 'vcx', 'vcx', 'vcx', 'vcx', ...
                      'st', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'lgn', 'st', ...
                      'pl', 'pl', 'hp', 'hp', 'hp', 'hp', 'hp', 'pl', 'pl', 'pl'};
cfg.channels.selected_labels = {'hp', 'pl'};
cfg.channels.selected_site1 = 23:27;
cfg.channels.selected_site2 = [21, 22, 28:30];
cfg.channels.selected_all = [cfg.channels.selected_site1, cfg.channels.selected_site2];

%% session selection
cfg.sessions = struct([]);

cfg.sessions(1).session_id = 1:10;
cfg.sessions(1).include = false;
cfg.sessions(1).selected_channels = [];
cfg.sessions(1).notes = 'spont; excluded because no data files in current G root';

%% BOLD ROI naming
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'ele_HP';
cfg.bold.role_map.hp = 'HP';

end
