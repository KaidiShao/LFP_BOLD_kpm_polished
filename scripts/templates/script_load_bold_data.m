this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

cfg = cfg_E10gb1();

% Add cfg.bold fields manually before running, for example:
% cfg.bold.input_file = 'D:\path\to\bold_roiTs.mat';
% cfg.bold.input_var = 'roiTs';
% cfg.bold.dx = 2;                   % TR in seconds
% cfg.bold.variable_mode = 'voxels'; % 'voxels', 'roi_mean', or 'observables'
% cfg.bold.selected_regions = [];    % leave empty for all regions
% cfg.bold.region_labels = {};
% Optional fallback only if roiTs has no session-length metadata:
% cfg.bold.session_lengths = [];

D = load_bold_dataset(cfg);

X = D.data;
session_ids = D.session_ids;
session_lengths = D.session_lengths;
border_idx = D.border_idx;
selected_variables = D.selected_variables;
variable_labels = D.variable_labels;
