this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

cfg = cfg_E10gb1();

% Add cfg.bold fields manually before running, for example:
% cfg.bold.data_subfolder = 'roits';    % default
% cfg.bold.input_var = 'roiTs';         % default
% cfg.bold.variable_mode = 'voxels';    % 'voxels' or 'roi_mean'
% cfg.bold.selected_region_names = {};  % exact roiTs{1,r}.name values; {} keeps all regions
% Exact semantic ROI names such as cfg.bold.role_map.hp / .elehp are best
% defined in the dataset cfg_*.m file rather than in scripts.

D = io_raw.load_bold_dataset(cfg);

X = D.data;
session_ids = D.session_ids;
session_lengths = D.session_lengths;
border_idx = D.border_idx;
selected_variables = D.selected_variables;
variable_labels = D.variable_labels;
