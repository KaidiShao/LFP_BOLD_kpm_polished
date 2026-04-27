this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

cfg = cfg_E10gb1();

% Add cfg.bold fields manually before running.
% cfg.bold.data_subfolder = 'roits';    % default
% cfg.bold.input_var = 'roiTs';         % default
% cfg.bold.variable_mode = 'voxels';    % 'voxels' or 'roi_mean'
% cfg.bold.selected_region_names = {};  % exact roiTs{1,r}.name values; {} keeps all regions
% Exact semantic ROI names such as cfg.bold.role_map.hp / .elehp are best
% defined in the dataset cfg_*.m file rather than in scripts.
D = load_bold_dataset(cfg);

params = struct();
params.demean = true;
params.zscore = false;
params.detrend = false;
params.notch_hz = [];       % set manually if a notch filter is needed
params.notch_q = 35;
params.bandpass_hz = [];    % set manually if bandpass-filtered BOLD is needed

P = preprocess_bold_sessions(D, params);
