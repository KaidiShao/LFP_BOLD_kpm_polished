this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

cfg = cfg_E10gb1();

% Add cfg.bold fields manually before running.
D = load_bold_dataset(cfg);

params = struct();
params.demean = true;
params.zscore = false;
params.detrend = false;
params.notch_hz = [];       % set manually if a notch filter is needed
params.notch_q = 35;
params.bandpass_hz = [];    % set manually if bandpass-filtered BOLD is needed

P = preprocess_bold_sessions(D, params);
