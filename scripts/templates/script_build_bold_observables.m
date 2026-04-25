this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

cfg = cfg_E10gb1();

% Add cfg.bold fields manually before running.
D = load_bold_dataset(cfg);

pre = struct();
pre.demean = true;
pre.notch_hz = [];       % set if you want filtered BOLD first
P = preprocess_bold_sessions(D, pre);

params = struct();
params.source = 'data';  % 'data' or 'filtered'
params.mode = 'identity'; % 'identity', 'roi_mean', 'slow_band_power', or 'svd'
params.n_components = 50;
params.precision = 'single';

% Old BOLD power observable branch:
% params.mode = 'slow_band_power';
% params.include_raw = true;
% params.band_power_transform = 'power'; % 'power' matches the old scripts; 'log_power' is optional

O = build_bold_observables(P, params);

snap = make_session_aware_snapshot_pairs(O, 1);
