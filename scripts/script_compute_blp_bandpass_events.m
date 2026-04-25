cfg = cfg_F12m01();

D = load_blp_dataset(cfg);
output_root = get_project_processed_root();

params = struct();
params.passband = [2, 15; 30, 90; 90, 190];
params.band_labels = {'theta', 'gamma', 'ripple'};
params.L_start_range = [151, 101, 51];
params.L_extract_range = [301, 201, 101];
params.ThresRatio_range = [3.5, 4, 4];
params.input_normalization = 'zscore_per_channel';
params.force_recompute = false;

R = compute_blp_bandpass_events( ...
    D, cfg, output_root, params);
