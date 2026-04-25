cfg = cfg_F12m01();

output_root = get_project_processed_root();
event_input = [];
[R, source_event_file] = load_event_results(cfg, output_root, event_input);

params = struct();
params.bin_sec = 2;
params.smooth_sigma_sec = 2;
params.force_recompute = false;

E = compute_blp_event_density( ...
    cfg, ...
    output_root, ...
    R, ...
    params, ...
    source_event_file);
