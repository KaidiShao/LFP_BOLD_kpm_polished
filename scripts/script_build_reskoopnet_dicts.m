cfg = cfg_F12m01();

D = load_blp_dataset(cfg);

params = struct();
params.spec_mode = 'complex_split';   % 'abs' or 'complex_split'
params.low_full_max_hz = 50;
params.high_max_hz = 250;
params.high_group_size = 2;
params.chunk_size = 200000;
params.precision = 'single';

output_root = get_project_processed_root();

dict = build_reskoopnet_dicts( ...
    D, cfg, output_root, params);
