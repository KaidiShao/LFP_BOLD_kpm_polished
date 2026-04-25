close all force;

cfg = cfg_E10gb1();

output_root = get_project_processed_root();
event_input = [];
[R, source_event_file] = load_event_results(cfg, output_root, event_input);

params = struct();
params.min_channel_count = [];
params.require_region_presence = false;
params.required_regions = {'hp', 'pl'};
params.force_recompute = false;

C = compute_blp_consensus_states( ...
    cfg, ...
    output_root, ...
    R, ...
    params, ...
    source_event_file);

fprintf('Saved consensus-state result to:\n  %s\n', C.save_file);
