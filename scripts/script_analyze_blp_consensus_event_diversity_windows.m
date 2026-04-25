% Dataset-specific wrapper currently fixed to f12m01.
% For the canonical e10gb1 diversity workflow, use:
%   script_analyze_blp_consensus_event_diversity_windows_e10gb1.m
%   script_refresh_e10gb1_diversity_top_windows_and_ref_figs.m

close all force;

cfg = cfg_F12m01();

output_root = get_project_processed_root();
consensus_input = [];
[C, source_consensus_file] = load_consensus_state_results(cfg, output_root, consensus_input);

params = struct();
params.window_length_samples = 5000;
params.keep_partial_window = false;
params.top_k = 10;
params.save_csv = true;
params.force_recompute = false;

W = analyze_blp_consensus_event_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    params, ...
    source_consensus_file);

disp(W.top_windows_table);
fprintf('Saved per-window event-diversity summary to:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('Saved all-window CSV to:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('Saved top-window CSV to:\n  %s\n', W.top_csv_file);
end
