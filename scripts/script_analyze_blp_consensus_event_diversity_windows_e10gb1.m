% Canonical e10gb1 diversity-analysis entry point.
% This script always recomputes the diversity ranking from the current
% consensus-state result using:
%   - window_length_samples = 6000
%   - window_mode = 'global'

close all force;

cfg = cfg_E10gb1();

output_root = get_project_processed_root();
consensus_input = [];
[C, source_consensus_file] = load_consensus_state_results(cfg, output_root, consensus_input);

params = struct();
params.window_length_samples = 6000;
params.window_mode = 'global';
params.keep_partial_window = false;
params.top_k = 10;
params.save_csv = true;
params.force_recompute = true;

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
