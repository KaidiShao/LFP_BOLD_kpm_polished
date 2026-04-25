% Canonical e10gb1 consensus-state-diversity analysis entry point.
% This script always recomputes the state-diversity ranking from the current
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
params.top_k = 30;
params.save_csv = true;
params.force_recompute = true;

W = analyze_blp_consensus_state_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    params, ...
    source_consensus_file);

disp(W.top_windows_table);
fprintf('Saved per-window consensus-state-diversity summary to:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('Saved all-window CSV to:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('Saved top-window CSV to:\n  %s\n', W.top_csv_file);
end
