close all force;

cfg = cfg_F12m01();

output_root = 'D:\DataPons_processed\';
consensus_input = [];
[C, source_consensus_file] = load_consensus_state_results(cfg, output_root, consensus_input);

params = struct();
params.save_csv = true;
params.force_recompute = false;

S = summarize_blp_consensus_state_types( ...
    cfg, ...
    output_root, ...
    C, ...
    params, ...
    source_consensus_file);

disp(S.summary_table);
fprintf('Saved consensus-state type summary to:\n  %s\n', S.save_file);
if ~isempty(S.csv_file)
    fprintf('Saved CSV summary to:\n  %s\n', S.csv_file);
end
