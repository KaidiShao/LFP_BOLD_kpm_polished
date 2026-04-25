% Canonical end-to-end refresh for e10gb1 consensus-state-diversity windows.
% Order:
%   1. recompute 6000-sample global-window state-diversity ranking
%   2. redraw top-window trace/spectrogram figures
%   3. export matching reference Koopman figures

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

run(fullfile(this_script_dir, 'script_analyze_blp_consensus_state_diversity_windows_e10gb1.m'));
run(fullfile(this_script_dir, 'script_plot_top_consensus_state_diversity_windows_e10gb1.m'));
run(fullfile(this_script_dir, 'script_export_e10gb1_state_diversity_ref_figs_6000.m'));
