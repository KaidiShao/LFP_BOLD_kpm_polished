% Fixed-dataset example for state-diversity plotting.
cfg_name = 'E10gb1';
window_length_samples = 6000;
window_mode = 'global';
skip_existing = false;
run('scripts/script_plot_top_consensus_state_diversity_windows.m');
