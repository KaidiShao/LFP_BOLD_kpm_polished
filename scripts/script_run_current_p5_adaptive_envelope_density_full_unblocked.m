% Batch wrapper for current P5 adaptive RMS-envelope density products.
%
% This script intentionally writes to *_rmsenv_adaptive condition tags so the
% existing P5 abs-magnitude outputs remain intact.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

dataset_filter = {};
condition_filter = {};
method_filter = {'svd', 'nmf', 'mds', 'umap'};
component_count_sweep = 3:8;

activity_suffix = 'rmsenv_adaptive';
density_value_transform = 'abs';
lfp_activity_transform = 'rms_envelope';
lfp_activity_window_policy = 'eigenvalue_adaptive_rms_envelope';
envelope_enable = true;
envelope_policy = 'eigenvalue_adaptive_rms';
envelope_alpha = 0.35;
envelope_min_window_sec = 0.03;
envelope_max_window_sec = 1.0;
envelope_fallback_window_sec = 0.10;

make_state_space_plot = false;
make_consensus_state_space_plot = false;
make_spectrum_diagnostics = false;
make_top30_window_plots = false;
dry_run = false;
skip_existing_density = true;

run(fullfile(repo_root, 'scripts', ...
    'script_run_current_p5_adaptive_envelope_density_from_existing_reductions_unblocked.m'));
