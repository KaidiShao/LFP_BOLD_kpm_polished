% Batch wrapper for current P5 adaptive RMS-envelope peak statistics.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'e10gw1', 'f12m01'};
variant_filter = { ...
    'abs_projected_vlambda', ...
    'complex_split_projected_vlambda'};
target_variant_suffix = 'rmsenv_adaptive';
exclude_dataset_variant_pairs = {'e10fV1|complex_split_projected_vlambda'};

peak_mode = 'max_abs';
peak_activity_suffix = 'rmsenv_adaptive';
density_value_transform = 'abs';
lfp_activity_transform = 'rms_envelope';
lfp_activity_window_policy = 'eigenvalue_adaptive_rms_envelope';
envelope_enable = true;
envelope_policy = 'eigenvalue_adaptive_rms';
envelope_alpha = 0.35;
envelope_min_window_sec = 0.03;
envelope_max_window_sec = 1.0;
envelope_fallback_window_sec = 0.10;

force_recompute = false;
save_figures = false;
continue_on_error = true;

run(fullfile(repo_root, 'scripts', ...
    'script_run_current_pipeline5_peak_state_all_datasets.m'));
