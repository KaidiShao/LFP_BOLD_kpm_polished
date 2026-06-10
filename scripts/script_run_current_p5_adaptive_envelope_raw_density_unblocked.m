% Batch wrapper for P5 adaptive RMS-envelope raw density only.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

dataset_filter = {};
condition_filter = {};
method_filter = {'svd', 'nmf', 'mds', 'umap'};
component_count_sweep = 3:8;
activity_suffix = 'rmsenv_adaptive';
run_raw_density = true;
run_dimred_density = false;
skip_existing_density = true;
dry_run = false;

run(fullfile(repo_root, 'scripts', ...
    'script_run_current_p5_adaptive_envelope_density_from_existing_reductions_unblocked.m'));
