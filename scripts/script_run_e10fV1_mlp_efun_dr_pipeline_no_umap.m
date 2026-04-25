% Run e10fV1 complete MLP ResKoopNet outputs through efun DR, excluding UMAP.
%
% This script reuses script_run_all_complete_mlp_efun_dr_pipeline.m with a
% narrow override so only e10fV1 conditions and stable methods are run.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

processed_root = get_project_processed_root();

efun_batch_override = struct();
efun_batch_override.dataset_filter = {'e10fV1'};
efun_batch_override.output_run_filter = { ...
    'realabs_kv', ...
    'realabs_vlambda', ...
    'cs02_kv', ...
    'cs02_vlambda'};
efun_batch_override.method_filter = {'svd', 'logsvd', 'nmf', 'mds'};
efun_batch_override.master_manifest_file = fullfile(processed_root, ...
    'e10fV1', 'efun', 'e10fV1_dr_no_umap.csv');
efun_batch_override.skip_unavailable_methods = true;
efun_batch_override.continue_on_error = true; %#ok<STRNU>

run(fullfile(repo_root, 'scripts', ...
    'script_run_all_complete_mlp_efun_dr_pipeline.m'));

clear efun_batch_override
