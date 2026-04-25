% Run e10fV1 complete MLP ResKoopNet efun DR pipeline without UMAP.
%
% This wrapper reuses script_run_all_complete_mlp_efun_dr_pipeline.m and only
% narrows the run/method filters. Outputs stay under:
%   E:\DataPons_processed\e10fV1\efun\<condition>\<method>\

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
    'e10fV1_efun_mlp_dr_no_umap.csv');
efun_batch_override.sweep_name = 'dr_no_umap';
efun_batch_override.continue_on_error = true;
efun_batch_override.skip_unavailable_methods = true;
efun_batch_override.close_figures = true;

fprintf('Running e10fV1 MLP efun DR pipeline without UMAP.\n');
fprintf('Conditions: %s\n', strjoin(string(efun_batch_override.output_run_filter), ', '));
fprintf('Methods: %s\n', strjoin(string(efun_batch_override.method_filter), ', '));
fprintf('Manifest:\n  %s\n\n', efun_batch_override.master_manifest_file);

run(fullfile(this_script_dir, 'script_run_all_complete_mlp_efun_dr_pipeline.m'));

clear efun_batch_override
