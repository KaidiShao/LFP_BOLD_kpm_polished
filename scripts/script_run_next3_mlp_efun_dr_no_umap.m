% Run the next three available complete MLP ResKoopNet efun DR pipelines.
%
% Conditions:
%   e10gb1 / vlambda
%   e10gh1 / cs01_kv
%   e10gh1 / cs01_vlambda
%
% UMAP is intentionally excluded for this batch. Outputs stay under:
%   E:\DataPons_processed\<dataset>\efun\<condition>\<method>\

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

processed_root = io_project.get_project_processed_root();

efun_batch_override = struct();
efun_batch_override.dataset_filter = {'e10gb1', 'e10gh1'};
efun_batch_override.output_run_filter = { ...
    'vlambda', ...
    'cs01_kv', ...
    'cs01_vlambda'};
efun_batch_override.run_key_filter = { ...
    'e10gb1/vlambda', ...
    'e10gh1/cs01_kv', ...
    'e10gh1/cs01_vlambda'};
efun_batch_override.method_filter = {'svd', 'logsvd', 'nmf', 'mds'};
efun_batch_override.master_manifest_file = fullfile(processed_root, ...
    'next3_efun_mlp_dr_no_umap.csv');
efun_batch_override.sweep_name = 'dr_no_umap';
efun_batch_override.continue_on_error = true;
efun_batch_override.skip_unavailable_methods = true;
efun_batch_override.close_figures = true;

fprintf('Running next-three MLP efun DR pipeline without UMAP.\n');
fprintf('Datasets: %s\n', strjoin(string(efun_batch_override.dataset_filter), ', '));
fprintf('Conditions: %s\n', strjoin(string(efun_batch_override.output_run_filter), ', '));
fprintf('Methods: %s\n', strjoin(string(efun_batch_override.method_filter), ', '));
fprintf('Manifest:\n  %s\n\n', efun_batch_override.master_manifest_file);

run(fullfile(this_script_dir, 'script_run_all_complete_mlp_efun_dr_pipeline.m'));

clear efun_batch_override
