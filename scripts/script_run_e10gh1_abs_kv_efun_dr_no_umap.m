% Run E10.gH1 abs/projected_kv eigenfunction dimension reduction.
%
% This uses the completed fixed E10.gH1 abs-kv MLP run and intentionally
% excludes UMAP for now. It reuses the shared batch runner so the outputs
% land in the same efun folder layout as the other datasets.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
end
set(groot, 'defaultFigureVisible', 'off');

efun_batch_override = struct();
efun_batch_override.dataset_filter = {'e10gh1'};
efun_batch_override.run_key_filter = {'e10gh1/abs_kv'};
efun_batch_override.method_filter = {'svd', 'logsvd', 'nmf', 'mds'};
efun_batch_override.sweep_name = 'dr_no_umap';
efun_batch_override.master_manifest_file = fullfile( ...
    get_project_processed_root(), 'e10gh1_abs_kv_efun_mlp_dr_no_umap.csv');

efun_batch_override.continue_on_error = true;
efun_batch_override.skip_unavailable_methods = true;
efun_batch_override.close_figures = true;

efun_batch_override.make_overview = true;
efun_batch_override.make_state_space = true;
efun_batch_override.make_consensus_state_space = true;
efun_batch_override.make_spectrum_diagnostics = true;
efun_batch_override.make_top30 = true;
efun_batch_override.make_peak_analysis = true;

efun_batch_override.top30_n_windows = 30;
efun_batch_override.top30_force_recompute_norm_cache = false;

fprintf('Running E10.gH1 abs-kv eigenfunction DR without UMAP.\n');
fprintf('Manifest:\n  %s\n', efun_batch_override.master_manifest_file);

run(fullfile(this_script_dir, 'script_run_all_complete_mlp_efun_dr_pipeline.m'));

clear efun_batch_override
