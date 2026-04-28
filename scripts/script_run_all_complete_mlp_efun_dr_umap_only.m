% Run UMAP-only eigenfunction dimension reduction for all complete MLP
% ResKoopNet outputs currently listed in the shared batch runner.
%
% This wrapper intentionally runs only the UMAP method so previously
% completed SVD/logSVD/NMF/MDS outputs do not get recomputed.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
    % Older MATLAB releases do not expose this graphics default.
end
set(groot, 'defaultFigureVisible', 'off');

efun_batch_override = struct();
efun_batch_override.method_filter = {'umap'};
efun_batch_override.sweep_name = 'dr_umap';
efun_batch_override.master_manifest_file = fullfile( ...
    io_project.get_project_processed_root(), 'efun_mlp_dr_all_umap.csv');

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

fprintf('Running ALL complete MLP efun DR pipelines with UMAP only.\n');
fprintf('Manifest:\n  %s\n\n', efun_batch_override.master_manifest_file);

run(fullfile(this_script_dir, 'script_run_all_complete_mlp_efun_dr_pipeline.m'));

clear efun_batch_override
