% Run pre-reduction MLP/EDMD postprocessing on completed E10 conditions.
%
% Windows are always taken from each dataset's saved consensus-state-diversity
% top-30 global windows. Outputs are written under each dataset's
% postprocessing folder. F12m01 is intentionally excluded.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

params = struct();
params.dataset_stems = {'e10fV1', 'e10gb1', 'e10gh1'};
params.exclude_dataset_stems = {'f12m01'};
params.autodl_root = 'E:\autodl_results';
params.processed_root = get_project_processed_root();
params.output_folder_name = 'mlp_top_state_diversity_postprocessing';

params.n_top_windows = 30;
params.window_length_samples = 6000;
params.window_mode = 'global';

params.max_basis = 30;
params.timescale_max_modes_all = 30;
params.timescale_max_modes_sel = 20;
params.timescale_max_lag = 200;

params.do_timescale = true;
params.do_deconv = true;
params.deconv_method = 'koopman_residual';
params.deconv_lambda_source = 'edmd';
params.deconv_plot_normalize_scope = 'global';
params.do_deconv_window_norm = true;
params.deconv_window_plot_normalize_scope = 'window';
params.deconv_max_modes_all = 30;
params.deconv_max_modes_sel = 20;

params.do_spike_correlation = true;
params.spike_corr_channels = 'all';
params.spike_corr_max_modes = 20;
params.spike_corr_top_n_rows = 100;
params.spike_corr_progress_every = 50;
params.spike_corr_skip_existing = true;
params.spike_corr_save_overview = true;
params.spike_corr_overview_features = {'abs_mean', 'abs_rms'};
params.spike_corr_save_under_run = false;
params.spike_corr_save_dir_name = 'mlp_fulltime_spike_residual_correlation';

params.save_png = true;
params.save_fig = false;
params.skip_existing = true;
params.close_figures = true;
params.headless = true;
params.resolution = 180;

manifest = run_mlp_top_state_diversity_postprocessing_pipeline(params);

fprintf('\nFinished completed MLP top-window postprocessing.\n');
fprintf('Rows: %d\n', height(manifest.table));
fprintf('CSV manifest:\n  %s\n', manifest.csv_file);
fprintf('MAT manifest:\n  %s\n', manifest.mat_file);
