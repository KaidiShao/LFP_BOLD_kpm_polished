% Run the currently missing completed-model postprocessing conditions.
%
% This intentionally skips the old event-diversity branch. It runs the
% top-state-diversity-window postprocessing stages:
%   main postprocess, timescale diagnostics, deconv residuals,
%   window-normalized deconv residuals, and full-time spike correlation.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

params = struct();
params.dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01'};
params.exclude_dataset_stems = {};
params.autodl_root = 'E:\autodl_results';
params.processed_root = io_project.get_project_processed_root();
params.output_folder_name = 'mlp_top_state_diversity_postprocessing';

% Remaining conditions after the 2026-04-23 E10/F12 postprocessing runs.
% E10gH1 abs-vlambda was exported later, so it still needs this stage.
% Keep abs-kv in the list: its top-window plots already exist, but the
% spike-correlation stage still needs to be generated. Existing window PNGs
% are skipped by run_mlp_top_state_diversity_postprocessing_pipeline.
% Format: dataset|observable_mode|residual_form
params.condition_key_filter = { ...
    'e10gh1|abs|projected_vlambda', ...
    'f12m01|abs|projected_kv', ...
    'f12m01|abs|projected_vlambda', ...
    'f12m01|complex_split|projected_kv', ...
    'f12m01|complex_split|projected_vlambda'};

params.n_top_windows = 30;
params.window_length_samples = 6000;
params.window_mode = 'global';

params.max_basis = 30;
params.timescale_max_modes_all = 30;
params.timescale_max_modes_sel = 20;
params.timescale_max_lag = 200;

params.do_main_plot = true;
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

fprintf('Running missing completed-model postprocessing conditions.\n');
fprintf('Condition filter:\n');
for i = 1:numel(params.condition_key_filter)
    fprintf('  - %s\n', params.condition_key_filter{i});
end

manifest = run_mlp_top_state_diversity_postprocessing_pipeline(params);

fprintf('\nFinished missing completed-model postprocessing.\n');
fprintf('Window rows: %d\n', height(manifest.table));
if isfield(manifest, 'spike_table')
    fprintf('Spike rows: %d\n', height(manifest.spike_table));
end
fprintf('CSV manifest:\n  %s\n', manifest.csv_file);
fprintf('MAT manifest:\n  %s\n', manifest.mat_file);
if isfield(manifest, 'spike_csv_file') && ~isempty(manifest.spike_csv_file)
    fprintf('Spike manifest:\n  %s\n', manifest.spike_csv_file);
end
