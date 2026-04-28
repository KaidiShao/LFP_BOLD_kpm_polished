% Run raw-eigenfunction thresholded density for all complete MLP ResKoopNet runs.
%
% This script scans thresholds on the original eigenfunctions only. It does
% not run eigenfunction dimension reduction and does not compute density on
% dimension-reduced components.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
end
set(groot, 'defaultFigureVisible', 'off');

params = struct();
params.autodl_root = 'E:\autodl_results';
params.processed_root = io_project.get_project_processed_root();

% Empty dataset_stems means: discover every dataset folder under autodl_root.
% Set this to e.g. {'e10fV1','e10gb1','e10gh1'} to restrict the run.
params.dataset_stems = {};
params.exclude_dataset_stems = {};

params.output_folder_name = 'efun_den';
params.legacy_output_folder_names = {'eigenfunction_density_scan'};
params.threshold_ratios = [0.1:0.1:0.9, 0.99];
params.threshold_mode = 'quantile';

params.window_sec = 2;
params.step_sec = 2;
params.feature_variant = 'abs';
params.value_transform = 'none';
params.abs_thresh = 0.01;
params.sort_by = 'modulus';
params.sort_dir = 'descend';
params.max_modes = Inf;

params.require_session_metadata = true;
params.density_denominator = 'window_samples';
params.output_class = 'single';
params.smooth_density = false;
params.smooth_window_sec = 0;

params.save_results = true;
params.make_figure = true;
params.save_figure = true;
params.max_plot_modes = 100;
params.figure_resolution = 200;

params.skip_existing = true;
params.continue_on_error = true;
params.headless = true;
params.close_figures = true;
params.verbose = true;
params.progress_every = 10;
params.progress_timing = true;

manifest = run_reskoopnet_eigenfunction_density_batch(params);

fprintf('\nFinished raw-eigenfunction density threshold scan.\n');
fprintf('Rows: %d\n', height(manifest.table));
fprintf('CSV manifest:\n  %s\n', manifest.csv_file);
fprintf('MAT manifest:\n  %s\n', manifest.mat_file);
