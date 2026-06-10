function params = build_blp_eigenfunction_postprocessing_params()
%BUILD_BLP_EIGENFUNCTION_POSTPROCESSING_PARAMS Default params for pipeline 6.

params = struct();
params.dataset_stems = {'e10fV1', 'e10gb1', 'e10gh1'};
params.exclude_dataset_stems = {'f12m01'};
params.autodl_root = 'E:\autodl_results';
params.processed_root = io_project.get_project_processed_root();
params.output_root = [];
params.output_folder_name = io_project.get_pipeline_stage_name(6, 'top_state_diversity_postprocessing');
params.filename_pattern = '*_outputs_*.mat';
params.variable_name = 'EDMD_outputs';
params.run_name_filter = {};
params.condition_key_filter = {};
params.invalid_run_registry_file = fullfile( ...
    fileparts(fileparts(fileparts(mfilename('fullpath')))), ...
    'configs', 'current_invalid_blp_mlp_runs.csv');
params.include_invalid_runs = false;

params.n_top_windows = 30;
params.window_length_samples = 6000;
params.window_mode = 'global';

params.abs_thresh = 0.01;
params.sort_by = 'modulus';
params.sort_dir = 'descend';
params.max_basis = 30;
params.timescale_max_modes_all = Inf;
params.timescale_max_modes_sel = params.max_basis;
params.timescale_max_lag = [];
params.timescale_xlim_time = 20;
params.timescale_match_empirical_scale = true;
params.timescale_similarity_max_samples = 100000;

params.do_main_plot = true;
params.do_timescale = true;
params.do_deconv = true;
params.deconv_method = 'koopman_residual';
params.deconv_lambda_source = 'edmd';
params.deconv_plot_normalize_scope = 'global';
params.deconv_normalize_exclude_idx = 1;
params.do_deconv_window_norm = true;
params.deconv_window_plot_normalize_scope = 'window';
params.deconv_max_modes_all = Inf;
params.deconv_max_modes_sel = params.max_basis;

params.do_spkt_cross_correlation = true;
params.spkt_cross_channels = 'all';
params.spkt_cross_max_modes = 20;
params.spkt_cross_top_n_rows = 100;
params.spkt_cross_progress_every = 50;
params.spkt_cross_skip_existing = true;
params.spkt_cross_save_overview = true;
params.spkt_cross_overview_features = {'abs_mean', 'abs_rms'};
params.spkt_cross_overview_resolution = 220;
params.spkt_cross_overview_top_n = 30;
params.spkt_cross_verbose = true;
params.spkt_cross_save_under_run = false;
params.spkt_cross_save_dir_name = io_project.get_pipeline_stage_name(6, 'spkt_residual_cross_correlation');

params.do_mua_cross_correlation = true;
params.mua_cross_channels = 'selected';
params.mua_cross_band = 'last';
params.mua_cross_pairings = {'abs_abs', 'raw_real', 'raw_imag'};
params.mua_cross_max_modes = 20;
params.mua_cross_top_n_rows = 100;
params.mua_cross_progress_every = 50;
params.mua_cross_skip_existing = true;
params.mua_cross_save_overview = true;
params.mua_cross_overview_resolution = 220;
params.mua_cross_verbose = true;
params.mua_cross_save_under_run = false;
params.mua_cross_save_dir_name = io_project.get_pipeline_stage_name(6, 'mua_residual_cross_correlation');

params.save_png = true;
params.save_fig = false;
params.skip_existing = true;
params.close_figures = true;
params.headless = true;
params.resolution = 180;
end
