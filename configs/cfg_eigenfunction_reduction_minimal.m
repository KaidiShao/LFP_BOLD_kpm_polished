function cfg = cfg_eigenfunction_reduction_minimal()
%CFG_EIGENFUNCTION_REDUCTION_MINIMAL Minimal config for eigenfunction reduction.

this_cfg_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_cfg_dir);

dataset_name = 'e10gb1';
observable_tag = 'abs';
source_run_name = 'mlp_obs_e10gb1_260415_shuffle_plateau_projected_kv_abs';
output_run_name = 'mlp260415_kvabs';
processed_root = get_project_processed_root();
output_root = fullfile(processed_root, dataset_name, ...
    'efun', output_run_name);
figure_output_dir = fullfile(output_root, 'fig');
result_output_dir = fullfile(output_root, 'mat');
observable_file = fullfile(processed_root, dataset_name, ...
    'reskoopnet_dictionary', ...
    sprintf('%s_low50_high250_g2_%s_single.mat', dataset_name, observable_tag));

cfg = struct();
cfg.repo_root = repo_root;
cfg.dataset = struct();
cfg.dataset.name = dataset_name;
cfg.dataset.processed_root = processed_root;

cfg.output = struct();
cfg.output.root = output_root;
cfg.output.figure_dir = figure_output_dir;
cfg.output.result_dir = result_output_dir;
cfg.output.source_run_name = source_run_name;
cfg.output.output_run_name = output_run_name;

cfg.source = struct();
cfg.source.mode = 'chunk_dir';  % 'chunk_dir' | 'mat_file'
cfg.source.data_dir = fullfile('E:\autodl_results', dataset_name, ...
    'mlp', 'outputs', source_run_name);
cfg.source.edmd_file = '';
cfg.source.concat = struct();
cfg.source.concat.filename_pattern = '*_outputs_*.mat';
cfg.source.concat.variable_name = 'EDMD_outputs';
cfg.source.concat.concat_fields = {'efuns'};
cfg.source.concat.concat_dim = 1;
cfg.source.concat.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form', ...
    'dt', 'dx', 'sampling_period', 'sample_period', 'fs', 'sampling_frequency', ...
    'session_dx', 'session_fs', 'session_ids', 'session_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
cfg.source.concat.allow_missing_chunks = false;
cfg.source.concat.verbose = true;
cfg.source.concat.progress_every = 10;
cfg.source.concat.progress_timing = true;

cfg.input = struct();
cfg.input.dt = [];
cfg.input.observable_file = observable_file;

cfg.selection = struct();
cfg.selection.abs_thresh = 0.01;
cfg.selection.sort_by = 'modulus';   % 'modulus' | 'real'
cfg.selection.sort_dir = 'descend';  % 'ascend' | 'descend'
cfg.selection.max_modes = Inf;

cfg.feature = struct();
cfg.feature.family = 'eigenfunction';
cfg.feature.variant = 'abs';         % 'abs' | 'real'
cfg.feature.normalization = 'maxabs_per_mode';

cfg.path = struct();
cfg.path.kind = 'time';              % 'time' | 'spectrum'

cfg.path.time = struct();
cfg.path.time.method = 'SVD';        % 'SVD' | 'logSVD' | 'NMF'
cfg.path.time.n_components = 4;
cfg.path.time.center_modes = true;
cfg.path.time.log_epsilon = 1e-10;
cfg.path.time.nmf_max_iter = 1000;
cfg.path.time.nmf_nonnegative_strategy = 'shift_global';  % 'shift_global' | 'clip_zero'

cfg.path.spectrum = struct();
cfg.path.spectrum.method = 'MDS';    % 'MDS' | 'UMAP' | 'tSNE' | 'diffusion_map'
cfg.path.spectrum.output_dim = 3;
cfg.path.spectrum.n_components = 4;
cfg.path.spectrum.distance = 'corr_abs';
cfg.path.spectrum.downsample_step = 10;
cfg.path.spectrum.cluster_method = 'kmeans';   % 'kmeans' | 'GMM'
cfg.path.spectrum.cluster_replicates = 10;
cfg.path.spectrum.gmm_replicates = 3;
cfg.path.spectrum.gmm_regularization = 1e-5;
cfg.path.spectrum.mds_max_iter = 1000;
cfg.path.spectrum.tsne_perplexity = 30;
cfg.path.spectrum.umap_dir = '';
cfg.path.spectrum.sign_alignment = 'dominant_cluster';  % 'none' | 'dominant_cluster'
cfg.path.spectrum.reconstruction_weight_method = 'least_squares';  % 'least_squares' | 'nonnegative'

cfg.summary = struct();
cfg.summary.smooth = struct();
cfg.summary.smooth.enable = true;
cfg.summary.smooth.method = 'movmean';
cfg.summary.smooth.window = 10;

cfg.viz = struct();
cfg.viz.overview = struct();
cfg.viz.overview.enable = true;
cfg.viz.overview.max_modes = 50;
cfg.viz.overview.window_idx = [];
cfg.viz.overview.component_source = 'smooth_if_available';  % 'raw' | 'smooth' | 'smooth_if_available'
cfg.viz.overview.plot_raw_under_smooth = true;
cfg.viz.overview.figure_position = [120, 80, 1180, 820];
cfg.viz.overview.background_color = [0, 0, 0];
cfg.viz.overview.axes_color = [0, 0, 0];
cfg.viz.overview.grid_color = [0.75, 0.75, 0.75];
cfg.viz.overview.text_color = [1, 1, 1];
cfg.viz.overview.raw_line_color = [0.55, 0.55, 0.55];
cfg.viz.overview.component_line_color = [0.96, 0.96, 0.96];
cfg.viz.overview.event_windows = [];
cfg.viz.overview.event_colors = [
    0.96, 0.84, 0.62;
    0.70, 0.80, 0.97;
    0.88, 0.72, 0.80;
    0.85, 0.72, 0.95];
cfg.viz.overview.event_alpha = 0.38;
cfg.viz.overview.title = 'Eigenfunctions And Temporal Components';
cfg.viz.overview.save_figure = true;
cfg.viz.overview.save_dir = figure_output_dir;
cfg.viz.overview.save_tag = 'ov';

cfg.viz.spectrum = struct();
cfg.viz.spectrum.enable = true;
cfg.viz.spectrum.figure_position = [80, 120, 1260, 360];
cfg.viz.spectrum.background_color = [0, 0, 0];
cfg.viz.spectrum.axes_color = [0, 0, 0];
cfg.viz.spectrum.grid_color = [0.75, 0.75, 0.75];
cfg.viz.spectrum.text_color = [1, 1, 1];
cfg.viz.spectrum.distance_colormap = turbo(256);
cfg.viz.spectrum.embedding_index_colormap = parula(256);
cfg.viz.spectrum.cluster_colormap = [];
cfg.viz.spectrum.marker_size = 28;
cfg.viz.spectrum.draw_unit_circle = true;
cfg.viz.spectrum.use_bilinear_evalues = false;
cfg.viz.spectrum.title = 'Spectrum Path Diagnostics';
cfg.viz.spectrum.save_figure = true;
cfg.viz.spectrum.save_dir = figure_output_dir;
cfg.viz.spectrum.save_tag = 'spec';

cfg.viz.state_space = struct();
cfg.viz.state_space.enable = true;
cfg.viz.state_space.window_idx = [];
cfg.viz.state_space.component_source = 'smooth_if_available';  % 'raw' | 'smooth' | 'smooth_if_available'
cfg.viz.state_space.figure_position = [120, 100, 1040, 860];
cfg.viz.state_space.background_color = [0, 0, 0];
cfg.viz.state_space.axes_color = [0, 0, 0];
cfg.viz.state_space.grid_color = [0.82, 0.82, 0.82];
cfg.viz.state_space.text_color = [1, 1, 1];
cfg.viz.state_space.line_width = 1.15;
cfg.viz.state_space.downsample_step = 10;
cfg.viz.state_space.max_plot_points = Inf;
cfg.viz.state_space.view_angle = [37, 24];
cfg.viz.state_space.time_colormap = parula(256);
cfg.viz.state_space.value_colormap = turbo(256);
cfg.viz.state_space.title = 'State-Space Trajectory';
cfg.viz.state_space.save_figure = true;
cfg.viz.state_space.figure_resolution = 180;
cfg.viz.state_space.verbose = true;
cfg.viz.state_space.save_dir = figure_output_dir;
cfg.viz.state_space.save_tag = 'ss';

cfg.viz.state_space_consensus = struct();
cfg.viz.state_space_consensus.enable = true;
cfg.viz.state_space_consensus.consensus_input = [];
cfg.viz.state_space_consensus.window_idx = [];
cfg.viz.state_space_consensus.component_source = 'smooth_if_available';  % 'raw' | 'smooth' | 'smooth_if_available'
cfg.viz.state_space_consensus.state_start_idx = 1;
cfg.viz.state_space_consensus.baseline_code = 0;
cfg.viz.state_space_consensus.baseline_label = 'baseline';
cfg.viz.state_space_consensus.figure_position = [140, 90, 1120, 860];
cfg.viz.state_space_consensus.background_color = [0, 0, 0];
cfg.viz.state_space_consensus.axes_color = [0, 0, 0];
cfg.viz.state_space_consensus.grid_color = [0.78, 0.78, 0.78];
cfg.viz.state_space_consensus.text_color = [1, 1, 1];
cfg.viz.state_space_consensus.baseline_line_width = 0.75;
cfg.viz.state_space_consensus.event_line_width = 2.0;
cfg.viz.state_space_consensus.patch_line_width = 1.15;
cfg.viz.state_space_consensus.downsample_step = 10;
cfg.viz.state_space_consensus.max_plot_points = Inf;
cfg.viz.state_space_consensus.single_sample_marker_size = 18;
cfg.viz.state_space_consensus.view_angle = [37, 24];
cfg.viz.state_space_consensus.state_colormap = [];
cfg.viz.state_space_consensus.show_state_number_labels = true;
cfg.viz.state_space_consensus.label_baseline = false;
cfg.viz.state_space_consensus.state_label_marker_size = 13;
cfg.viz.state_space_consensus.title = 'State-Space Trajectory By Consensus State';
cfg.viz.state_space_consensus.save_figure = true;
cfg.viz.state_space_consensus.figure_resolution = 200;
cfg.viz.state_space_consensus.verbose = true;
cfg.viz.state_space_consensus.save_dir = figure_output_dir;
cfg.viz.state_space_consensus.save_tag = 'ssc';

cfg.thresholded_density = struct();
cfg.thresholded_density.enable = true;
cfg.thresholded_density.window_sec = 2;
cfg.thresholded_density.step_sec = 2;
cfg.thresholded_density.threshold_ratio = 0.7;
cfg.thresholded_density.threshold_mode = 'quantile';  % 'quantile' | 'maxfrac' | 'meanplusstd'
cfg.thresholded_density.value_transform = 'none';        % feature.variant already controls abs/real
cfg.thresholded_density.density_denominator = 'window_samples';
cfg.thresholded_density.smooth_density = false;
cfg.thresholded_density.smooth_window_sec = 0;
cfg.thresholded_density.output_class = 'single';
cfg.thresholded_density.require_session_metadata = true;
cfg.thresholded_density.observable_file = observable_file;
cfg.thresholded_density.save_results = true;
cfg.thresholded_density.make_figure = true;
cfg.thresholded_density.save_figure = true;
cfg.thresholded_density.close_figure = true;
cfg.thresholded_density.save_dir = '';
cfg.thresholded_density.figure_dir = '';
cfg.thresholded_density.save_stem = sprintf('%s_thresholded_density', dataset_name);
cfg.thresholded_density.save_tag = output_run_name;
cfg.thresholded_density.max_plot_modes = 100;
cfg.thresholded_density.title = 'Thresholded Eigenfunction Density';

cfg.thresholded_events = struct();
cfg.thresholded_events.enable = false;
cfg.thresholded_events.threshold_ratio = cfg.thresholded_density.threshold_ratio;
cfg.thresholded_events.threshold_mode = cfg.thresholded_density.threshold_mode;
cfg.thresholded_events.value_transform = cfg.thresholded_density.value_transform;
cfg.thresholded_events.event_detector = 'find_peak_loc';  % 'find_peak_loc' | 'runs'
cfg.thresholded_events.event_score = 'area_above_threshold';
cfg.thresholded_events.min_duration_sec = 0.03;
cfg.thresholded_events.merge_gap_sec = 0.02;
cfg.thresholded_events.find_peak_loc_window_sec = [];
cfg.thresholded_events.find_peak_loc_window_samples = [];
cfg.thresholded_events.find_peak_loc_drop_first = false;
cfg.thresholded_events.find_peak_loc_post_nms = false;
cfg.thresholded_events.bin_sec = 2;
cfg.thresholded_events.step_sec = 2;
cfg.thresholded_events.smooth_rate = true;
cfg.thresholded_events.smooth_sigma_sec = 2;
cfg.thresholded_events.output_class = 'single';
cfg.thresholded_events.require_session_metadata = true;
cfg.thresholded_events.observable_file = observable_file;
cfg.thresholded_events.nms = struct();
cfg.thresholded_events.nms.enable = true;
cfg.thresholded_events.nms.mode = 'eigen_period';  % 'eigen_period' | 'fixed'
cfg.thresholded_events.nms.period_fraction = 1.0;
cfg.thresholded_events.nms.fixed_sec = 0.25;
cfg.thresholded_events.nms.default_sec = 0.25;
cfg.thresholded_events.nms.min_sec = 0.05;
cfg.thresholded_events.nms.max_sec = 2.0;
cfg.thresholded_events.nms.angle_eps = 1e-4;
cfg.thresholded_events.save_results = false;
cfg.thresholded_events.make_figure = false;
cfg.thresholded_events.save_figure = false;
cfg.thresholded_events.close_figure = true;
cfg.thresholded_events.save_dir = '';
cfg.thresholded_events.figure_dir = '';
cfg.thresholded_events.save_stem = sprintf('%s_thresholded_events', dataset_name);
cfg.thresholded_events.save_tag = output_run_name;
cfg.thresholded_events.max_plot_modes = 80;
cfg.thresholded_events.max_plot_events = 5000;
cfg.thresholded_events.title = 'Thresholded Eigenfunction Events';

cfg.dimred_thresholded_density = struct();
cfg.dimred_thresholded_density.enable = false;
cfg.dimred_thresholded_density.window_sec = cfg.thresholded_density.window_sec;
cfg.dimred_thresholded_density.step_sec = cfg.thresholded_density.step_sec;
cfg.dimred_thresholded_density.threshold_ratio = cfg.thresholded_density.threshold_ratio;
cfg.dimred_thresholded_density.threshold_mode = cfg.thresholded_density.threshold_mode;
cfg.dimred_thresholded_density.value_transform = 'none';
cfg.dimred_thresholded_density.density_denominator = cfg.thresholded_density.density_denominator;
cfg.dimred_thresholded_density.smooth_density = cfg.thresholded_density.smooth_density;
cfg.dimred_thresholded_density.smooth_window_sec = cfg.thresholded_density.smooth_window_sec;
cfg.dimred_thresholded_density.output_class = cfg.thresholded_density.output_class;
cfg.dimred_thresholded_density.require_session_metadata = cfg.thresholded_density.require_session_metadata;
cfg.dimred_thresholded_density.observable_file = observable_file;
cfg.dimred_thresholded_density.save_results = false;
cfg.dimred_thresholded_density.make_figure = false;
cfg.dimred_thresholded_density.save_figure = false;
cfg.dimred_thresholded_density.close_figure = true;
cfg.dimred_thresholded_density.save_dir = '';
cfg.dimred_thresholded_density.figure_dir = '';
cfg.dimred_thresholded_density.save_stem = sprintf('%s_dimred_thresholded_density', dataset_name);
cfg.dimred_thresholded_density.save_tag = output_run_name;
cfg.dimred_thresholded_density.max_plot_modes = 40;
cfg.dimred_thresholded_density.title = 'Thresholded Eigenfunction Component Density';

cfg.dimred_thresholded_events = struct();
cfg.dimred_thresholded_events.enable = false;
cfg.dimred_thresholded_events.threshold_ratio = cfg.thresholded_events.threshold_ratio;
cfg.dimred_thresholded_events.threshold_mode = cfg.thresholded_events.threshold_mode;
cfg.dimred_thresholded_events.value_transform = cfg.thresholded_events.value_transform;
cfg.dimred_thresholded_events.event_detector = cfg.thresholded_events.event_detector;
cfg.dimred_thresholded_events.event_score = cfg.thresholded_events.event_score;
cfg.dimred_thresholded_events.min_duration_sec = cfg.thresholded_events.min_duration_sec;
cfg.dimred_thresholded_events.merge_gap_sec = cfg.thresholded_events.merge_gap_sec;
cfg.dimred_thresholded_events.find_peak_loc_window_sec = cfg.thresholded_events.find_peak_loc_window_sec;
cfg.dimred_thresholded_events.find_peak_loc_window_samples = cfg.thresholded_events.find_peak_loc_window_samples;
cfg.dimred_thresholded_events.find_peak_loc_drop_first = cfg.thresholded_events.find_peak_loc_drop_first;
cfg.dimred_thresholded_events.find_peak_loc_post_nms = cfg.thresholded_events.find_peak_loc_post_nms;
cfg.dimred_thresholded_events.bin_sec = cfg.thresholded_events.bin_sec;
cfg.dimred_thresholded_events.step_sec = cfg.thresholded_events.step_sec;
cfg.dimred_thresholded_events.smooth_rate = cfg.thresholded_events.smooth_rate;
cfg.dimred_thresholded_events.smooth_sigma_sec = cfg.thresholded_events.smooth_sigma_sec;
cfg.dimred_thresholded_events.output_class = cfg.thresholded_events.output_class;
cfg.dimred_thresholded_events.require_session_metadata = cfg.thresholded_events.require_session_metadata;
cfg.dimred_thresholded_events.observable_file = observable_file;
cfg.dimred_thresholded_events.component_evalue_weight_transform = 'abs';
cfg.dimred_thresholded_events.component_evalue_angle_statistic = 'weighted_median_abs';
cfg.dimred_thresholded_events.component_evalue_top_modes = 5;
cfg.dimred_thresholded_events.nms = cfg.thresholded_events.nms;
cfg.dimred_thresholded_events.save_results = false;
cfg.dimred_thresholded_events.make_figure = false;
cfg.dimred_thresholded_events.save_figure = false;
cfg.dimred_thresholded_events.close_figure = true;
cfg.dimred_thresholded_events.save_dir = '';
cfg.dimred_thresholded_events.figure_dir = '';
cfg.dimred_thresholded_events.save_stem = sprintf('%s_dimred_thresholded_events', dataset_name);
cfg.dimred_thresholded_events.save_tag = output_run_name;
cfg.dimred_thresholded_events.max_plot_modes = 40;
cfg.dimred_thresholded_events.max_plot_events = cfg.thresholded_events.max_plot_events;
cfg.dimred_thresholded_events.title = 'Thresholded Eigenfunction Component Events';

cfg.save = struct();
cfg.save.enable = true;
cfg.save.dir = result_output_dir;
cfg.save.file_stem = sprintf('%s_efun', dataset_name);
cfg.save.tag = 'min';
cfg.save.payload = 'compact';  % 'compact' saves derived outputs only; use 'full' for debugging.
cfg.save.v7_3 = true;
end
