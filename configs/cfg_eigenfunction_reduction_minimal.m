function cfg = cfg_eigenfunction_reduction_minimal()
%CFG_EIGENFUNCTION_REDUCTION_MINIMAL Minimal config for eigenfunction reduction.

this_cfg_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_cfg_dir);
results_root = get_project_results_root(repo_root);

cfg = struct();
cfg.repo_root = repo_root;

cfg.source = struct();
cfg.source.mode = 'chunk_dir';  % 'chunk_dir' | 'mat_file'
cfg.source.data_dir = 'E:\autodl_results\initial_point_test1_KV';
cfg.source.edmd_file = '';
cfg.source.concat = struct();
cfg.source.concat.filename_pattern = '*_outputs_*.mat';
cfg.source.concat.variable_name = 'EDMD_outputs';
cfg.source.concat.concat_fields = {'efuns'};
cfg.source.concat.concat_dim = 1;
cfg.source.concat.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form'};
cfg.source.concat.allow_missing_chunks = false;
cfg.source.concat.verbose = true;
cfg.source.concat.progress_every = 50;

cfg.input = struct();
cfg.input.dt = [];

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
cfg.path.kind = 'spectrum';          % 'time' | 'spectrum'

cfg.path.time = struct();
cfg.path.time.method = 'SVD';        % 'SVD' | 'logSVD' | 'NMF' | 'NMF2'
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
cfg.viz.overview.save_figure = false;
cfg.viz.overview.save_dir = fullfile(results_root, 'eigenfunction_reduction');
cfg.viz.overview.save_tag = 'overview';

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
cfg.viz.spectrum.save_figure = false;
cfg.viz.spectrum.save_dir = fullfile(results_root, 'eigenfunction_reduction');
cfg.viz.spectrum.save_tag = 'spectrum_diag';

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
cfg.viz.state_space.view_angle = [37, 24];
cfg.viz.state_space.time_colormap = parula(256);
cfg.viz.state_space.value_colormap = turbo(256);
cfg.viz.state_space.title = 'State-Space Trajectory';
cfg.viz.state_space.save_figure = false;
cfg.viz.state_space.save_dir = fullfile(results_root, 'eigenfunction_reduction');
cfg.viz.state_space.save_tag = 'state_space';
cfg.save = struct();
cfg.save.enable = true;
cfg.save.dir = fullfile(results_root, 'eigenfunction_reduction');
cfg.save.file_stem = 'eigenfunction_reduction_result';
cfg.save.tag = 'minimal';
cfg.save.v7_3 = true;
end
