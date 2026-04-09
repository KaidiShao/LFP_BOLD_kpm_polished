function cfg = cfg_eigenfunction_reduction_minimal()
%CFG_EIGENFUNCTION_REDUCTION_MINIMAL Minimal config for eigenfunction reduction.

this_cfg_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_cfg_dir);

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

cfg.save = struct();
cfg.save.enable = true;
cfg.save.dir = fullfile(repo_root, 'results', 'eigenfunction_reduction');
cfg.save.file_stem = 'eigenfunction_reduction_result';
cfg.save.tag = 'minimal';
cfg.save.v7_3 = true;
end
