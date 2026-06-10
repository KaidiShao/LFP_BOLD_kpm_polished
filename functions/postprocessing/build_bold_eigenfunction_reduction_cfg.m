function [cfg, method_tag, feature_tag] = build_bold_eigenfunction_reduction_cfg( ...
        candidate, feature_name, spec, params, repo_root)
%BUILD_BOLD_EIGENFUNCTION_REDUCTION_CFG Build one P9 method config.

if nargin < 5
    repo_root = '';
end

method_tag = build_blp_eigenfunction_reduction_method_tag(spec);
feature_tag = local_feature_tag(feature_name);

run_root = fullfile( ...
    io_project.get_pipeline_stage_dir(params.processed_root, ...
    candidate.dataset_stem, 9, 'bold_eigenfunction_reduction'), ...
    candidate.run_tag);
feature_root = fullfile(run_root, feature_tag);
method_root = fullfile(feature_root, method_tag);

cfg = struct();
cfg.repo_root = repo_root;

cfg.dataset = struct();
cfg.dataset.name = candidate.dataset_stem;
cfg.dataset.dataset_id = candidate.dataset_id;
cfg.dataset.processed_root = params.processed_root;

cfg.output = struct();
cfg.output.root = method_root;
cfg.output.run_root = run_root;
cfg.output.feature_root = feature_root;
cfg.output.figure_dir = fullfile(method_root, 'fig');
cfg.output.result_dir = fullfile(method_root, 'mat');
cfg.output.run_tag = candidate.run_tag;
cfg.output.run_name = candidate.run_name;
cfg.output.feature_tag = feature_tag;
cfg.output.method_tag = method_tag;

cfg.source = struct();
cfg.source.bold_post_file = candidate.bold_post_file;

cfg.feature = struct();
cfg.feature.name = char(string(feature_name));
cfg.feature.normalization = params.feature_normalization;

cfg.selection = params.selection;

cfg.path = struct();
cfg.path.kind = spec.kind;
cfg.path.time = struct();
cfg.path.time.method = 'SVD';
cfg.path.time.n_components = 4;
cfg.path.time.center_modes = true;
cfg.path.time.log_epsilon = 1e-10;
cfg.path.time.nmf_max_iter = 1000;
cfg.path.time.nmf_nonnegative_strategy = 'shift_global';
cfg.path.spectrum = struct();
cfg.path.spectrum.method = 'MDS';
cfg.path.spectrum.output_dim = 3;
cfg.path.spectrum.n_components = 4;
cfg.path.spectrum.distance = 'corr_abs';
cfg.path.spectrum.downsample_step = 10;
cfg.path.spectrum.cluster_method = 'kmeans';
cfg.path.spectrum.cluster_replicates = 10;
cfg.path.spectrum.gmm_replicates = 3;
cfg.path.spectrum.gmm_regularization = 1e-5;
cfg.path.spectrum.mds_max_iter = 1000;
cfg.path.spectrum.tsne_perplexity = 30;
cfg.path.spectrum.umap_dir = find_default_umap_dir(repo_root);
cfg.path.spectrum.sign_alignment = 'dominant_cluster';
cfg.path.spectrum.reconstruction_weight_method = 'least_squares';

switch lower(spec.kind)
    case 'time'
        cfg.path.time.method = spec.method;
        cfg.path.time = local_apply_options(cfg.path.time, spec.options);
    case 'spectrum'
        cfg.path.spectrum.method = spec.method;
        cfg.path.spectrum = local_apply_options(cfg.path.spectrum, spec.options);
    otherwise
        error('Unsupported reduction kind: %s.', spec.kind);
end

cfg.save = struct();
cfg.save.enable = true;
cfg.save.dir = cfg.output.result_dir;
cfg.save.file_stem = sprintf('%s_%s_%s', ...
    candidate.dataset_stem, candidate.run_tag, feature_tag);
cfg.save.tag = method_tag;
cfg.save.payload = params.save_payload;
cfg.save.v7_3 = params.save_v7_3;

cfg.plot = params.plot;
cfg.plot.enable = params.make_summary_plot;
cfg.plot.save_dir = cfg.output.figure_dir;
cfg.plot.save_tag = sprintf('%s_%s', feature_tag, method_tag);

cfg.progress = struct();
cfg.progress.verbose = true;
end


function target = local_apply_options(target, options)
if isempty(options) || isempty(fieldnames(options))
    return;
end
names = fieldnames(options);
for i = 1:numel(names)
    target.(names{i}) = options.(names{i});
end
end


function tag = local_feature_tag(feature_name)
tag = lower(char(string(feature_name)));
tag = regexprep(tag, '[^\w\-]+', '_');
end
