function [cfg, method_tag] = build_blp_eigenfunction_reduction_cfg( ...
        dataset_cfg, run_info, spec, params, repo_root, ...
        EDMD_outputs, concat_info, source_info)
%BUILD_BLP_EIGENFUNCTION_REDUCTION_CFG Build one explicit per-method config.

cfg = build_blp_eigenfunction_reduction_defaults();
method_tag = build_blp_eigenfunction_reduction_method_tag(spec);
condition_tag = resolve_blp_eigenfunction_condition_output_tag(run_info, params);

dataset_name = char(run_info.dataset);
stage_root = io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'eigenfunction_reduction');
condition_root = fullfile(stage_root, condition_tag);
method_root = fullfile(condition_root, method_tag);
figure_dir = fullfile(method_root, 'fig');
result_dir = fullfile(method_root, 'mat');
raw_td_root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'raw_thresholded_density'), ...
    condition_tag, method_tag);
raw_te_root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'raw_thresholded_events'), ...
    condition_tag, method_tag);
dimred_td_root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'dimred_thresholded_density'), ...
    condition_tag, method_tag);
dimred_te_root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'dimred_thresholded_events'), ...
    condition_tag, method_tag);

observable_file = local_observable_file(params.processed_root, dataset_name, run_info);

cfg.repo_root = repo_root;
cfg.dataset.name = dataset_name;
cfg.dataset.processed_root = params.processed_root;
cfg.output.root = method_root;
cfg.output.figure_dir = figure_dir;
cfg.output.result_dir = result_dir;
cfg.output.source_run_name = run_info.run_name;
cfg.output.output_run_name = condition_tag;
cfg.output.condition_tag = condition_tag;

cfg.source.mode = 'preloaded';
cfg.source.data_dir = run_info.output_dir;
cfg.source.preloaded_EDMD_outputs = EDMD_outputs;
cfg.source.preloaded_concat_info = concat_info;
cfg.source.preloaded_source_info = source_info;
cfg.source.concat.verbose = false;

cfg.input.observable_file = observable_file;
cfg.feature.variant = params.feature_variant;
cfg.feature.normalization = params.feature_normalization;

cfg.path.kind = spec.kind;
switch lower(spec.kind)
    case 'time'
        cfg.path.time.method = spec.method;
        cfg.path.time = local_apply_options(cfg.path.time, spec.options);
    case 'spectrum'
        cfg.path.spectrum.method = spec.method;
        cfg.path.spectrum = local_apply_options(cfg.path.spectrum, spec.options);
    otherwise
        error('Unsupported reduction kind: %s', spec.kind);
end

if isempty(cfg.path.spectrum.umap_dir)
    cfg.path.spectrum.umap_dir = find_default_umap_dir(repo_root);
end

cfg.save.enable = true;
cfg.save.dir = result_dir;
cfg.save.file_stem = sprintf('%s_%s_efun', dataset_name, condition_tag);
cfg.save.tag = method_tag;
cfg.save.payload = params.save_payload;
cfg.save.v7_3 = params.save_v7_3;

cfg.viz.overview.enable = false;
cfg.viz.overview.save_dir = figure_dir;
cfg.viz.overview.save_tag = [method_tag '_ov'];

cfg.viz.state_space.enable = params.make_state_space_plot;
cfg.viz.state_space.save_dir = figure_dir;
cfg.viz.state_space.save_tag = [method_tag '_ss'];

cfg.viz.state_space_consensus.enable = params.make_consensus_state_space_plot;
cfg.viz.state_space_consensus.save_dir = figure_dir;
cfg.viz.state_space_consensus.save_tag = [method_tag '_ssc'];

cfg.viz.spectrum.enable = params.make_spectrum_diagnostics;
cfg.viz.spectrum.save_dir = figure_dir;
cfg.viz.spectrum.save_tag = [method_tag '_spec'];

cfg.thresholded_density.enable = false;
cfg.thresholded_density.observable_file = observable_file;
cfg.thresholded_density.threshold_mode = params.threshold_mode;
cfg.thresholded_density.threshold_ratio = params.threshold_ratio;
cfg.thresholded_density.threshold_ratio_sweep = params.threshold_ratio_sweep;
cfg.thresholded_density.save_results = false;
cfg.thresholded_density.make_figure = false;
cfg.thresholded_density.save_figure = false;
cfg.thresholded_density.save_dir = fullfile(raw_td_root, 'mat');
cfg.thresholded_density.figure_dir = fullfile(raw_td_root, 'fig');
cfg.thresholded_density.save_stem = sprintf('%s_thresholded_density', dataset_name);
cfg.thresholded_density.save_tag = [condition_tag '_' method_tag];

cfg.thresholded_events.enable = false;
cfg.thresholded_events.observable_file = observable_file;
cfg.thresholded_events.threshold_mode = params.threshold_mode;
cfg.thresholded_events.threshold_ratio = params.threshold_ratio;
cfg.thresholded_events.threshold_ratio_sweep = [];
cfg.thresholded_events.save_results = false;
cfg.thresholded_events.make_figure = false;
cfg.thresholded_events.save_figure = false;
cfg.thresholded_events.save_dir = fullfile(raw_te_root, 'mat');
cfg.thresholded_events.figure_dir = fullfile(raw_te_root, 'fig');
cfg.thresholded_events.save_stem = sprintf('%s_thresholded_events', dataset_name);
cfg.thresholded_events.save_tag = [condition_tag '_' method_tag];

cfg.dimred_thresholded_density.enable = params.do_dimred_thresholded_density;
cfg.dimred_thresholded_density.observable_file = observable_file;
cfg.dimred_thresholded_density.threshold_mode = params.threshold_mode;
cfg.dimred_thresholded_density.threshold_ratio = params.threshold_ratio;
cfg.dimred_thresholded_density.threshold_ratio_sweep = params.threshold_ratio_sweep;
cfg.dimred_thresholded_density.save_results = local_stage_save_flag( ...
    params, 'save_dimred_thresholded_density_results', ...
    params.do_dimred_thresholded_density);
cfg.dimred_thresholded_density.make_figure = params.do_dimred_thresholded_density;
cfg.dimred_thresholded_density.save_figure = params.do_dimred_thresholded_density;
cfg.dimred_thresholded_density.save_dir = fullfile(dimred_td_root, 'mat');
cfg.dimred_thresholded_density.figure_dir = fullfile(dimred_td_root, 'fig');
cfg.dimred_thresholded_density.save_stem = sprintf('%s_dimred_thresholded_density', dataset_name);
cfg.dimred_thresholded_density.save_tag = [condition_tag '_' method_tag];

cfg.dimred_thresholded_events.enable = params.do_dimred_thresholded_events;
cfg.dimred_thresholded_events.observable_file = observable_file;
cfg.dimred_thresholded_events.threshold_mode = params.threshold_mode;
cfg.dimred_thresholded_events.threshold_ratio = params.threshold_ratio;
cfg.dimred_thresholded_events.threshold_ratio_sweep = [];
cfg.dimred_thresholded_events.save_results = false;
cfg.dimred_thresholded_events.make_figure = params.do_dimred_thresholded_events;
cfg.dimred_thresholded_events.save_figure = params.do_dimred_thresholded_events;
cfg.dimred_thresholded_events.save_dir = fullfile(dimred_te_root, 'mat');
cfg.dimred_thresholded_events.figure_dir = fullfile(dimred_te_root, 'fig');
cfg.dimred_thresholded_events.save_stem = sprintf('%s_dimred_thresholded_events', dataset_name);
cfg.dimred_thresholded_events.save_tag = [condition_tag '_' method_tag];

cfg.progress = struct();
cfg.progress.verbose = true;

if isfield(dataset_cfg, 'sessions')
    cfg.dataset.sessions = dataset_cfg.sessions;
end
end


function tf = local_stage_save_flag(params, field_name, default_value)
tf = logical(default_value);
if isfield(params, field_name) && ~isempty(params.(field_name))
    tf = logical(params.(field_name));
end
end


function target = local_apply_options(target, options)
if isempty(fieldnames(options))
    return;
end

names = fieldnames(options);
for i = 1:numel(names)
    target.(names{i}) = options.(names{i});
end
end


function observable_file = local_observable_file(processed_root, dataset_name, run_info)
observable_tag = char(run_info.observable_mode);
observable_file = fullfile(processed_root, dataset_name, ...
    io_project.get_pipeline_stage_name(1, 'dictionary'), ...
    sprintf('%s_low50_high250_g2_%s_single.mat', dataset_name, observable_tag));

if exist(observable_file, 'file') ~= 2
    error('Observable metadata file was not found:\n  %s', observable_file);
end
end
