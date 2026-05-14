function cfg = build_blp_raw_eigenfunction_threshold_cfg( ...
        dataset_cfg, run_info, params, repo_root, ...
        EDMD_outputs, concat_info, source_info)
%BUILD_BLP_RAW_EIGENFUNCTION_THRESHOLD_CFG Build run-level raw threshold config.

cfg = build_blp_eigenfunction_reduction_defaults();
condition_tag = resolve_blp_eigenfunction_condition_output_tag(run_info, params);
dataset_name = char(run_info.dataset);

raw_td_root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'raw_thresholded_density'), ...
    condition_tag);
raw_te_root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'raw_thresholded_events'), ...
    condition_tag);

observable_file = local_observable_file(params.processed_root, dataset_name, run_info);

cfg.repo_root = repo_root;
cfg.dataset.name = dataset_name;
cfg.dataset.processed_root = params.processed_root;
cfg.output.root = fullfile(io_project.get_pipeline_stage_dir( ...
    params.processed_root, dataset_name, 5, 'eigenfunction_reduction'), condition_tag);
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

cfg.thresholded_density.enable = params.do_thresholded_density;
cfg.thresholded_density.observable_file = observable_file;
cfg.thresholded_density.threshold_mode = params.threshold_mode;
cfg.thresholded_density.threshold_ratio = params.threshold_ratio;
cfg.thresholded_density.threshold_ratio_sweep = params.threshold_ratio_sweep;
cfg.thresholded_density.save_results = local_stage_save_flag( ...
    params, 'save_thresholded_density_results', params.do_thresholded_density);
cfg.thresholded_density.make_figure = params.do_thresholded_density;
cfg.thresholded_density.save_figure = params.do_thresholded_density;
cfg.thresholded_density.save_dir = fullfile(raw_td_root, 'mat');
cfg.thresholded_density.figure_dir = fullfile(raw_td_root, 'fig');
cfg.thresholded_density.save_stem = sprintf('%s_thresholded_density', dataset_name);
cfg.thresholded_density.save_tag = condition_tag;

cfg.thresholded_events.enable = params.do_thresholded_events;
cfg.thresholded_events.observable_file = observable_file;
cfg.thresholded_events.threshold_mode = params.threshold_mode;
cfg.thresholded_events.threshold_ratio = params.threshold_ratio;
cfg.thresholded_events.threshold_ratio_sweep = [];
cfg.thresholded_events.save_results = false;
cfg.thresholded_events.make_figure = params.do_thresholded_events;
cfg.thresholded_events.save_figure = params.do_thresholded_events;
cfg.thresholded_events.save_dir = fullfile(raw_te_root, 'mat');
cfg.thresholded_events.figure_dir = fullfile(raw_te_root, 'fig');
cfg.thresholded_events.save_stem = sprintf('%s_thresholded_events', dataset_name);
cfg.thresholded_events.save_tag = condition_tag;

cfg.dimred_thresholded_density.enable = false;
cfg.dimred_thresholded_events.enable = false;
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


function observable_file = local_observable_file(processed_root, dataset_name, run_info)
observable_tag = char(run_info.observable_mode);
observable_file = fullfile(processed_root, dataset_name, ...
    io_project.get_pipeline_stage_name(1, 'dictionary'), ...
    sprintf('%s_low50_high250_g2_%s_single.mat', dataset_name, observable_tag));

if exist(observable_file, 'file') ~= 2
    error('Observable metadata file was not found:\n  %s', observable_file);
end
end
