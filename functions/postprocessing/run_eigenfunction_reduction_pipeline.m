function [result, EDMD_outputs, concat_info, source_info] = run_eigenfunction_reduction_pipeline(cfg)
%RUN_EIGENFUNCTION_REDUCTION_PIPELINE Stage orchestrator for one eigenfunction reduction run.

if nargin < 1 || isempty(cfg)
    error('cfg must be provided.');
end

cfg = apply_blp_eigenfunction_reduction_defaults(cfg);
progress_enabled = local_progress_enabled(cfg);
pipeline_tic = tic;

stage_tic = local_stage_start(progress_enabled, 'Loading EDMD source');
[EDMD_outputs, concat_info, source_info] = io_edmd.load_edmd_source(cfg.source);
local_stage_done(progress_enabled, 'Loaded EDMD source', stage_tic);

stage_tic = local_stage_start(progress_enabled, 'Preparing eigenfunction inputs');
prep = prepare_blp_eigenfunction_reduction_inputs(EDMD_outputs, cfg);
local_stage_done(progress_enabled, 'Prepared eigenfunction inputs', stage_tic);
local_print_prepared_input_summary(progress_enabled, prep, EDMD_outputs);

[reduction_label, method_name] = local_reduction_stage_label(cfg);
stage_tic = local_stage_start(progress_enabled, reduction_label);
[reduction, method_name] = run_blp_eigenfunction_reduction_method(prep, cfg);
local_stage_done(progress_enabled, 'Finished dimension reduction', stage_tic);

stage_tic = local_stage_start(progress_enabled, 'Packing reduction result');
result = build_blp_eigenfunction_reduction_result( ...
    prep, reduction, method_name, cfg, source_info, concat_info);
local_stage_done(progress_enabled, 'Packed reduction result', stage_tic);

if cfg.thresholded_density.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded eigenfunction density');
    [result.thresholded_density, td_artifacts] = run_blp_thresholded_density_stage( ...
        prep, result, cfg.thresholded_density);
    result.artifacts.thresholded_density_mat_file = td_artifacts.mat_file;
    result.artifacts.thresholded_density_figure_file = td_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded eigenfunction density', stage_tic);
end

if cfg.thresholded_events.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded eigenfunction events');
    [result.thresholded_events, te_artifacts] = run_blp_thresholded_events_stage( ...
        prep, result, cfg.thresholded_events);
    result.artifacts.thresholded_events_mat_file = te_artifacts.mat_file;
    result.artifacts.thresholded_events_figure_file = te_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded eigenfunction events', stage_tic);
end

if cfg.dimred_thresholded_density.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded component density');
    [result.dimred_thresholded_density, dtd_artifacts] = ...
        run_blp_dimred_thresholded_density_stage( ...
        prep, result, cfg.dimred_thresholded_density);
    result.artifacts.dimred_thresholded_density_mat_file = dtd_artifacts.mat_file;
    result.artifacts.dimred_thresholded_density_figure_file = dtd_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded component density', stage_tic);
end

if cfg.dimred_thresholded_events.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded component events');
    [result.dimred_thresholded_events, dte_artifacts] = ...
        run_blp_dimred_thresholded_events_stage( ...
        prep, result, cfg.dimred_thresholded_events);
    result.artifacts.dimred_thresholded_events_mat_file = dte_artifacts.mat_file;
    result.artifacts.dimred_thresholded_events_figure_file = dte_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded component events', stage_tic);
end

if cfg.save.enable
    stage_tic = local_stage_start(progress_enabled, 'Saving compact reduction result');
    result = save_blp_eigenfunction_reduction_result(result, cfg);
    local_stage_done(progress_enabled, 'Saved reduction result', stage_tic);
end

local_stage_done(progress_enabled, 'Pipeline method finished', pipeline_tic);
end


function tf = local_progress_enabled(cfg)
tf = true;
if isfield(cfg, 'progress') && isfield(cfg.progress, 'verbose') && ...
        ~isempty(cfg.progress.verbose)
    tf = logical(cfg.progress.verbose);
elseif isfield(cfg, 'source') && isfield(cfg.source, 'concat') && ...
        isfield(cfg.source.concat, 'verbose') && ~isempty(cfg.source.concat.verbose)
    tf = logical(cfg.source.concat.verbose);
end
end


function stage_tic = local_stage_start(enabled, label)
stage_tic = tic;
if enabled
    fprintf('[stage] %s...\n', label);
end
end


function local_stage_done(enabled, label, stage_tic)
if enabled
    fprintf('[stage] %s in %s.\n', label, local_format_seconds(toc(stage_tic)));
end
end


function local_print_prepared_input_summary(enabled, prep, EDMD_outputs)
if ~enabled
    return;
end

raw_size = size(prep.efun_raw_time_by_mode);
feature_size = size(prep.efun_feature_time_by_mode);
n_original_modes = numel(EDMD_outputs.evalues);
n_selected_modes = numel(prep.evalues_discrete);

if isempty(prep.dt)
    dt_text = 'missing';
    duration_text = 'unknown';
else
    dt_text = sprintf('%.10g sec', prep.dt);
    duration_text = local_format_duration_for_summary( ...
        (raw_size(1) - 1) * double(prep.dt));
end

fprintf(['[summary] Eigenfunction length: T=%d samples | selected modes=%d/%d | ', ...
    'raw=[%s] | feature=[%s]\n'], ...
    raw_size(1), n_selected_modes, n_original_modes, ...
    num2str(raw_size), num2str(feature_size));
fprintf('[summary] Time axis: dt=%s | duration=%s | dt_source=%s.%s\n', ...
    dt_text, duration_text, prep.dt_source.source, prep.dt_source.field);
end


function txt = local_format_duration_for_summary(seconds_in)
if ~isfinite(seconds_in)
    txt = 'unknown';
    return;
end

seconds_in = max(0, seconds_in);
days = floor(seconds_in / 86400);
hours = floor(mod(seconds_in, 86400) / 3600);
minutes = floor(mod(seconds_in, 3600) / 60);
seconds = floor(mod(seconds_in, 60));

if days > 0
    txt = sprintf('%dd %02d:%02d:%02d', days, hours, minutes, seconds);
elseif hours > 0
    txt = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
else
    txt = sprintf('%02d:%02d', minutes, seconds);
end
end


function [label, method_name] = local_reduction_stage_label(cfg)
switch lower(cfg.path.kind)
    case 'time'
        method_name = cfg.path.time.method;
    case 'spectrum'
        method_name = cfg.path.spectrum.method;
    otherwise
        error('Unknown cfg.path.kind = %s. Use ''time'' or ''spectrum''.', cfg.path.kind);
end

label = sprintf('Running dimension reduction: %s/%s', cfg.path.kind, method_name);
end


function txt = local_format_seconds(seconds_in)
if ~isfinite(seconds_in)
    txt = '--';
    return;
end

seconds_in = max(0, seconds_in);
if seconds_in < 1
    txt = sprintf('%.3fs', seconds_in);
elseif seconds_in < 60
    txt = sprintf('%.1fs', seconds_in);
else
    hours = floor(seconds_in / 3600);
    minutes = floor(mod(seconds_in, 3600) / 60);
    seconds = floor(mod(seconds_in, 60));
    if hours > 0
        txt = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
    else
        txt = sprintf('%02d:%02d', minutes, seconds);
    end
end
end

% Legacy in-file defaults helpers remain parked below.
% The canonical pipeline 5 path now uses apply_blp_eigenfunction_reduction_defaults.


function cfg = local_apply_defaults(cfg)
if ~isfield(cfg, 'feature') || ~isfield(cfg.feature, 'variant')
    error('cfg.feature.variant must be provided.');
end

if ~isfield(cfg, 'feature') || ~isfield(cfg.feature, 'normalization') || isempty(cfg.feature.normalization)
    cfg.feature.normalization = 'maxabs_per_mode';
end

if ~isfield(cfg, 'input') || ~isfield(cfg.input, 'dt')
    cfg.input.dt = [];
end

if ~isfield(cfg.input, 'observable_file')
    cfg.input.observable_file = '';
end

if ~isfield(cfg, 'selection') || ~isfield(cfg.selection, 'abs_thresh')
    cfg.selection.abs_thresh = 0.01;
end

if ~isfield(cfg.selection, 'sort_by') || isempty(cfg.selection.sort_by)
    cfg.selection.sort_by = 'modulus';
end

if ~isfield(cfg.selection, 'sort_dir') || isempty(cfg.selection.sort_dir)
    cfg.selection.sort_dir = 'descend';
end

if ~isfield(cfg.selection, 'max_modes') || isempty(cfg.selection.max_modes)
    cfg.selection.max_modes = Inf;
end

if ~isfield(cfg, 'summary') || ~isfield(cfg.summary, 'smooth')
    cfg.summary.smooth = struct();
end

if ~isfield(cfg.summary.smooth, 'enable') || isempty(cfg.summary.smooth.enable)
    cfg.summary.smooth.enable = false;
end

if ~isfield(cfg.summary.smooth, 'method') || isempty(cfg.summary.smooth.method)
    cfg.summary.smooth.method = 'movmean';
end

if ~isfield(cfg.summary.smooth, 'window') || isempty(cfg.summary.smooth.window)
    cfg.summary.smooth.window = 10;
end

if ~isfield(cfg, 'save') || ~isfield(cfg.save, 'enable')
    cfg.save.enable = false;
end

if ~isfield(cfg.save, 'dir') || isempty(cfg.save.dir)
    cfg.save.dir = local_default_save_dir(cfg);
end

if ~isfield(cfg.save, 'file_stem') || isempty(cfg.save.file_stem)
    cfg.save.file_stem = local_default_file_stem(cfg);
end

if ~isfield(cfg.save, 'tag')
    cfg.save.tag = '';
end

if ~isfield(cfg.save, 'payload') || isempty(cfg.save.payload)
    cfg.save.payload = 'compact';
end

if ~isfield(cfg.save, 'v7_3') || isempty(cfg.save.v7_3)
    cfg.save.v7_3 = true;
end

cfg = local_apply_thresholded_density_defaults(cfg);
cfg = local_apply_thresholded_events_defaults(cfg);
cfg = local_apply_dimred_thresholded_density_defaults(cfg);
cfg = local_apply_dimred_thresholded_events_defaults(cfg);
end


function save_dir = local_default_save_dir(cfg)
source_run_name = local_source_run_name(cfg);
dataset_name = local_dataset_name(cfg);
stage_root = io_project.get_pipeline_stage_dir( ...
    io_project.get_project_processed_root(), dataset_name, 5, 'eigenfunction_reduction');
save_dir = fullfile(stage_root, source_run_name, 'mat');
end


function file_stem = local_default_file_stem(cfg)
if isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    file_stem = sprintf('%s_efun', cfg.dataset.name);
else
    file_stem = 'efun';
end
end


function source_run_name = local_source_run_name(cfg)
source_run_name = 'unspecified_source';

if isfield(cfg, 'output') && isfield(cfg.output, 'output_run_name') && ...
        ~isempty(cfg.output.output_run_name)
    source_run_name = cfg.output.output_run_name;
elseif isfield(cfg, 'output') && isfield(cfg.output, 'source_run_name') && ...
        ~isempty(cfg.output.source_run_name)
    source_run_name = cfg.output.source_run_name;
elseif isfield(cfg, 'source') && isfield(cfg.source, 'data_dir') && ...
        ~isempty(cfg.source.data_dir)
    [~, source_run_name] = fileparts(cfg.source.data_dir);
elseif isfield(cfg, 'source') && isfield(cfg.source, 'edmd_file') && ...
        ~isempty(cfg.source.edmd_file)
    [~, source_run_name] = fileparts(cfg.source.edmd_file);
end

source_run_name = regexprep(char(source_run_name), '[^\w\-]+', '_');
end


function cfg = local_apply_thresholded_density_defaults(cfg)
if ~isfield(cfg, 'thresholded_density') || ~isstruct(cfg.thresholded_density)
    cfg.thresholded_density = struct();
end

td = cfg.thresholded_density;

if ~isfield(td, 'enable') || isempty(td.enable)
    td.enable = false;
end

if ~isfield(td, 'window_sec') || isempty(td.window_sec)
    td.window_sec = 2;
end

if ~isfield(td, 'step_sec') || isempty(td.step_sec)
    td.step_sec = td.window_sec;
end

if ~isfield(td, 'threshold_ratio') || isempty(td.threshold_ratio)
    td.threshold_ratio = 0.7;
end

if ~isfield(td, 'threshold_mode') || isempty(td.threshold_mode)
    td.threshold_mode = 'meanplusstd';
end

if ~isfield(td, 'value_transform') || isempty(td.value_transform)
    td.value_transform = 'none';
end

if ~isfield(td, 'density_denominator') || isempty(td.density_denominator)
    td.density_denominator = 'window_samples';
end

if ~isfield(td, 'output_class') || isempty(td.output_class)
    td.output_class = 'single';
end

if ~isfield(td, 'smooth_density') || isempty(td.smooth_density)
    td.smooth_density = false;
end

if ~isfield(td, 'smooth_window_sec') || isempty(td.smooth_window_sec)
    td.smooth_window_sec = 0;
end

if ~isfield(td, 'require_session_metadata') || isempty(td.require_session_metadata)
    td.require_session_metadata = false;
end

if ~isfield(td, 'observable_file') || isempty(td.observable_file)
    if isfield(cfg, 'input') && isfield(cfg.input, 'observable_file')
        td.observable_file = cfg.input.observable_file;
    else
        td.observable_file = '';
    end
end

if ~isfield(td, 'save_results') || isempty(td.save_results)
    td.save_results = td.enable;
end

if ~isfield(td, 'make_figure') || isempty(td.make_figure)
    td.make_figure = td.enable;
end

if ~isfield(td, 'save_figure') || isempty(td.save_figure)
    td.save_figure = td.enable;
end

if ~isfield(td, 'close_figure') || isempty(td.close_figure)
    td.close_figure = true;
end

if ~isfield(td, 'save_v7_3') || isempty(td.save_v7_3)
    td.save_v7_3 = true;
end

output_root = local_output_root(cfg);
if ~isfield(td, 'save_dir') || isempty(td.save_dir)
    td.save_dir = fullfile(output_root, 'thresholded_density', 'mat');
end

if ~isfield(td, 'figure_dir') || isempty(td.figure_dir)
    td.figure_dir = fullfile(output_root, 'thresholded_density', 'fig');
end

if ~isfield(td, 'save_stem') || isempty(td.save_stem)
    td.save_stem = sprintf('%s_thresholded_density', local_dataset_name(cfg));
end

if ~isfield(td, 'save_tag')
    td.save_tag = '';
end

if ~isfield(td, 'figure_visible') || isempty(td.figure_visible)
    td.figure_visible = 'off';
end

if ~isfield(td, 'figure_resolution') || isempty(td.figure_resolution)
    td.figure_resolution = 200;
end

if ~isfield(td, 'max_plot_modes') || isempty(td.max_plot_modes)
    td.max_plot_modes = 100;
end

if ~isfield(td, 'title') || isempty(td.title)
    td.title = 'Thresholded Eigenfunction Density';
end

if ~isfield(td, 'colormap') || isempty(td.colormap)
    td.colormap = 'turbo';
end

if ~isfield(td, 'background_color') || isempty(td.background_color)
    td.background_color = [0, 0, 0];
end

if ~isfield(td, 'axes_color') || isempty(td.axes_color)
    td.axes_color = td.background_color;
end

if ~isfield(td, 'text_color') || isempty(td.text_color)
    td.text_color = [1, 1, 1];
end

if ~isfield(td, 'grid_color') || isempty(td.grid_color)
    td.grid_color = [0.75, 0.75, 0.75];
end

cfg.thresholded_density = td;
end


function cfg = local_apply_thresholded_events_defaults(cfg)
if ~isfield(cfg, 'thresholded_events') || ~isstruct(cfg.thresholded_events)
    cfg.thresholded_events = struct();
end

te = cfg.thresholded_events;
td = cfg.thresholded_density;

if ~isfield(te, 'enable') || isempty(te.enable)
    te.enable = false;
end

if ~isfield(te, 'threshold_ratio') || isempty(te.threshold_ratio)
    te.threshold_ratio = td.threshold_ratio;
end

if ~isfield(te, 'threshold_mode') || isempty(te.threshold_mode)
    te.threshold_mode = td.threshold_mode;
end

if ~isfield(te, 'value_transform') || isempty(te.value_transform)
    te.value_transform = td.value_transform;
end

if ~isfield(te, 'event_detector') || isempty(te.event_detector)
    te.event_detector = 'find_peak_loc';
end

if ~isfield(te, 'event_score') || isempty(te.event_score)
    te.event_score = 'area_above_threshold';
end

if ~isfield(te, 'min_duration_sec') || isempty(te.min_duration_sec)
    te.min_duration_sec = 0.03;
end

if ~isfield(te, 'merge_gap_sec') || isempty(te.merge_gap_sec)
    te.merge_gap_sec = 0.02;
end

if ~isfield(te, 'find_peak_loc_window_sec')
    te.find_peak_loc_window_sec = [];
end

if ~isfield(te, 'find_peak_loc_window_samples')
    te.find_peak_loc_window_samples = [];
end

if ~isfield(te, 'find_peak_loc_drop_first') || isempty(te.find_peak_loc_drop_first)
    te.find_peak_loc_drop_first = false;
end

if ~isfield(te, 'find_peak_loc_post_nms') || isempty(te.find_peak_loc_post_nms)
    te.find_peak_loc_post_nms = false;
end

if ~isfield(te, 'bin_sec') || isempty(te.bin_sec)
    te.bin_sec = td.window_sec;
end

if ~isfield(te, 'step_sec') || isempty(te.step_sec)
    te.step_sec = te.bin_sec;
end

if ~isfield(te, 'smooth_rate') || isempty(te.smooth_rate)
    te.smooth_rate = true;
end

if ~isfield(te, 'smooth_sigma_sec') || isempty(te.smooth_sigma_sec)
    te.smooth_sigma_sec = te.bin_sec;
end

if ~isfield(te, 'output_class') || isempty(te.output_class)
    te.output_class = 'single';
end

if ~isfield(te, 'require_session_metadata') || isempty(te.require_session_metadata)
    te.require_session_metadata = true;
end

if ~isfield(te, 'observable_file') || isempty(te.observable_file)
    if isfield(cfg, 'input') && isfield(cfg.input, 'observable_file')
        te.observable_file = cfg.input.observable_file;
    else
        te.observable_file = '';
    end
end

if ~isfield(te, 'nms') || ~isstruct(te.nms)
    te.nms = struct();
end
if ~isfield(te.nms, 'enable') || isempty(te.nms.enable)
    te.nms.enable = true;
end
if ~isfield(te.nms, 'mode') || isempty(te.nms.mode)
    te.nms.mode = 'eigen_period';
end
if ~isfield(te.nms, 'period_fraction') || isempty(te.nms.period_fraction)
    te.nms.period_fraction = 1.0;
end
if ~isfield(te.nms, 'fixed_sec') || isempty(te.nms.fixed_sec)
    te.nms.fixed_sec = 0.25;
end
if ~isfield(te.nms, 'default_sec') || isempty(te.nms.default_sec)
    te.nms.default_sec = 0.25;
end
if ~isfield(te.nms, 'min_sec') || isempty(te.nms.min_sec)
    te.nms.min_sec = 0.05;
end
if ~isfield(te.nms, 'max_sec') || isempty(te.nms.max_sec)
    te.nms.max_sec = 2.0;
end
if ~isfield(te.nms, 'angle_eps') || isempty(te.nms.angle_eps)
    te.nms.angle_eps = 1e-4;
end

if ~isfield(te, 'save_results') || isempty(te.save_results)
    te.save_results = te.enable;
end

if ~isfield(te, 'make_figure') || isempty(te.make_figure)
    te.make_figure = te.enable;
end

if ~isfield(te, 'save_figure') || isempty(te.save_figure)
    te.save_figure = te.enable;
end

if ~isfield(te, 'close_figure') || isempty(te.close_figure)
    te.close_figure = true;
end

if ~isfield(te, 'save_v7_3') || isempty(te.save_v7_3)
    te.save_v7_3 = true;
end

output_root = local_output_root(cfg);
if ~isfield(te, 'save_dir') || isempty(te.save_dir)
    te.save_dir = fullfile(output_root, 'thresholded_events', 'mat');
end

if ~isfield(te, 'figure_dir') || isempty(te.figure_dir)
    te.figure_dir = fullfile(output_root, 'thresholded_events', 'fig');
end

if ~isfield(te, 'save_stem') || isempty(te.save_stem)
    te.save_stem = sprintf('%s_thresholded_events', local_dataset_name(cfg));
end

if ~isfield(te, 'save_tag')
    te.save_tag = '';
end

if ~isfield(te, 'figure_visible') || isempty(te.figure_visible)
    te.figure_visible = 'off';
end

if ~isfield(te, 'figure_resolution') || isempty(te.figure_resolution)
    te.figure_resolution = 200;
end

if ~isfield(te, 'max_plot_modes') || isempty(te.max_plot_modes)
    te.max_plot_modes = 80;
end

if ~isfield(te, 'max_plot_events') || isempty(te.max_plot_events)
    te.max_plot_events = 5000;
end

if ~isfield(te, 'title') || isempty(te.title)
    te.title = 'Thresholded Eigenfunction Events';
end

if ~isfield(te, 'colormap') || isempty(te.colormap)
    te.colormap = 'turbo';
end

if ~isfield(te, 'background_color') || isempty(te.background_color)
    te.background_color = [0, 0, 0];
end

if ~isfield(te, 'axes_color') || isempty(te.axes_color)
    te.axes_color = te.background_color;
end

if ~isfield(te, 'text_color') || isempty(te.text_color)
    te.text_color = [1, 1, 1];
end

if ~isfield(te, 'grid_color') || isempty(te.grid_color)
    te.grid_color = [0.75, 0.75, 0.75];
end

cfg.thresholded_events = te;
end


function cfg = local_apply_dimred_thresholded_density_defaults(cfg)
if ~isfield(cfg, 'dimred_thresholded_density') || ...
        ~isstruct(cfg.dimred_thresholded_density)
    cfg.dimred_thresholded_density = struct();
end

dtd = cfg.dimred_thresholded_density;
td = cfg.thresholded_density;

if ~isfield(dtd, 'enable') || isempty(dtd.enable)
    dtd.enable = false;
end

dtd = local_copy_missing_fields(dtd, td, { ...
    'window_sec', 'step_sec', 'window_samples', 'step_samples', ...
    'threshold_ratio', 'threshold_mode', 'value_transform', ...
    'density_denominator', 'output_class', 'smooth_density', ...
    'smooth_window_sec', 'smooth_window_bins', 'save_v7_3', ...
    'require_session_metadata', 'observable_file', ...
    'session_start_idx', 'session_end_idx', 'session_lengths', ...
    'session_ids', 'session_dx', ...
    'figure_visible', 'figure_resolution', 'colormap', ...
    'background_color', 'axes_color', 'text_color', 'grid_color'});

if ~isfield(dtd, 'save_results') || isempty(dtd.save_results)
    dtd.save_results = dtd.enable;
end
if ~isfield(dtd, 'make_figure') || isempty(dtd.make_figure)
    dtd.make_figure = dtd.enable;
end
if ~isfield(dtd, 'save_figure') || isempty(dtd.save_figure)
    dtd.save_figure = dtd.enable;
end
if ~isfield(dtd, 'close_figure') || isempty(dtd.close_figure)
    dtd.close_figure = true;
end

output_root = local_output_root(cfg);
if ~isfield(dtd, 'save_dir') || isempty(dtd.save_dir)
    dtd.save_dir = fullfile(output_root, 'dimred_thresholded_density', 'mat');
end
if ~isfield(dtd, 'figure_dir') || isempty(dtd.figure_dir)
    dtd.figure_dir = fullfile(output_root, 'dimred_thresholded_density', 'fig');
end
if ~isfield(dtd, 'save_stem') || isempty(dtd.save_stem)
    dtd.save_stem = sprintf('%s_dimred_thresholded_density', ...
        local_dataset_name(cfg));
end
if ~isfield(dtd, 'save_tag')
    dtd.save_tag = '';
end
if ~isfield(dtd, 'max_plot_modes') || isempty(dtd.max_plot_modes)
    dtd.max_plot_modes = 40;
end
if ~isfield(dtd, 'title') || isempty(dtd.title)
    dtd.title = 'Thresholded Eigenfunction Component Density';
end

cfg.dimred_thresholded_density = dtd;
end


function cfg = local_apply_dimred_thresholded_events_defaults(cfg)
if ~isfield(cfg, 'dimred_thresholded_events') || ...
        ~isstruct(cfg.dimred_thresholded_events)
    cfg.dimred_thresholded_events = struct();
end

dte = cfg.dimred_thresholded_events;
te = cfg.thresholded_events;

if ~isfield(dte, 'enable') || isempty(dte.enable)
    dte.enable = false;
end

dte = local_copy_missing_fields(dte, te, { ...
    'threshold_ratio', 'threshold_mode', 'value_transform', ...
    'event_detector', 'event_score', 'min_duration_sec', 'merge_gap_sec', ...
    'min_duration_samples', 'merge_gap_samples', 'bin_sec', ...
    'step_sec', 'bin_samples', 'step_samples', 'smooth_rate', ...
    'smooth_sigma_sec', 'output_class', 'require_session_metadata', ...
    'find_peak_loc_window_sec', 'find_peak_loc_window_samples', ...
    'find_peak_loc_drop_first', 'find_peak_loc_post_nms', ...
    'observable_file', 'save_v7_3', 'figure_visible', ...
    'figure_resolution', 'max_plot_events', 'colormap', ...
    'background_color', 'axes_color', 'text_color', 'grid_color'});

if ~isfield(dte, 'component_evalue_weight_transform') || ...
        isempty(dte.component_evalue_weight_transform)
    dte.component_evalue_weight_transform = 'abs';
end
if ~isfield(dte, 'component_evalue_angle_statistic') || ...
        isempty(dte.component_evalue_angle_statistic)
    dte.component_evalue_angle_statistic = 'weighted_median_abs';
end
if ~isfield(dte, 'component_evalue_top_modes') || ...
        isempty(dte.component_evalue_top_modes)
    dte.component_evalue_top_modes = 5;
end

if ~isfield(dte, 'nms') || ~isstruct(dte.nms)
    dte.nms = struct();
end
if isfield(te, 'nms') && isstruct(te.nms)
    dte.nms = local_copy_missing_fields(dte.nms, te.nms, { ...
        'enable', 'mode', 'period_fraction', 'fixed_sec', ...
        'default_sec', 'min_sec', 'max_sec', 'angle_eps'});
end
if ~isfield(dte.nms, 'enable') || isempty(dte.nms.enable)
    dte.nms.enable = true;
end
if ~isfield(dte.nms, 'mode') || isempty(dte.nms.mode)
    dte.nms.mode = 'eigen_period';
end
if ~isfield(dte.nms, 'period_fraction') || isempty(dte.nms.period_fraction)
    dte.nms.period_fraction = 1.0;
end
if ~isfield(dte.nms, 'fixed_sec') || isempty(dte.nms.fixed_sec)
    dte.nms.fixed_sec = 0.25;
end
if ~isfield(dte.nms, 'default_sec') || isempty(dte.nms.default_sec)
    dte.nms.default_sec = 0.25;
end
if ~isfield(dte.nms, 'min_sec') || isempty(dte.nms.min_sec)
    dte.nms.min_sec = 0.05;
end
if ~isfield(dte.nms, 'max_sec') || isempty(dte.nms.max_sec)
    dte.nms.max_sec = 2.0;
end
if ~isfield(dte.nms, 'angle_eps') || isempty(dte.nms.angle_eps)
    dte.nms.angle_eps = 1e-4;
end

if ~isfield(dte, 'save_results') || isempty(dte.save_results)
    dte.save_results = dte.enable;
end
if ~isfield(dte, 'make_figure') || isempty(dte.make_figure)
    dte.make_figure = dte.enable;
end
if ~isfield(dte, 'save_figure') || isempty(dte.save_figure)
    dte.save_figure = dte.enable;
end
if ~isfield(dte, 'close_figure') || isempty(dte.close_figure)
    dte.close_figure = true;
end

output_root = local_output_root(cfg);
if ~isfield(dte, 'save_dir') || isempty(dte.save_dir)
    dte.save_dir = fullfile(output_root, 'dimred_thresholded_events', 'mat');
end
if ~isfield(dte, 'figure_dir') || isempty(dte.figure_dir)
    dte.figure_dir = fullfile(output_root, 'dimred_thresholded_events', 'fig');
end
if ~isfield(dte, 'save_stem') || isempty(dte.save_stem)
    dte.save_stem = sprintf('%s_dimred_thresholded_events', ...
        local_dataset_name(cfg));
end
if ~isfield(dte, 'save_tag')
    dte.save_tag = '';
end
if ~isfield(dte, 'max_plot_modes') || isempty(dte.max_plot_modes)
    dte.max_plot_modes = 40;
end
if ~isfield(dte, 'title') || isempty(dte.title)
    dte.title = 'Thresholded Eigenfunction Component Events';
end

cfg.dimred_thresholded_events = dte;
end


function dst = local_copy_missing_fields(dst, src, names)
for i = 1:numel(names)
    name = names{i};
    if (~isfield(dst, name) || isempty(dst.(name))) && isfield(src, name)
        dst.(name) = src.(name);
    end
end
end


function output_root = local_output_root(cfg)
if isfield(cfg, 'output') && isfield(cfg.output, 'root') && ...
        ~isempty(cfg.output.root)
    output_root = cfg.output.root;
elseif isfield(cfg, 'save') && isfield(cfg.save, 'dir') && ...
        ~isempty(cfg.save.dir)
    output_root = fileparts(cfg.save.dir);
else
    output_root = pwd;
end
end


function dataset_name = local_dataset_name(cfg)
dataset_name = 'efun';
if isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    dataset_name = char(cfg.dataset.name);
end
dataset_name = regexprep(dataset_name, '[^\w\-]+', '_');
end
