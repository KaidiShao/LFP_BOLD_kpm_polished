function [result, prep, B] = run_bold_eigenfunction_reduction_pipeline(cfg)
%RUN_BOLD_EIGENFUNCTION_REDUCTION_PIPELINE Run one P9 BOLD efun DR method.

if nargin < 1 || isempty(cfg)
    error('cfg must be provided.');
end

progress_enabled = local_progress_enabled(cfg);
t_pipeline = tic;

stage_tic = local_stage_start(progress_enabled, 'Preparing BOLD eigenfunction inputs');
[prep, B] = prepare_bold_eigenfunction_reduction_inputs(cfg.source.bold_post_file, cfg);
local_stage_done(progress_enabled, 'Prepared BOLD eigenfunction inputs', stage_tic);

switch lower(cfg.path.kind)
    case 'time'
        method_name = cfg.path.time.method;
        stage_tic = local_stage_start(progress_enabled, ...
            sprintf('Running dimension reduction: %s/%s', cfg.path.kind, method_name));
        reduction = reduce_eigenfunction_time_path( ...
            prep.bold_feature_time_by_mode, cfg.path.time);
        local_stage_done(progress_enabled, 'Finished dimension reduction', stage_tic);

    case 'spectrum'
        spectrum_cfg = cfg.path.spectrum;
        spectrum_cfg.evalues_discrete = prep.evalues_discrete;
        spectrum_cfg.evalues_bilinear = prep.evalues_bilinear;
        spectrum_cfg.feature_variant = prep.feature.variant;
        method_name = cfg.path.spectrum.method;
        stage_tic = local_stage_start(progress_enabled, ...
            sprintf('Running dimension reduction: %s/%s', cfg.path.kind, method_name));
        reduction = reduce_eigenfunction_spectrum_path( ...
            prep.bold_feature_time_by_mode, spectrum_cfg);
        local_stage_done(progress_enabled, 'Finished dimension reduction', stage_tic);

    otherwise
        error('Unknown cfg.path.kind = %s. Use ''time'' or ''spectrum''.', cfg.path.kind);
end

result = build_bold_eigenfunction_reduction_result(prep, reduction, method_name, cfg, B);

if isfield(cfg, 'plot') && isstruct(cfg.plot) && ...
        isfield(cfg.plot, 'enable') && cfg.plot.enable
    stage_tic = local_stage_start(progress_enabled, 'Writing BOLD observable/efun/component summary plot');
    plot_params = cfg.plot;
    if ~isfield(plot_params, 'save_dir') || isempty(plot_params.save_dir)
        plot_params.save_dir = cfg.output.figure_dir;
    end
    if ~isfield(plot_params, 'save_tag') || isempty(plot_params.save_tag)
        plot_params.save_tag = sprintf('%s_%s', cfg.output.feature_tag, cfg.output.method_tag);
    end
    [~, plot_info] = plot_bold_observable_efun_dimred_summary(B, result, plot_params);
    result.artifacts.summary_png_file = local_get_field(plot_info, 'png_file', '');
    result.artifacts.summary_fig_file = local_get_field(plot_info, 'fig_file', '');
    result.plot_info = plot_info;
    local_stage_done(progress_enabled, 'Wrote summary plot', stage_tic);
end

if cfg.save.enable
    stage_tic = local_stage_start(progress_enabled, 'Saving P9 reduction result');
    if exist(cfg.save.dir, 'dir') ~= 7
        mkdir(cfg.save.dir);
    end
    save_path = local_build_save_path(cfg);
    result.artifacts.result_mat_file = save_path;
    save_vars = struct('result', result);
    if cfg.save.v7_3
        save(save_path, '-struct', 'save_vars', 'result', '-v7.3');
    else
        save(save_path, '-struct', 'save_vars', 'result');
    end
    local_stage_done(progress_enabled, 'Saved P9 reduction result', stage_tic);
end

local_stage_done(progress_enabled, 'Pipeline 9 method finished', t_pipeline);
end


function tf = local_progress_enabled(cfg)
tf = true;
if isfield(cfg, 'progress') && isfield(cfg.progress, 'verbose') && ...
        ~isempty(cfg.progress.verbose)
    tf = logical(cfg.progress.verbose);
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
    fprintf('[stage] %s in %.1fs.\n', label, toc(stage_tic));
end
end


function save_path = local_build_save_path(cfg)
stem = cfg.save.file_stem;
tag = cfg.save.tag;
if isempty(tag)
    file_name = sprintf('%s.mat', stem);
else
    file_name = sprintf('%s_%s.mat', stem, tag);
end
save_path = fullfile(cfg.save.dir, file_name);
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
