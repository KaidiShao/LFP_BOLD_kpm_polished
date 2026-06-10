function stage = run_bold_reskoopnet_postprocessing_stages(ctx, params)
%RUN_BOLD_RESKOOPNET_POSTPROCESSING_STAGES Execute the shared pipeline 7 stage sequence.

if nargin < 1 || ~isstruct(ctx)
    error('ctx must be a struct.');
end
if nargin < 2 || ~isstruct(params)
    error('params must be a struct.');
end
if ~isfield(ctx, 'skip_existing') || ctx.skip_existing
    error('ctx must be a prepared non-skipped run context.');
end

t_run = tic;

stage = struct();
stage.ctx = ctx;
stage.run_info = ctx.run_info;
stage.EDMD_outputs = ctx.EDMD_outputs;
stage.concat_info = ctx.concat_info;
stage.source_info = ctx.source_info;
stage.observable_file = ctx.observable_file;
stage.obs_meta = ctx.obs_meta;
stage.session = ctx.session;
stage.dt = ctx.dt;
stage.time_vec = ctx.time_vec;
stage.session_border_t = ctx.session_border_t;

t_stage = tic;
local_log_stage(params, 'main postprocess started');
stage.post_opts = params.post_opts;
stage.post_opts.dt = ctx.dt;
stage.post_opts.do_plot = params.make_main_plot;
stage.post_opts.time_vec = ctx.time_vec;
stage.post_opts.session_border = ctx.session_border_t;
stage.post_opts.draw_border = params.draw_session_borders;
[stage.EDMD_outputs, stage.fig_main] = ...
    postprocess_EDMD_outputs(stage.EDMD_outputs, stage.post_opts);
local_log_stage(params, 'main postprocess finished in %.1f sec', toc(t_stage));
stage.main_png = '';
if params.make_main_plot
    t_stage = tic;
    local_log_stage(params, 'main figure save started');
    stage.main_png = save_bold_reskoopnet_postprocessing_figure( ...
        stage.fig_main, ctx.fig_dir, [ctx.run_info.run_name, '_efuns.png'], params);
    local_log_stage(params, 'main figure save finished in %.1f sec', toc(t_stage));
end

stage.deconv_cfg = params.deconv;
stage.fig_deconv = [];
stage.deconv = struct();
stage.deconv_png = '';
if params.compute_deconv || params.make_deconv_plot
    t_stage = tic;
    local_log_stage(params, 'deconv stage started');
    stage.deconv_cfg.dt = ctx.dt;
    stage.deconv_cfg.do_plot = params.make_deconv_plot;
    stage.deconv_cfg.time_vec = ctx.time_vec;
    stage.deconv_cfg.session_border = ctx.session_border_t;
    stage.deconv_cfg.draw_border = params.draw_session_borders;
    stage.deconv_cfg.save_path = '';
    [stage.EDMD_outputs, stage.fig_deconv, stage.deconv] = ...
        postprocess_EDMD_outputs_deconv_efuns(stage.EDMD_outputs, stage.deconv_cfg);
    local_log_stage(params, 'deconv stage finished in %.1f sec', toc(t_stage));
    if params.make_deconv_plot
        t_stage = tic;
        local_log_stage(params, 'deconv figure save started');
        stage.deconv_png = save_bold_reskoopnet_postprocessing_figure( ...
            stage.fig_deconv, ctx.fig_dir, [ctx.run_info.run_name, '_deconv_efuns.png'], params);
        local_log_stage(params, 'deconv figure save finished in %.1f sec', toc(t_stage));
    end
end

stage.timescale_cfg = local_build_timescale_cfg(params, stage.EDMD_outputs, ctx);
stage.timescale_png = '';
if params.make_timescale_plot
    t_stage = tic;
    local_log_stage(params, 'timescale stage started');
    [stage.fig_timescale, stage.timescale_info] = ...
        postprocess_EDMD_outputs_timescale(stage.EDMD_outputs, stage.timescale_cfg);
    local_log_stage(params, 'timescale stage finished in %.1f sec', toc(t_stage));
    t_stage = tic;
    local_log_stage(params, 'timescale figure save started');
    stage.timescale_png = save_bold_reskoopnet_postprocessing_figure( ...
        stage.fig_timescale, ctx.fig_dir, [ctx.run_info.run_name, '_timescale.png'], params);
    local_log_stage(params, 'timescale figure save finished in %.1f sec', toc(t_stage));
else
    stage.fig_timescale = [];
    stage.timescale_info = struct();
end

stage.artifacts = struct();
stage.artifacts.main_png = stage.main_png;
stage.artifacts.deconv_png = stage.deconv_png;
stage.artifacts.timescale_png = stage.timescale_png;
t_stage = tic;
local_log_stage(params, 'BOLD_POST save started');
[stage.BOLD_POST, stage.save_result] = save_bold_reskoopnet_postprocessing_output( ...
    ctx, stage.EDMD_outputs, stage.deconv, stage.timescale_info, ...
    params, stage.artifacts, NaN);
local_log_stage(params, 'BOLD_POST save finished in %.1f sec', toc(t_stage));

stage.act_result = struct();
stage.intrinsic_activation_dir = '';
stage.n_intrinsic_maps = 0;
if isfield(params, 'intrinsic_activation') && isstruct(params.intrinsic_activation) && ...
        isfield(params.intrinsic_activation, 'enabled') && params.intrinsic_activation.enabled
    t_stage = tic;
    local_log_stage(params, 'intrinsic activation maps started');
    act_params = params.intrinsic_activation;
    act_params.processed_root = params.processed_root;
    act_params.datapons_root = params.datapons_root;
    act_params.output_root = ctx.out_root;
    act_params.skip_existing = params.skip_existing && act_params.skip_existing;
    act_params.update_bold_post = false;
    [stage.act_result, stage.BOLD_POST] = export_bold_intrinsic_activation_maps( ...
        stage.BOLD_POST, act_params);
    BOLD_POST = stage.BOLD_POST;
    save_bold_post_mat(ctx.post_file, BOLD_POST);
    stage.intrinsic_activation_dir = stage.act_result.activation_dir;
    stage.n_intrinsic_maps = stage.act_result.n_requested_maps;
    local_log_stage(params, 'intrinsic activation maps finished in %.1f sec', toc(t_stage));
end

stage.roi_summary_result = struct();
stage.intrinsic_roi_summary_dir = '';
stage.intrinsic_roi_summary_png = '';
stage.n_intrinsic_roi_summary_modes = 0;
if isfield(params, 'intrinsic_roi_summary') && isstruct(params.intrinsic_roi_summary) && ...
        isfield(params.intrinsic_roi_summary, 'enabled') && params.intrinsic_roi_summary.enabled
    t_stage = tic;
    local_log_stage(params, 'intrinsic ROI summary started');
    roi_params = params.intrinsic_roi_summary;
    roi_params.processed_root = params.processed_root;
    roi_params.datapons_root = params.datapons_root;
    roi_params.output_root = ctx.out_root;
    roi_params.skip_existing = params.skip_existing && roi_params.skip_existing;
    roi_params.update_bold_post = false;
    [stage.roi_summary_result, stage.BOLD_POST] = export_bold_intrinsic_roi_bar_summaries( ...
        stage.BOLD_POST, roi_params);
    BOLD_POST = stage.BOLD_POST;
    save_bold_post_mat(ctx.post_file, BOLD_POST);
    stage.intrinsic_roi_summary_dir = stage.roi_summary_result.summary_dir;
    stage.intrinsic_roi_summary_png = stage.roi_summary_result.png_file;
    stage.n_intrinsic_roi_summary_modes = stage.roi_summary_result.n_selected_modes;
    local_log_stage(params, 'intrinsic ROI summary finished in %.1f sec', toc(t_stage));
end

stage.result = stage.save_result;
stage.result.runtime_sec = toc(t_run);
stage.result.intrinsic_activation_dir = stage.intrinsic_activation_dir;
stage.result.n_intrinsic_maps = stage.n_intrinsic_maps;
if isfield(stage.act_result, 'info_file')
    stage.result.intrinsic_activation_info_file = stage.act_result.info_file;
end
stage.result.intrinsic_roi_summary_dir = stage.intrinsic_roi_summary_dir;
stage.result.intrinsic_roi_summary_png = stage.intrinsic_roi_summary_png;
stage.result.n_intrinsic_roi_summary_modes = stage.n_intrinsic_roi_summary_modes;
if isfield(stage.roi_summary_result, 'info_file')
    stage.result.intrinsic_roi_summary_info_file = stage.roi_summary_result.info_file;
end
end


function cfg_ts = local_build_timescale_cfg(params, EDMD_outputs, ctx)
cfg_ts = struct();
if isfield(params, 'timescale') && isstruct(params.timescale)
    cfg_ts = params.timescale;
end

cfg_ts.max_modes_all = local_param_value(params, cfg_ts, ...
    'timescale_max_modes_all', 'max_modes_all', Inf);
cfg_ts.max_modes_sel = local_param_value(params, cfg_ts, ...
    'timescale_max_modes_sel', 'max_modes_sel', Inf);
cfg_ts.maxLag = local_param_value(params, cfg_ts, ...
    'timescale_max_lag', 'maxLag', []);
cfg_ts.xlim_time = local_param_value(params, cfg_ts, ...
    'timescale_xlim_time', 'xlim_time', 240);
cfg_ts.match_empirical_scale_to_theoretical = local_param_value(params, cfg_ts, ...
    'timescale_match_empirical_scale', ...
    'match_empirical_scale_to_theoretical', true);

cfg_ts.max_modes_all = local_cap_modes(cfg_ts.max_modes_all, ...
    local_count_all_modes(EDMD_outputs));
cfg_ts.max_modes_sel = local_cap_modes(cfg_ts.max_modes_sel, ...
    local_count_selected_modes(EDMD_outputs));
cfg_ts.dt = ctx.dt;
cfg_ts.t_plot = 1:size(EDMD_outputs.efuns, 1);
cfg_ts.title_prefix = sprintf('BOLD EDMD timescale diagnostics | %s | %s', ...
    ctx.run_info.dataset_id, ctx.run_info.run_name);
end


function value = local_param_value(params, cfg_ts, top_name, nested_name, default_value)
if isfield(params, top_name) && ~isempty(params.(top_name))
    value = params.(top_name);
elseif isfield(cfg_ts, nested_name)
    value = cfg_ts.(nested_name);
else
    value = default_value;
end
end


function n = local_count_all_modes(EDMD_outputs)
if isfield(EDMD_outputs, 'original_sorted') && ...
        isstruct(EDMD_outputs.original_sorted) && ...
        isfield(EDMD_outputs.original_sorted, 'efuns') && ...
        ~isempty(EDMD_outputs.original_sorted.efuns)
    n = size(EDMD_outputs.original_sorted.efuns, 2);
else
    n = local_count_selected_modes(EDMD_outputs);
end
end


function n = local_count_selected_modes(EDMD_outputs)
if isfield(EDMD_outputs, 'evalues') && ~isempty(EDMD_outputs.evalues)
    n = numel(EDMD_outputs.evalues);
else
    n = size(EDMD_outputs.efuns, 2);
end
end


function value = local_cap_modes(value, n_available)
if isempty(value)
    value = n_available;
elseif isfinite(value)
    value = min(double(value), n_available);
end
end


function local_log_stage(params, message, varargin)
if isfield(params, 'verbose') && ~params.verbose
    return;
end
fprintf(['[pipeline7] ', message, '\n'], varargin{:});
end
