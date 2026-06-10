function result = run_one_bold_reskoopnet_postprocessing_core(run_info, params)
%RUN_ONE_BOLD_RESKOOPNET_POSTPROCESSING_CORE Thin compatibility entry for one pipeline 7 run.
%
% Canonical interactive entry:
%   scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m
%
% The actual stage sequence now lives in:
%   prepare_bold_reskoopnet_postprocessing_run
%   run_bold_reskoopnet_postprocessing_stages
%
% This wrapper exists for batch runners and compatibility wrappers that
% still need a function-level entry point for one completed BOLD run.

if nargin < 1 || ~isstruct(run_info)
    error('run_info must be a struct.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

ctx = prepare_bold_reskoopnet_postprocessing_run(run_info, params);
if ctx.skip_existing
    result = local_make_result(ctx.run_info, 'skipped_existing', '', ...
        ctx.post_file, '', '', '', 0);
    result.observable_file = '';
    result.out_root = ctx.out_root;
    [main_result, attempted_main] = backfill_bold_main_postprocessing_plot(ctx.post_file, params);
    if attempted_main
        result.main_png = local_get_field(main_result, 'main_png', '');
        result.message = local_get_field(main_result, 'message', '');
        if strcmp(local_get_field(main_result, 'status', ''), 'ok')
            result.status = 'ok_main_backfill';
        end
    end
    [deconv_result, attempted_deconv] = backfill_bold_deconv_efun_plot(ctx.post_file, params);
    if attempted_deconv
        result.deconv_png = local_get_field(deconv_result, 'deconv_png', '');
        deconv_message = local_get_field(deconv_result, 'message', '');
        if isempty(result.message)
            result.message = deconv_message;
        elseif ~isempty(deconv_message)
            result.message = strjoin({result.message, deconv_message}, ' | ');
        end
        if strcmp(local_get_field(deconv_result, 'status', ''), 'ok')
            result.status = 'ok_deconv_backfill';
        end
    end
    [timescale_result, attempted_timescale] = backfill_bold_timescale_plot(ctx.post_file, params);
    if attempted_timescale
        result.timescale_png = local_get_field(timescale_result, 'timescale_png', '');
        timescale_message = local_get_field(timescale_result, 'message', '');
        if isempty(result.message)
            result.message = timescale_message;
        elseif ~isempty(timescale_message)
            result.message = strjoin({result.message, timescale_message}, ' | ');
        end
        if strcmp(local_get_field(timescale_result, 'status', ''), 'ok')
            result.status = 'ok_timescale_backfill';
        end
    end
    [act_result, attempted] = backfill_bold_intrinsic_activation_maps(ctx.post_file, params);
    if attempted
        result.intrinsic_activation_dir = local_get_field(act_result, 'activation_dir', '');
        result.n_intrinsic_maps = local_get_field(act_result, 'n_requested_maps', 0);
        result.intrinsic_activation_info_file = local_get_field(act_result, 'info_file', '');
        act_message = local_get_field(act_result, 'message', '');
        if isempty(result.message)
            result.message = act_message;
        elseif ~isempty(act_message)
            result.message = strjoin({result.message, act_message}, ' | ');
        end
        if strcmp(local_get_field(act_result, 'status', ''), 'ok')
            result.status = 'ok_intrinsic_only';
        end
    end
    [roi_result, roi_attempted] = backfill_bold_intrinsic_roi_bar_summaries(ctx.post_file, params);
    if roi_attempted
        result.intrinsic_roi_summary_dir = local_get_field(roi_result, 'summary_dir', '');
        result.intrinsic_roi_summary_png = local_get_field(roi_result, 'png_file', '');
        result.intrinsic_roi_summary_info_file = local_get_field(roi_result, 'info_file', '');
        result.n_intrinsic_roi_summary_modes = local_get_field(roi_result, 'n_selected_modes', 0);
        if strcmp(local_get_field(roi_result, 'status', ''), 'ok')
            result.status = 'ok_intrinsic_only';
        end
    end
    return;
end

stage = run_bold_reskoopnet_postprocessing_stages(ctx, params);
result = stage.result;
end


function params = local_apply_defaults(params)
defaults = build_bold_reskoopnet_postprocessing_params();
params = local_merge_defaults(defaults, params);
end


function out = local_merge_defaults(defaults, overrides)
out = defaults;
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    value = overrides.(name);
    if isstruct(value) && isfield(defaults, name) && isstruct(defaults.(name)) && ...
            isscalar(value) && isscalar(defaults.(name))
        out.(name) = local_merge_defaults(defaults.(name), value);
    elseif ~isempty(value)
        out.(name) = value;
    else
        out.(name) = defaults.(name);
    end
end
end


function result = local_make_result(run_info, status, message, post_file, ...
        main_png, deconv_png, timescale_png, runtime_sec)
result = struct();
result.run_info = run_info;
result.dataset_stem = run_info.dataset_stem;
result.run_name = run_info.run_name;
result.observable_mode = local_get_field(run_info, 'observable_mode', '');
result.residual_form = local_get_field(run_info, 'residual_form', '');
result.status = status;
result.message = message;
result.post_file = post_file;
result.main_png = main_png;
result.deconv_png = deconv_png;
result.timescale_png = timescale_png;
result.runtime_sec = runtime_sec;
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
