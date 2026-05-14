function [BOLD_POST, result] = save_bold_reskoopnet_postprocessing_output( ...
        ctx, EDMD_outputs, deconv, timescale_info, params, artifacts, runtime_sec)
%SAVE_BOLD_RESKOOPNET_POSTPROCESSING_OUTPUT Save the canonical pipeline 7 BOLD_POST artifact.

if nargin < 7 || isempty(runtime_sec)
    runtime_sec = NaN;
end

main_png = local_get_artifact_path(artifacts, 'main_png');
deconv_png = local_get_artifact_path(artifacts, 'deconv_png');
timescale_png = local_get_artifact_path(artifacts, 'timescale_png');

BOLD_POST = struct();
BOLD_POST.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
BOLD_POST.run_info = ctx.run_info;
BOLD_POST.source_info = ctx.source_info;
BOLD_POST.concat_info = ctx.concat_info;
BOLD_POST.observable_file = ctx.observable_file;
BOLD_POST.observable = ctx.obs_meta;
BOLD_POST.session = ctx.session;
BOLD_POST.dt = ctx.dt;
BOLD_POST.time_vec = ctx.time_vec;
BOLD_POST.session_border_t = ctx.session_border_t;
BOLD_POST.EDMD_outputs = EDMD_outputs;
BOLD_POST.deconv = deconv;
BOLD_POST.timescale_info = timescale_info;
BOLD_POST.params = params;
BOLD_POST.source_file = ctx.post_file;
BOLD_POST.artifacts = struct('post_file', ctx.post_file, ...
    'main_png', main_png, 'deconv_png', deconv_png, ...
    'timescale_png', timescale_png);

save_bold_post_mat(ctx.post_file, BOLD_POST);

result = struct();
result.run_info = ctx.run_info;
result.dataset_stem = ctx.run_info.dataset_stem;
result.run_name = ctx.run_info.run_name;
result.observable_mode = local_get_field(ctx.run_info, 'observable_mode', '');
result.residual_form = local_get_field(ctx.run_info, 'residual_form', '');
result.status = 'ok';
result.message = '';
result.post_file = ctx.post_file;
result.main_png = main_png;
result.deconv_png = deconv_png;
result.timescale_png = timescale_png;
result.runtime_sec = runtime_sec;
result.observable_file = ctx.observable_file;
result.out_root = ctx.out_root;
result.timescale_info = timescale_info;
end


function path_out = local_get_artifact_path(artifacts, name)
if isstruct(artifacts) && isfield(artifacts, name) && ~isempty(artifacts.(name))
    path_out = artifacts.(name);
else
    path_out = '';
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
