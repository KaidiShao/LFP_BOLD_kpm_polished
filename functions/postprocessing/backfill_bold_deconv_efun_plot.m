function [deconv_result, attempted] = backfill_bold_deconv_efun_plot(post_file, params)
%BACKFILL_BOLD_DECONV_EFUN_PLOT Regenerate the pipeline 7 deconv efun plot from an existing BOLD_POST.

attempted = false;
deconv_result = struct();

if nargin < 1 || isempty(post_file)
    return;
end
if nargin < 2 || ~isstruct(params) || ...
        ~isfield(params, 'make_deconv_plot') || ~params.make_deconv_plot
    return;
end
if exist(post_file, 'file') ~= 2
    return;
end

attempted = true;
fig_deconv = [];
try
    S = load(post_file, 'BOLD_POST');
    if ~isfield(S, 'BOLD_POST') || ~isstruct(S.BOLD_POST) || ...
            ~isfield(S.BOLD_POST, 'EDMD_outputs') || ~isstruct(S.BOLD_POST.EDMD_outputs)
        deconv_result.status = 'missing_payload';
        deconv_result.message = 'BOLD_POST.EDMD_outputs is missing.';
        return;
    end

    BOLD_POST = S.BOLD_POST;
    EDMD_outputs = BOLD_POST.EDMD_outputs;

    deconv_cfg = local_get_field(params, 'deconv', struct());
    deconv_cfg.do_plot = true;
    deconv_cfg.dt = local_resolve_dt(BOLD_POST, params);
    deconv_cfg.time_vec = local_get_field(BOLD_POST, 'time_vec', []);
    deconv_cfg.session_border = local_get_field(BOLD_POST, 'session_border_t', []);
    deconv_cfg.draw_border = local_get_field(params, 'draw_session_borders', false);
    deconv_cfg.save_path = '';

    [EDMD_outputs, fig_deconv, deconv] = ...
        postprocess_EDMD_outputs_deconv_efuns(EDMD_outputs, deconv_cfg);

    out_root = local_resolve_output_root(post_file);
    fig_dir = fullfile(out_root, 'fig');
    if exist(fig_dir, 'dir') ~= 7
        mkdir(fig_dir);
    end

    run_info = local_get_field(BOLD_POST, 'run_info', struct());
    run_name = local_get_field(run_info, 'run_name', '');
    if isempty(run_name)
        [~, run_name] = fileparts(post_file);
        run_name = regexprep(run_name, '_bold_post$', '');
    end

    deconv_png = save_bold_reskoopnet_postprocessing_figure( ...
        fig_deconv, fig_dir, [run_name, '_deconv_efuns.png'], params);

    if ~isempty(fig_deconv) && isvalid(fig_deconv)
        close(fig_deconv);
    end

    BOLD_POST.EDMD_outputs = EDMD_outputs;
    BOLD_POST.deconv = deconv;
    if ~isfield(BOLD_POST, 'artifacts') || ~isstruct(BOLD_POST.artifacts)
        BOLD_POST.artifacts = struct();
    end
    BOLD_POST.artifacts.post_file = post_file;
    BOLD_POST.artifacts.deconv_png = deconv_png;
    save_bold_post_mat(post_file, BOLD_POST);

    deconv_result.status = 'ok';
    deconv_result.message = '';
    deconv_result.post_file = post_file;
    deconv_result.deconv_png = deconv_png;
    deconv_result.output_root = out_root;
catch ME
    deconv_result.status = 'error';
    deconv_result.message = ME.message;
    if ~isempty(fig_deconv) && isvalid(fig_deconv)
        close(fig_deconv);
    end
end
end


function out_root = local_resolve_output_root(post_file)
mat_dir = fileparts(char(string(post_file)));
out_root = fileparts(mat_dir);
end


function dt = local_resolve_dt(BOLD_POST, params)
dt = local_get_field(BOLD_POST, 'dt', []);
if isempty(dt) && isfield(BOLD_POST, 'EDMD_outputs') && isstruct(BOLD_POST.EDMD_outputs)
    edmd = BOLD_POST.EDMD_outputs;
    fields = {'dt', 'dx', 'sampling_period', 'sample_period'};
    for i = 1:numel(fields)
        name = fields{i};
        if isfield(edmd, name) && ~isempty(edmd.(name))
            dt = double(edmd.(name));
            break;
        end
    end
    if isempty(dt) && isfield(edmd, 'fs') && ~isempty(edmd.fs)
        dt = 1 / double(edmd.fs);
    end
end
if isempty(dt)
    dt = local_get_field(params, 'default_dt', []);
end
if isempty(dt) || ~isfinite(dt(1)) || dt(1) <= 0
    dt = [];
else
    dt = dt(1);
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
