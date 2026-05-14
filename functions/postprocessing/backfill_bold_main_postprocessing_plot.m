function [main_result, attempted] = backfill_bold_main_postprocessing_plot(post_file, params)
%BACKFILL_BOLD_MAIN_POSTPROCESSING_PLOT Regenerate the pipeline 7 main efun plot from an existing BOLD_POST.

attempted = false;
main_result = struct();

if nargin < 1 || isempty(post_file)
    return;
end
if nargin < 2 || ~isstruct(params) || ...
        ~isfield(params, 'make_main_plot') || ~params.make_main_plot
    return;
end
if exist(post_file, 'file') ~= 2
    return;
end

attempted = true;
try
    S = load(post_file, 'BOLD_POST');
    if ~isfield(S, 'BOLD_POST') || ~isstruct(S.BOLD_POST) || ...
            ~isfield(S.BOLD_POST, 'EDMD_outputs') || ~isstruct(S.BOLD_POST.EDMD_outputs)
        main_result.status = 'missing_payload';
        main_result.message = 'BOLD_POST.EDMD_outputs is missing.';
        return;
    end

    BOLD_POST = S.BOLD_POST;

    post_opts = local_get_field(params, 'post_opts', struct());
    post_opts.do_plot = true;
    post_opts.dt = local_resolve_dt(BOLD_POST, params);
    post_opts.time_vec = local_get_field(BOLD_POST, 'time_vec', []);
    post_opts.session_border = local_get_field(BOLD_POST, 'session_border_t', []);
    post_opts.draw_border = local_get_field(params, 'draw_session_borders', false);

    [EDMD_outputs, fig_main] = postprocess_EDMD_outputs(BOLD_POST.EDMD_outputs, post_opts);

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

    main_png = save_bold_reskoopnet_postprocessing_figure( ...
        fig_main, fig_dir, [run_name, '_efuns.png'], params);

    if ~isempty(fig_main) && isvalid(fig_main)
        close(fig_main);
    end

    BOLD_POST.EDMD_outputs = EDMD_outputs;
    if ~isfield(BOLD_POST, 'artifacts') || ~isstruct(BOLD_POST.artifacts)
        BOLD_POST.artifacts = struct();
    end
    BOLD_POST.artifacts.post_file = post_file;
    BOLD_POST.artifacts.main_png = main_png;
    save_bold_post_mat(post_file, BOLD_POST);

    main_result.status = 'ok';
    main_result.message = '';
    main_result.post_file = post_file;
    main_result.main_png = main_png;
    main_result.output_root = out_root;
catch ME
    main_result.status = 'error';
    main_result.message = ME.message;
    if exist('fig_main', 'var') && ~isempty(fig_main) && isvalid(fig_main)
        close(fig_main);
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
