function [timescale_result, attempted] = backfill_bold_timescale_plot(post_file, params)
%BACKFILL_BOLD_TIMESCALE_PLOT Regenerate the pipeline 7 timescale plot from an existing BOLD_POST.

attempted = false;
timescale_result = struct();

if nargin < 1 || isempty(post_file)
    return;
end
if nargin < 2 || ~isstruct(params) || ...
        ~isfield(params, 'make_timescale_plot') || ~params.make_timescale_plot
    return;
end
if exist(post_file, 'file') ~= 2
    return;
end

attempted = true;
fig_timescale = [];
try
    S = load(post_file, 'BOLD_POST');
    if ~isfield(S, 'BOLD_POST') || ~isstruct(S.BOLD_POST) || ...
            ~isfield(S.BOLD_POST, 'EDMD_outputs') || ~isstruct(S.BOLD_POST.EDMD_outputs)
        timescale_result.status = 'missing_payload';
        timescale_result.message = 'BOLD_POST.EDMD_outputs is missing.';
        return;
    end

    BOLD_POST = S.BOLD_POST;
    EDMD_outputs = BOLD_POST.EDMD_outputs;
    cfg_ts = local_build_timescale_cfg(params, BOLD_POST, EDMD_outputs);
    [fig_timescale, timescale_info] = postprocess_EDMD_outputs_timescale(EDMD_outputs, cfg_ts);

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

    timescale_png = save_bold_reskoopnet_postprocessing_figure( ...
        fig_timescale, fig_dir, [run_name, '_timescale.png'], params);

    BOLD_POST.timescale_info = timescale_info;
    if ~isfield(BOLD_POST, 'artifacts') || ~isstruct(BOLD_POST.artifacts)
        BOLD_POST.artifacts = struct();
    end
    BOLD_POST.artifacts.post_file = post_file;
    BOLD_POST.artifacts.timescale_png = timescale_png;
    save_bold_post_mat(post_file, BOLD_POST);

    timescale_result.status = 'ok';
    timescale_result.message = '';
    timescale_result.post_file = post_file;
    timescale_result.timescale_png = timescale_png;
    timescale_result.output_root = out_root;

    if ~isempty(fig_timescale) && isvalid(fig_timescale)
        close(fig_timescale);
    end
catch ME
    timescale_result.status = 'error';
    timescale_result.message = ME.message;
    if ~isempty(fig_timescale) && isvalid(fig_timescale)
        close(fig_timescale);
    end
end
end


function cfg_ts = local_build_timescale_cfg(params, BOLD_POST, EDMD_outputs)
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

dt = local_resolve_dt(BOLD_POST, params);
if ~isempty(dt)
    cfg_ts.dt = dt;
end
cfg_ts.t_plot = 1:size(EDMD_outputs.efuns, 1);
cfg_ts.max_modes_all = local_cap_modes(cfg_ts.max_modes_all, ...
    local_count_all_modes(EDMD_outputs));
cfg_ts.max_modes_sel = local_cap_modes(cfg_ts.max_modes_sel, ...
    local_count_selected_modes(EDMD_outputs));

run_info = local_get_field(BOLD_POST, 'run_info', struct());
dataset_id = local_get_field(run_info, 'dataset_id', '');
run_name = local_get_field(run_info, 'run_name', '');
cfg_ts.title_prefix = sprintf('BOLD EDMD timescale diagnostics | %s | %s', ...
    dataset_id, run_name);
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


function out_root = local_resolve_output_root(post_file)
mat_dir = fileparts(char(string(post_file)));
out_root = fileparts(mat_dir);
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
