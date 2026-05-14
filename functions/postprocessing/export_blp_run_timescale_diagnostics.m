function timescale_png = export_blp_run_timescale_diagnostics(run_info, params, EDMD_full)
%EXPORT_BLP_RUN_TIMESCALE_DIAGNOSTICS Export one condition-level timescale figure.

stage_dir = local_timescale_stage_dir(params, run_info.dataset);
file_stub = sprintf('%s_%s_timescale_diagnostics', ...
    run_info.dataset, local_filename_safe(run_info.run_name));
timescale_png = fullfile(stage_dir, [file_stub '.png']);

if ~params.do_timescale
    timescale_png = "";
    return;
end

if params.skip_existing && exist(timescale_png, 'file') == 2
    fprintf('  Reusing run-level timescale diagnostics.\n');
    timescale_png = string(timescale_png);
    return;
end

fprintf('  Run-level timescale diagnostics.\n');

if nargin < 3 || isempty(EDMD_full)
    EDMD_full = load_blp_run_edmd_outputs(run_info, params);
end
dt = local_resolve_edmd_dt(EDMD_full);

opts = struct();
opts.abs_thresh = params.abs_thresh;
opts.sort_by = params.sort_by;
opts.sort_dir = params.sort_dir;
opts.max_basis = params.max_basis;
if ~isempty(dt)
    opts.dt = dt;
end
EDMD_post = local_build_timescale_edmd_view(EDMD_full, opts, params);

cfg_ts = struct();
if ~isempty(dt)
    cfg_ts.dt = dt;
end
cfg_ts.t_plot = local_timescale_similarity_indices( ...
    local_edmd_view_num_samples(EDMD_post), params);
cfg_ts.max_modes_all = min(params.timescale_max_modes_all, ...
    numel(EDMD_post.original_sorted.evalues));
cfg_ts.max_modes_sel = min(params.timescale_max_modes_sel, numel(EDMD_post.evalues));
if isfield(params, 'timescale_max_lag')
    cfg_ts.maxLag = params.timescale_max_lag;
end
if isfield(params, 'timescale_xlim_time') && ~isempty(params.timescale_xlim_time)
    cfg_ts.xlim_time = params.timescale_xlim_time;
end
if isfield(params, 'timescale_match_empirical_scale') && ...
        ~isempty(params.timescale_match_empirical_scale)
    cfg_ts.match_empirical_scale_to_theoretical = ...
        params.timescale_match_empirical_scale;
end
cfg_ts.title_prefix = sprintf('EDMD timescale diagnostics | %s | %s', ...
    run_info.dataset, run_info.run_name);

[fig_ts, ~] = postprocess_EDMD_outputs_timescale(EDMD_post, cfg_ts);
if ~isempty(fig_ts) && isgraphics(fig_ts)
    set(fig_ts, 'Name', cfg_ts.title_prefix, 'NumberTitle', 'off');
end

timescale_png = string(local_save_fig(fig_ts, stage_dir, file_stub, params));
local_close_fig(fig_ts, params);
end


function stage_dir = local_timescale_stage_dir(params, dataset)
if isempty(params.output_root)
    stage_dir = io_project.get_pipeline_stage_dir( ...
        params.processed_root, dataset, 6, 'figures_timescale_diagnostics');
else
    stage_dir = fullfile( ...
        params.output_root, dataset, ...
        io_project.get_pipeline_stage_name(6, 'figures_timescale_diagnostics'));
end
end


function png_file = local_save_fig(fig_handle, stage_dir, file_stub, params)
if exist(stage_dir, 'dir') ~= 7
    mkdir(stage_dir);
end

png_file = fullfile(stage_dir, [file_stub '.png']);
fig_file = fullfile(stage_dir, [file_stub, '.fig']);

if isempty(fig_handle) || ~isgraphics(fig_handle)
    png_file = '';
    return;
end

drawnow;
if params.save_png
    exportgraphics(fig_handle, png_file, 'Resolution', params.resolution);
end
if params.save_fig
    savefig(fig_handle, fig_file);
end
end


function local_close_fig(fig_handle, params)
if params.close_figures && ~isempty(fig_handle) && isgraphics(fig_handle)
    close(fig_handle);
end
end


function dt = local_resolve_edmd_dt(EDMD_outputs)
dt = [];
dt_fields = {'dt', 'dx', 'sampling_period', 'sample_period'};
for i = 1:numel(dt_fields)
    field_name = dt_fields{i};
    if isfield(EDMD_outputs, field_name) && ~isempty(EDMD_outputs.(field_name))
        value = double(EDMD_outputs.(field_name));
        if isscalar(value) && isfinite(value) && value > 0
            dt = value;
            return;
        end
    end
end

fs_fields = {'fs', 'sampling_frequency'};
for i = 1:numel(fs_fields)
    field_name = fs_fields{i};
    if isfield(EDMD_outputs, field_name) && ~isempty(EDMD_outputs.(field_name))
        value = double(EDMD_outputs.(field_name));
        if isscalar(value) && isfinite(value) && value > 0
            dt = 1 / value;
            return;
        end
    end
end

warning('export_blp_run_timescale_diagnostics:MissingDt', ...
    'No dt/dx/fs found in EDMD outputs. Timescale diagnostics will use dt=1.');
end


function EDMD_view = local_build_timescale_edmd_view(EDMD_outputs, opts, params)
%LOCAL_BUILD_TIMESCALE_EDMD_VIEW Sort/select only modes needed by timescale plots.
%   Avoids postprocess_EDMD_outputs because that function normalizes the
%   whole full-basis eigenfunction matrix before timescale mode caps apply.

if ~isfield(EDMD_outputs, 'evalues') || ~isfield(EDMD_outputs, 'efuns')
    error('EDMD_outputs must contain evalues and efuns.');
end

evalues0 = EDMD_outputs.evalues(:);
efuns0 = EDMD_outputs.efuns;
K = numel(evalues0);
if size(efuns0, 2) ~= K
    error('Dimension mismatch: size(efuns,2)=%d but numel(evalues)=%d.', ...
        size(efuns0, 2), K);
end

switch lower(opts.sort_by)
    case {'modulus','abs'}
        key_full = abs(evalues0);
    case {'real','realpart'}
        key_full = real(evalues0);
    otherwise
        error('Unknown opts.sort_by = %s. Use ''modulus'' or ''real''.', opts.sort_by);
end
[~, ord_full] = sort(key_full, opts.sort_dir);
mask_sorted_full = abs(evalues0(ord_full)) > opts.abs_thresh;
if ~any(mask_sorted_full)
    error('Threshold removed all eigen-basis: abs_thresh=%g. Nothing left after masking.', ...
        opts.abs_thresh);
end
ord_thresholded_full = ord_full(mask_sorted_full);

if isfield(params, 'timescale_max_modes_all') && ...
        ~isempty(params.timescale_max_modes_all) && ...
        isfinite(params.timescale_max_modes_all)
    n_all = min(numel(ord_thresholded_full), max(0, floor(params.timescale_max_modes_all)));
else
    n_all = numel(ord_thresholded_full);
end
if n_all < 1
    error('timescale_max_modes_all must keep at least one mode.');
end

ord_all = ord_thresholded_full(1:n_all);
evalues_sorted_all = evalues0(ord_all);

ord_selected_full = ord_thresholded_full;
if isfinite(opts.max_basis)
    n_sel = min(numel(ord_selected_full), max(0, floor(opts.max_basis)));
else
    n_sel = numel(ord_selected_full);
end
if isfield(params, 'timescale_max_modes_sel') && ...
        ~isempty(params.timescale_max_modes_sel) && ...
        isfinite(params.timescale_max_modes_sel)
    n_sel = min(n_sel, max(0, floor(params.timescale_max_modes_sel)));
end
if n_sel < 1
    error('No selected modes remain for timescale diagnostics.');
end

ord_sel = ord_selected_full(1:n_sel);
evalues_sel = evalues0(ord_sel);

EDMD_view = struct();
EDMD_view.evalues = evalues_sel;
EDMD_view.source_efuns = efuns0;
EDMD_view.efun_col_indices = ord_sel(:);
if isfield(EDMD_outputs, 'kpm_modes') && size(EDMD_outputs.kpm_modes, 1) == K
    EDMD_view.kpm_modes = EDMD_outputs.kpm_modes(ord_sel, :);
end
EDMD_view.original_sorted = struct();
EDMD_view.original_sorted.evalues = evalues_sorted_all;
EDMD_view.original_sorted.source_efuns = efuns0;
EDMD_view.original_sorted.efun_col_indices = ord_all(:);
if isfield(EDMD_outputs, 'kpm_modes') && size(EDMD_outputs.kpm_modes, 1) == K
    EDMD_view.original_sorted.kpm_modes = EDMD_outputs.kpm_modes(ord_all, :);
end
EDMD_view.original_sorted.sort_by = opts.sort_by;
EDMD_view.original_sorted.sort_dir = opts.sort_dir;
EDMD_view.original_sorted.abs_thresh = opts.abs_thresh;
EDMD_view.original_sorted.max_basis = opts.max_basis;
EDMD_view.original_sorted.mode_cap = n_all;
EDMD_view.masked = struct();
EDMD_view.masked.evalues = evalues_sel;
EDMD_view.masked.efun_col_indices = ord_sel(:);
EDMD_view.masked.abs_thresh = opts.abs_thresh;
EDMD_view.masked.sort_by = opts.sort_by;
EDMD_view.masked.sort_dir = opts.sort_dir;
EDMD_view.masked.max_basis = opts.max_basis;
EDMD_view.idx_final_in_original = ord_sel(:);
end


function n_samples = local_edmd_view_num_samples(EDMD_view)
if isfield(EDMD_view, 'source_efuns') && ~isempty(EDMD_view.source_efuns)
    n_samples = size(EDMD_view.source_efuns, 1);
elseif isfield(EDMD_view, 'efuns') && ~isempty(EDMD_view.efuns)
    n_samples = size(EDMD_view.efuns, 1);
elseif isfield(EDMD_view, 'original_sorted') && ...
        isfield(EDMD_view.original_sorted, 'source_efuns') && ...
        ~isempty(EDMD_view.original_sorted.source_efuns)
    n_samples = size(EDMD_view.original_sorted.source_efuns, 1);
elseif isfield(EDMD_view, 'original_sorted') && ...
        isfield(EDMD_view.original_sorted, 'efuns') && ...
        ~isempty(EDMD_view.original_sorted.efuns)
    n_samples = size(EDMD_view.original_sorted.efuns, 1);
else
    error('Unable to resolve number of samples for timescale diagnostics.');
end
end


function idx = local_timescale_similarity_indices(n_samples, params)
max_samples = inf;
if isfield(params, 'timescale_similarity_max_samples') && ...
        ~isempty(params.timescale_similarity_max_samples)
    max_samples = double(params.timescale_similarity_max_samples);
end

if isfinite(max_samples) && max_samples > 0 && n_samples > max_samples
    step = max(1, ceil(n_samples / max_samples));
    idx = 1:step:n_samples;
else
    idx = 1:n_samples;
end
end


function out = local_filename_safe(name_in)
out = regexprep(char(string(name_in)), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled_run';
end
end
