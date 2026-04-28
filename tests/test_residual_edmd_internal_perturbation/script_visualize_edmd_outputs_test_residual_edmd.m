this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

source_cfg = struct();
source_cfg.mode = 'chunk_dir';   % 'chunk_dir' | 'mat_file'
source_cfg.data_dir = 'E:\autodl_results\initial_point_test1_KV';
source_cfg.edmd_file = '';
source_cfg.concat = struct();
source_cfg.concat.filename_pattern = '*_outputs_*.mat';
source_cfg.concat.variable_name = 'EDMD_outputs';
source_cfg.concat.concat_fields = {'efuns'};
source_cfg.concat.concat_dim = 1;
source_cfg.concat.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form'};
source_cfg.concat.allow_missing_chunks = false;
source_cfg.concat.verbose = true;
source_cfg.concat.progress_every = 50;

viz_cfg = struct();
% viz_cfg.window_start = 1;
% viz_cfg.window_length = 3000;
viz_cfg.window_start = 120001;
viz_cfg.window_length = 3000;
viz_cfg.dt = [];
viz_cfg.max_basis = 30;
viz_cfg.do_timescale = true;
viz_cfg.do_deconv = true;
viz_cfg.do_deconv_similarity = true;
viz_cfg.draw_chunk_borders = false;
viz_cfg.deconv_methods = {'koopman_residual', 'wiener'};
viz_cfg.deconv_lambda_sources = {'edmd', 'empirical_complex'};
viz_cfg.deconv_plot_normalize_scopes = {'global', 'window'};
viz_cfg.deconv_normalize_exclude_idx = 1;
viz_cfg.deconv_similarity_fields = {'abs_all', 'real_all'};
viz_cfg.timescale_max_modes_all = 30;
viz_cfg.timescale_max_modes_sel = 20;
viz_cfg.timescale_maxLag = [];
viz_cfg.timescale_xlim_time = [];
viz_cfg.timescale_title_prefix = 'EDMD timescale diagnostics';
viz_cfg.save_figures = true;
viz_cfg.save_png = true;
viz_cfg.save_fig = true;
viz_cfg.save_dir = this_script_dir;
viz_cfg.save_tag = '';

[EDMD_outputs, concat_info, source_info] = local_load_edmd_source(source_cfg);

fprintf('EDMD source mode: %s\n', source_info.mode);
fprintf('EDMD source path:\n  %s\n', source_info.path);
fprintf('Loaded %d time samples and %d modes.\n', ...
    size(EDMD_outputs.efuns, 1), size(EDMD_outputs.efuns, 2));

fig_timescale = [];
timescale_info = [];
deconv_results = struct();
fig_deconv = struct();
fig_deconv_similarity = [];
deconv_similarity_info = [];

save_dir = '';
if viz_cfg.save_figures
    save_dir = local_prepare_save_dir(viz_cfg);
    fprintf('Saving figures to:\n  %s\n', save_dir);
end

opts = struct();
opts.abs_thresh = 0.01;
opts.sort_by = 'modulus';
opts.sort_dir = 'descend';
opts.max_basis = viz_cfg.max_basis;
opts.do_plot = true;
opts.window_start = viz_cfg.window_start;
opts.max_plot_samples = viz_cfg.window_length;
opts.draw_border = false;
opts.session_border = [];

if ~isempty(viz_cfg.dt)
    opts.dt = viz_cfg.dt;
    opts.time_vec = (0:size(EDMD_outputs.efuns, 1)-1) * viz_cfg.dt;
else
    opts.time_vec = [];
end

if viz_cfg.draw_chunk_borders && ~isempty(concat_info) && isfield(concat_info, 'chunk_end_idx')
    opts.session_border = local_convert_chunk_borders(concat_info.chunk_end_idx, opts.time_vec);
    opts.draw_border = true;
end

[EDMD_outputs_post, fig_main] = postprocess_EDMD_outputs(EDMD_outputs, opts);
local_save_figure_if_requested(fig_main, save_dir, '01_postprocess_main', viz_cfg);

if viz_cfg.do_timescale
    cfg_ts = struct();
    if ~isempty(viz_cfg.dt)
        cfg_ts.dt = viz_cfg.dt;
    end
    cfg_ts.t_plot = opts.window_start:(opts.window_start + opts.max_plot_samples - 1);
    cfg_ts.max_modes_all = min(viz_cfg.timescale_max_modes_all, ...
        size(EDMD_outputs_post.original_sorted.efuns, 2));
    cfg_ts.max_modes_sel = min(viz_cfg.timescale_max_modes_sel, ...
        numel(EDMD_outputs_post.evalues));
    if ~isempty(viz_cfg.timescale_maxLag)
        cfg_ts.maxLag = viz_cfg.timescale_maxLag;
    end
    if ~isempty(viz_cfg.timescale_xlim_time)
        cfg_ts.xlim_time = viz_cfg.timescale_xlim_time;
    elseif ~isempty(viz_cfg.dt)
        cfg_ts.xlim_time = viz_cfg.window_length * viz_cfg.dt;
    end
    cfg_ts.title_prefix = viz_cfg.timescale_title_prefix;

    [fig_timescale, timescale_info] = postprocess_EDMD_outputs_timescale(EDMD_outputs_post, cfg_ts);
    EDMD_outputs_post.timescale_info = timescale_info;
    local_save_figure_if_requested(fig_timescale, save_dir, '02_timescale_diagnostics', viz_cfg);
end

if viz_cfg.do_deconv
    for i_method = 1:numel(viz_cfg.deconv_methods)
        method_name = viz_cfg.deconv_methods{i_method};
        method_field = matlab.lang.makeValidName(lower(method_name));

        for i_lambda = 1:numel(viz_cfg.deconv_lambda_sources)
            lambda_source = viz_cfg.deconv_lambda_sources{i_lambda};
            lambda_field = matlab.lang.makeValidName(lower(lambda_source));

            for i_scope = 1:numel(viz_cfg.deconv_plot_normalize_scopes)
                scope_name = viz_cfg.deconv_plot_normalize_scopes{i_scope};
                scope_field = matlab.lang.makeValidName(lower(scope_name));

                cfg = struct();
                cfg.do_plot = true;
                cfg.window_start = viz_cfg.window_start;
                cfg.max_plot_samples = viz_cfg.window_length;
                cfg.max_modes_all = min(30, numel(EDMD_outputs_post.evalues));
                cfg.max_modes_sel = min(20, numel(EDMD_outputs_post.evalues));
                cfg.remove_mean = false;
                cfg.time_vec = opts.time_vec;
                cfg.draw_border = opts.draw_border;
                cfg.session_border = opts.session_border;
                cfg.method = method_name;
                cfg.lambda_source = lambda_source;
                cfg.plot_normalize_scope = scope_name;
                cfg.normalize_exclude_idx = viz_cfg.deconv_normalize_exclude_idx;

                fprintf('Plotting deconv comparison: method=%s, lambda=%s, normalization=%s\n', ...
                    method_name, lambda_source, scope_name);

                [~, fig_this, deconv_this] = ...
                    postprocess_EDMD_outputs_deconv_efuns(EDMD_outputs_post, cfg);

                if ~isempty(fig_this) && isgraphics(fig_this)
                    set(fig_this, 'Name', sprintf('EDMD deconv | %s | %s | %s', ...
                        method_name, lambda_source, scope_name), ...
                        'NumberTitle', 'off');
                end

                file_stub = sprintf('deconv__%s__%s__%s', ...
                    local_filename_safe(method_name), ...
                    local_filename_safe(lambda_source), ...
                    local_filename_safe(scope_name));
                local_save_figure_if_requested(fig_this, save_dir, file_stub, viz_cfg);

                deconv_results.(method_field).(lambda_field).(scope_field) = deconv_this;
                fig_deconv.(method_field).(lambda_field).(scope_field) = fig_this;
            end
        end
    end

    if viz_cfg.do_deconv_similarity
        sim_cfg = struct();
        sim_cfg.method_order = local_make_valid_name_list(viz_cfg.deconv_methods);
        sim_cfg.lambda_source_order = local_make_valid_name_list(viz_cfg.deconv_lambda_sources);
        sim_cfg.scope_order = local_make_valid_name_list(viz_cfg.deconv_plot_normalize_scopes);
        sim_cfg.order_mode = 'scope_first';
        sim_cfg.compare_fields = viz_cfg.deconv_similarity_fields;
        sim_cfg.title_prefix = 'Deconv condition similarity';

        [fig_deconv_similarity, deconv_similarity_info] = ...
            compare_deconv_condition_similarity(deconv_results, sim_cfg);
        local_save_figure_if_requested(fig_deconv_similarity, save_dir, ...
            '03_deconv_condition_similarity', viz_cfg);

        local_print_similarity_summary(deconv_similarity_info, viz_cfg.deconv_similarity_fields);
    end
end


function [EDMD_outputs, concat_info, source_info] = local_load_edmd_source(source_cfg)
% Load EDMD outputs either from a chunk directory or from a single MAT file.

switch lower(source_cfg.mode)
    case 'chunk_dir'
        fprintf('Loading EDMD outputs directly from chunk directory:\n  %s\n', source_cfg.data_dir);
        [EDMD_outputs, concat_info] = io_edmd.load_and_concat_edmd_output_chunks( ...
            source_cfg.data_dir, source_cfg.concat);
        source_info = struct();
        source_info.mode = 'chunk_dir';
        source_info.path = source_cfg.data_dir;

    case 'mat_file'
        edmd_file = source_cfg.edmd_file;
        if isempty(edmd_file)
            edmd_file = local_find_default_edmd_file(source_cfg.data_dir);
        end

        fprintf('Loading EDMD outputs from MAT file:\n  %s\n', edmd_file);
        S = load(edmd_file);

        if ~isfield(S, 'EDMD_outputs')
            error('File %s does not contain variable EDMD_outputs.', edmd_file);
        end

        EDMD_outputs = S.EDMD_outputs;
        if isfield(S, 'concat_info')
            concat_info = S.concat_info;
        else
            concat_info = [];
        end

        source_info = struct();
        source_info.mode = 'mat_file';
        source_info.path = edmd_file;

    otherwise
        error('Unknown source_cfg.mode = %s. Use ''chunk_dir'' or ''mat_file''.', source_cfg.mode);
end
end


function edmd_file = local_find_default_edmd_file(data_dir)
% Prefer a concatenated EDMD file when MAT-file mode is used.

patterns = {'*_concat.mat', '*_outputs_*_to_*_concat.mat'};

for i = 1:numel(patterns)
    L = dir(fullfile(data_dir, patterns{i}));
    if isempty(L)
        continue;
    end

    [~, order] = sort([L.datenum], 'descend');
    L = L(order);
    edmd_file = fullfile(L(1).folder, L(1).name);
    return;
end

error(['No concatenated EDMD MAT file was found in %s. ', ...
    'Set source_cfg.edmd_file explicitly or use source_cfg.mode = ''chunk_dir'' instead.'], data_dir);
end


function borders = local_convert_chunk_borders(chunk_end_idx, time_vec)
% Convert chunk borders to the same x-axis units used by the heatmaps.

chunk_end_idx = chunk_end_idx(:);

if isempty(time_vec)
    borders = chunk_end_idx;
    return;
end

keep = chunk_end_idx >= 1 & chunk_end_idx <= numel(time_vec);
borders = time_vec(chunk_end_idx(keep));
end


function local_print_similarity_summary(sim_info, fields)
% Print a concise textual summary of cross-condition similarity.
for i = 1:numel(fields)
    field_name = fields{i};
    if ~isfield(sim_info, field_name)
        continue;
    end
    S = sim_info.(field_name).summary;
    fprintf('\nSimilarity summary for %s:\n', field_name);
    fprintf('  mean corr   = %.4f\n', S.mean_corr);
    fprintf('  median corr = %.4f\n', S.median_corr);
    fprintf('  mean MAD    = %.4f\n', S.mean_mad);
    fprintf('  median MAD  = %.4f\n', S.median_mad);
    fprintf('  most similar by corr: %s  <->  %s   (corr=%.4f, MAD=%.4f)\n', ...
        S.most_similar_by_corr.label_i, S.most_similar_by_corr.label_j, ...
        S.most_similar_by_corr.corr, S.most_similar_by_corr.mad);
    fprintf('  least similar by corr: %s  <->  %s   (corr=%.4f, MAD=%.4f)\n', ...
        S.least_similar_by_corr.label_i, S.least_similar_by_corr.label_j, ...
        S.least_similar_by_corr.corr, S.least_similar_by_corr.mad);
end
end


function out = local_make_valid_name_list(names_in)
% Convert a cellstr list to valid MATLAB field names while preserving order.
out = cell(size(names_in));
for i = 1:numel(names_in)
    out{i} = matlab.lang.makeValidName(names_in{i});
end
end


function save_dir = local_prepare_save_dir(viz_cfg)
% Create a deterministic output directory for all figures from this run.
base_dir = viz_cfg.save_dir;
if isempty(base_dir)
    error('viz_cfg.save_dir must be non-empty when viz_cfg.save_figures=true.');
end

tag = viz_cfg.save_tag;
if isempty(tag)
    save_dir = base_dir;
else
    save_dir = fullfile(base_dir, local_filename_safe(tag));
end
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end
end


function local_save_figure_if_requested(fig_handle, save_dir, file_stub, viz_cfg)
% Save a figure to disk in the configured formats when requested.
if ~viz_cfg.save_figures || isempty(fig_handle) || ~isgraphics(fig_handle)
    return;
end

if isempty(save_dir)
    return;
end

file_stub = local_filename_safe(file_stub);
base_path = fullfile(save_dir, file_stub);

if viz_cfg.save_png
    try
        local_delete_if_exists([base_path '.png']);
        exportgraphics(fig_handle, [base_path '.png'], 'Resolution', 200);
    catch ME
        warning('Failed to save PNG for %s: %s', file_stub, ME.message);
    end
end

if viz_cfg.save_fig
    local_savefig_lightweight(fig_handle, [base_path '.fig'], file_stub);
end
end


function safe_name = local_filename_safe(name_in)
% Convert a label to a filesystem-friendly filename fragment.
safe_name = regexprep(char(name_in), '[^a-zA-Z0-9_\-]+', '_');
safe_name = regexprep(safe_name, '_+', '_');
safe_name = regexprep(safe_name, '^_+|_+$', '');
if isempty(safe_name)
    safe_name = 'untitled';
end
end


function local_savefig_lightweight(fig_handle, fig_path, file_stub)
% Save a .fig after temporarily stripping heavy UserData from the figure.
orig_user_data = [];
has_user_data = false;

try
    if isprop(fig_handle, 'UserData')
        orig_user_data = get(fig_handle, 'UserData');
        has_user_data = true;
        set(fig_handle, 'UserData', []);
    end

    local_delete_if_exists(fig_path);
    savefig(fig_handle, fig_path);
catch ME
    warning('Failed to save FIG for %s: %s', file_stub, ME.message);
end

if has_user_data && isgraphics(fig_handle)
    try
        set(fig_handle, 'UserData', orig_user_data);
    catch
        % Best effort restore; failure here should not stop the script.
    end
end
end


function local_delete_if_exists(file_path)
% Remove an existing output file before overwriting it.
if exist(file_path, 'file') == 2
    delete(file_path);
end
end
