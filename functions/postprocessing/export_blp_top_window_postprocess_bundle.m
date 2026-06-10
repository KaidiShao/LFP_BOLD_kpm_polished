function [run_table, run_manifest_file] = export_blp_top_window_postprocess_bundle(run_info, window_cache, params, EDMD_full)
%EXPORT_BLP_TOP_WINDOW_POSTPROCESS_BUNDLE Export the three window-level top-window figures.

if nargin < 4
    EDMD_full = [];
end

top_windows = window_cache.W.top_windows_table;
n_windows = min(params.n_top_windows, height(top_windows));
top_windows = top_windows(1:n_windows, :);
window_file = window_cache.window_file;

run_output_root = fullfile(local_dataset_output_root(params, run_info.dataset), ...
    run_info.run_name);
if exist(run_output_root, 'dir') ~= 7
    mkdir(run_output_root);
end

chunk_index = [];
if isempty(EDMD_full)
    chunk_index = local_build_chunk_index(run_info.output_dir, params.filename_pattern);
end
row_cells = cell(n_windows, 1);

for i_win = 1:n_windows
    idx1 = double(top_windows.global_start_idx(i_win));
    idx2 = double(top_windows.global_end_idx(i_win));
    file_stub = local_window_file_stub(top_windows, i_win);
    title_suffix = local_window_title(top_windows, i_win, idx1, idx2);

    fprintf('  [%02d/%02d] %s samples [%d, %d]\n', ...
        i_win, n_windows, file_stub, idx1, idx2);

    out = struct();
    out.dataset = string(run_info.dataset);
    out.observable_mode = string(run_info.observable_mode);
    out.residual_form = string(run_info.residual_form);
    out.run_name = string(run_info.run_name);
    out.run_output_dir = string(run_info.output_dir);
    out.state_diversity_file = string(window_file);
    out.window_rank = double(top_windows.state_diversity_rank(i_win));
    out.global_start_idx = idx1;
    out.global_end_idx = idx2;
    out.postprocess_main_png = "";
    out.deconv_png = "";
    out.deconv_localwin_png = "";
    out.status = "ok";
    out.error_message = "";

    try
        out.postprocess_main_png = string(local_stage_png_path( ...
            run_output_root, '01_postprocess_main', file_stub));
        out.deconv_png = string(local_stage_png_path( ...
            run_output_root, '03_deconv', file_stub));
        out.deconv_localwin_png = string(local_stage_png_path( ...
            run_output_root, '04_deconv_localwin_norm', file_stub));

        if params.skip_existing && local_requested_window_outputs_exist(out, params)
            fprintf('    Skipping existing window outputs.\n');
            row_cells{i_win} = struct2table(out, 'AsArray', true);
            continue;
        end

        EDMD_window = local_load_edmd_window(run_info, chunk_index, EDMD_full, idx1, idx2, params.variable_name);

        [EDMD_post, fig_main] = local_run_main_postprocess( ...
            EDMD_window, idx1, idx2, title_suffix, params);
        if params.do_main_plot
            out.postprocess_main_png = string(local_save_fig( ...
                fig_main, run_output_root, '01_postprocess_main', file_stub, params));
        end
        local_close_fig(fig_main, params);

        if params.do_deconv
            [~, fig_deconv, ~] = local_run_deconv_postprocess( ...
                EDMD_post, idx1, idx2, title_suffix, params, params.deconv_plot_normalize_scope);
            out.deconv_png = string(local_save_fig( ...
                fig_deconv, run_output_root, '03_deconv', file_stub, params));
            local_close_fig(fig_deconv, params);
        end

        if params.do_deconv_window_norm
            [~, fig_deconv_local, ~] = local_run_deconv_postprocess( ...
                EDMD_post, idx1, idx2, title_suffix, params, params.deconv_window_plot_normalize_scope);
            out.deconv_localwin_png = string(local_save_fig( ...
                fig_deconv_local, run_output_root, '04_deconv_localwin_norm', file_stub, params));
            local_close_fig(fig_deconv_local, params);
        end
    catch ME
        out.status = "failed";
        out.error_message = string(local_single_line_error(ME));
        warning('Pipeline window failed: %s', out.error_message);
        close all force;
    end

    row_cells{i_win} = struct2table(out, 'AsArray', true);
end

run_table = vertcat(row_cells{:});
run_manifest_file = fullfile(run_output_root, ...
    sprintf('%s__top%d_state_diversity_postprocessing_manifest.csv', ...
    run_info.run_name, n_windows));
writetable(run_table, run_manifest_file);
fprintf('  Run manifest:\n    %s\n', run_manifest_file);
end


function dataset_output_root = local_dataset_output_root(params, dataset)
if isempty(params.output_root)
    dataset_output_root = io_project.get_pipeline_stage_dir( ...
        params.processed_root, dataset, 6, 'top_state_diversity_postprocessing');
else
    dataset_output_root = fullfile(params.output_root, dataset);
end
end


function chunk_index = local_build_chunk_index(output_dir, filename_pattern)
files = local_collect_chunk_files(output_dir, filename_pattern);
if isempty(files)
    error('No EDMD output chunks found in %s.', output_dir);
end

S = load(files(1).fullpath, 'EDMD_outputs');
if ~isfield(S, 'EDMD_outputs') || ~isfield(S.EDMD_outputs, 'efuns')
    error('First chunk does not contain EDMD_outputs.efuns: %s', files(1).fullpath);
end

base_chunk_len = size(S.EDMD_outputs.efuns, 1);
if base_chunk_len <= 0
    error('Invalid base chunk length in %s.', files(1).fullpath);
end

chunk_index = struct();
chunk_index.files = files;
chunk_index.base_chunk_len = base_chunk_len;
chunk_index.output_dir = output_dir;
end


function files = local_collect_chunk_files(output_dir, filename_pattern)
L = dir(fullfile(output_dir, filename_pattern));
files = repmat(struct('name', '', 'folder', '', 'fullpath', '', ...
    'chunk_id', NaN, 'datenum', NaN), numel(L), 1);
n_keep = 0;
for i = 1:numel(L)
    name = L(i).name;
    if contains(name, '_outputs_Psi_')
        continue;
    end
    tokens = regexp(name, '^(.*)_outputs_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tokens)
        continue;
    end
    n_keep = n_keep + 1;
    files(n_keep).name = name;
    files(n_keep).folder = L(i).folder;
    files(n_keep).fullpath = fullfile(L(i).folder, name);
    files(n_keep).chunk_id = str2double(tokens{2});
    files(n_keep).datenum = L(i).datenum;
end
files = files(1:n_keep);
if isempty(files)
    return;
end
[~, order] = sort([files.chunk_id]);
files = files(order);
end


function EDMD_window = local_load_edmd_window(run_info, chunk_index, EDMD_full, idx1, idx2, variable_name)
if ~isempty(EDMD_full)
    EDMD_window = local_slice_edmd_window_from_full(run_info, EDMD_full, idx1, idx2);
    return;
end

if idx1 < 1 || idx2 < idx1
    error('Invalid requested window [%d, %d].', idx1, idx2);
end

base_len = chunk_index.base_chunk_len;
start_chunk = floor((idx1 - 1) / base_len) + 1;
end_chunk = floor((idx2 - 1) / base_len) + 1;

if start_chunk < 1 || end_chunk > numel(chunk_index.files)
    error('Window [%d, %d] requires chunk %d..%d, but only %d chunks exist.', ...
        idx1, idx2, start_chunk, end_chunk, numel(chunk_index.files));
end

pieces = cell(end_chunk - start_chunk + 1, 1);
ref = [];
for chunk_id = start_chunk:end_chunk
    S = load(chunk_index.files(chunk_id).fullpath, variable_name);
    if ~isfield(S, variable_name)
        error('Chunk does not contain variable %s: %s', ...
            variable_name, chunk_index.files(chunk_id).fullpath);
    end
    current = S.(variable_name);
    if ~isfield(current, 'efuns')
        error('Chunk does not contain %s.efuns: %s', ...
            variable_name, chunk_index.files(chunk_id).fullpath);
    end
    if isempty(ref)
        ref = current;
    end

    chunk_start = (chunk_id - 1) * base_len + 1;
    chunk_end = chunk_start + size(current.efuns, 1) - 1;
    local_start = max(idx1, chunk_start) - chunk_start + 1;
    local_end = min(idx2, chunk_end) - chunk_start + 1;
    if local_start > local_end
        continue;
    end
    pieces{chunk_id - start_chunk + 1} = current.efuns(local_start:local_end, :);
end

pieces = pieces(~cellfun(@isempty, pieces));
if isempty(pieces)
    error('No samples were loaded for requested window [%d, %d].', idx1, idx2);
end

efuns = vertcat(pieces{:});
expected_len = idx2 - idx1 + 1;
if size(efuns, 1) ~= expected_len
    error('Loaded %d samples for window [%d, %d], expected %d.', ...
        size(efuns, 1), idx1, idx2, expected_len);
end

EDMD_window = ref;
EDMD_window.efuns = efuns;
EDMD_window.global_start_idx = idx1;
EDMD_window.global_end_idx = idx2;
end


function EDMD_window = local_slice_edmd_window_from_full(run_info, EDMD_full, idx1, idx2)
if idx1 < 1 || idx2 < idx1
    error('Invalid requested window [%d, %d].', idx1, idx2);
end
if ~isfield(EDMD_full, 'efuns') || isempty(EDMD_full.efuns)
    error('EDMD_full.efuns must be available for in-workspace window slicing.');
end
if idx2 > size(EDMD_full.efuns, 1)
    error('Requested window [%d, %d] exceeds loaded EDMD length %d for run %s.', ...
        idx1, idx2, size(EDMD_full.efuns, 1), run_info.run_name);
end

EDMD_window = EDMD_full;
EDMD_window.efuns = EDMD_full.efuns(idx1:idx2, :);
EDMD_window.global_start_idx = idx1;
EDMD_window.global_end_idx = idx2;
end


function [EDMD_post, fig_main] = local_run_main_postprocess(EDMD_window, idx1, idx2, title_suffix, params)
T = size(EDMD_window.efuns, 1);
opts = struct();
opts.abs_thresh = params.abs_thresh;
opts.sort_by = params.sort_by;
opts.sort_dir = params.sort_dir;
opts.max_basis = params.max_basis;
opts.do_plot = params.do_main_plot;
opts.window_start = 1;
opts.max_plot_samples = T;
opts.time_vec = (idx1:idx2).';
opts.draw_border = false;

[EDMD_post, fig_main] = postprocess_EDMD_outputs(EDMD_window, opts);
if ~isempty(fig_main) && isgraphics(fig_main)
    set(fig_main, 'Name', ['EDMD postprocess | ', title_suffix], ...
        'NumberTitle', 'off');
    sgtitle(['EDMD postprocess | ', title_suffix], 'Interpreter', 'none');
end
end


function [EDMD_with_deconv, fig_deconv, deconv] = local_run_deconv_postprocess( ...
    EDMD_post, idx1, idx2, title_suffix, params, normalize_scope)
if nargin < 6 || isempty(normalize_scope)
    normalize_scope = params.deconv_plot_normalize_scope;
end
T = size(EDMD_post.efuns, 1);
cfg_deconv = struct();
cfg_deconv.do_plot = true;
cfg_deconv.window_start = 1;
cfg_deconv.max_plot_samples = T;
cfg_deconv.max_modes_all = min(params.deconv_max_modes_all, local_n_thresholded_modes(EDMD_post));
cfg_deconv.max_modes_sel = min(params.deconv_max_modes_sel, numel(EDMD_post.evalues));
cfg_deconv.remove_mean = false;
cfg_deconv.time_vec = (idx1:idx2).';
cfg_deconv.draw_border = false;
cfg_deconv.method = params.deconv_method;
cfg_deconv.lambda_source = params.deconv_lambda_source;
cfg_deconv.plot_normalize_scope = normalize_scope;
cfg_deconv.normalize_exclude_idx = params.deconv_normalize_exclude_idx;

[EDMD_with_deconv, fig_deconv, deconv] = ...
    postprocess_EDMD_outputs_deconv_efuns(EDMD_post, cfg_deconv);
if ~isempty(fig_deconv) && isgraphics(fig_deconv)
    set(fig_deconv, 'Name', ['EDMD deconv | ', title_suffix], 'NumberTitle', 'off');
end
end


function n_modes = local_n_thresholded_modes(EDMD_post)
if isfield(EDMD_post, 'thresholded_sorted') && ...
        isfield(EDMD_post.thresholded_sorted, 'evalues') && ...
        ~isempty(EDMD_post.thresholded_sorted.evalues)
    n_modes = numel(EDMD_post.thresholded_sorted.evalues);
elseif isfield(EDMD_post, 'original_sorted') && ...
        isfield(EDMD_post.original_sorted, 'evalues') && ...
        ~isempty(EDMD_post.original_sorted.evalues)
    if isfield(EDMD_post.original_sorted, 'abs_thresh') && ...
            ~isempty(EDMD_post.original_sorted.abs_thresh)
        n_modes = nnz(abs(EDMD_post.original_sorted.evalues) > ...
            EDMD_post.original_sorted.abs_thresh);
    else
        n_modes = numel(EDMD_post.original_sorted.evalues);
    end
else
    n_modes = numel(EDMD_post.evalues);
end
end


function png_file = local_save_fig(fig_handle, run_output_root, stage_name, file_stub, params)
stage_dir = fullfile(run_output_root, stage_name);
if exist(stage_dir, 'dir') ~= 7
    mkdir(stage_dir);
end

png_file = local_stage_png_path(run_output_root, stage_name, file_stub);
fig_file = fullfile(stage_dir, [file_stub, '.fig']);

if params.skip_existing && exist(png_file, 'file') == 2
    return;
end

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


function png_file = local_stage_png_path(run_output_root, stage_name, file_stub)
png_file = fullfile(run_output_root, stage_name, [file_stub, '.png']);
end


function tf = local_requested_window_outputs_exist(out, params)
if ~params.save_png
    tf = false;
    return;
end

required_files = strings(0, 1);
if params.do_main_plot
    required_files(end+1, 1) = out.postprocess_main_png;
end
if params.do_deconv
    required_files(end+1, 1) = out.deconv_png;
end
if params.do_deconv_window_norm
    required_files(end+1, 1) = out.deconv_localwin_png;
end

if isempty(required_files)
    tf = false;
    return;
end

tf = true;
for i = 1:numel(required_files)
    if strlength(required_files(i)) == 0 || exist(char(required_files(i)), 'file') ~= 2
        tf = false;
        return;
    end
end
end


function local_close_fig(fig_handle, params)
if params.close_figures && ~isempty(fig_handle) && isgraphics(fig_handle)
    close(fig_handle);
end
end


function file_stub = local_window_file_stub(top_windows, i_win)
rank_value = double(top_windows.state_diversity_rank(i_win));
if ismember('global_window_idx', top_windows.Properties.VariableNames)
    file_stub = sprintf('rank_%02d_globalwin_%03d', ...
        rank_value, double(top_windows.global_window_idx(i_win)));
else
    file_stub = sprintf('rank_%02d_row_%03d', rank_value, i_win);
end
end


function title_text = local_window_title(top_windows, i_win, idx1, idx2)
rank_value = double(top_windows.state_diversity_rank(i_win));
if ismember('global_window_idx', top_windows.Properties.VariableNames)
    title_text = sprintf('rank %d | globalwin %d | samples [%d,%d]', ...
        rank_value, double(top_windows.global_window_idx(i_win)), idx1, idx2);
else
    title_text = sprintf('rank %d | samples [%d,%d]', rank_value, idx1, idx2);
end

if ismember('active_state_richness', top_windows.Properties.VariableNames)
    title_text = sprintf('%s | richness=%d', title_text, ...
        double(top_windows.active_state_richness(i_win)));
end
if ismember('normalized_state_entropy', top_windows.Properties.VariableNames)
    title_text = sprintf('%s | Hnorm=%.3f', title_text, ...
        double(top_windows.normalized_state_entropy(i_win)));
end
end


function msg = local_single_line_error(ME)
msg = getReport(ME, 'basic', 'hyperlinks', 'off');
msg = regexprep(msg, '\s+', ' ');
msg = strtrim(msg);
end
