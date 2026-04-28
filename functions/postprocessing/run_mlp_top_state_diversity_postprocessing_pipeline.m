function manifest = run_mlp_top_state_diversity_postprocessing_pipeline(params)
%RUN_MLP_TOP_STATE_DIVERSITY_POSTPROCESSING_PIPELINE
% Run pre-reduction EDMD/MLP postprocessing on consensus-state-diversity top windows.
%
% The pipeline intentionally uses only windows from the saved
% consensus-state-diversity top table. It does not pick arbitrary windows.

if nargin < 1
    params = struct();
end

params = local_apply_defaults(params);

if params.headless
    set(groot, 'defaultFigureVisible', 'off');
end
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
end

if exist(params.combined_manifest_root, 'dir') ~= 7
    mkdir(params.combined_manifest_root);
end

runs = local_discover_completed_runs(params);
if isempty(runs)
    error('No completed MLP output runs were found under %s.', params.autodl_root);
end

fprintf('Discovered %d completed unique MLP conditions.\n', numel(runs));
for i = 1:numel(runs)
    fprintf('  %s | %s | %s | %s\n', ...
        runs(i).dataset, runs(i).observable_mode, runs(i).residual_form, runs(i).run_name);
end

manifest_rows = table();
spike_manifest_rows = table();
dataset_cache = struct();

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n[%d/%d] Processing %s | %s | %s\n', ...
        i_run, numel(runs), run_info.dataset, ...
        run_info.observable_mode, run_info.residual_form);
    fprintf('Source EDMD chunks:\n  %s\n', run_info.output_dir);

    cfg = local_dataset_cfg(run_info.dataset);
    dataset_key = matlab.lang.makeValidName(run_info.dataset);
    if ~isfield(dataset_cache, dataset_key)
        dataset_cache.(dataset_key) = local_load_top_windows_for_dataset(cfg, params);
    end
    W = dataset_cache.(dataset_key).W;
    window_file = dataset_cache.(dataset_key).window_file;

    top_windows = W.top_windows_table;
    n_windows = min(params.n_top_windows, height(top_windows));
    top_windows = top_windows(1:n_windows, :);

    run_output_root = fullfile(local_dataset_output_root(params, run_info.dataset), ...
        run_info.run_name);
    if exist(run_output_root, 'dir') ~= 7
        mkdir(run_output_root);
    end

    chunk_index = local_build_chunk_index(run_info.output_dir, params.filename_pattern);

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
        out.timescale_png = "";
        out.deconv_png = "";
        out.deconv_localwin_png = "";
        out.status = "ok";
        out.error_message = "";

        try
            out.postprocess_main_png = string(local_stage_png_path( ...
                run_output_root, '01_postprocess_main', file_stub));
            out.timescale_png = string(local_stage_png_path( ...
                run_output_root, '02_timescale', file_stub));
            out.deconv_png = string(local_stage_png_path( ...
                run_output_root, '03_deconv', file_stub));
            out.deconv_localwin_png = string(local_stage_png_path( ...
                run_output_root, '04_deconv_localwin_norm', file_stub));

            if params.skip_existing && local_requested_window_outputs_exist(out, params)
                fprintf('    Skipping existing window outputs.\n');
                row_cells{i_win} = struct2table(out, 'AsArray', true);
                continue;
            end

            EDMD_window = local_load_edmd_window(chunk_index, idx1, idx2, params.variable_name);

            [EDMD_post, fig_main] = local_run_main_postprocess( ...
                EDMD_window, idx1, idx2, title_suffix, params);
            if params.do_main_plot
                out.postprocess_main_png = string(local_save_fig( ...
                    fig_main, run_output_root, '01_postprocess_main', ...
                    file_stub, params));
            end
            local_close_fig(fig_main, params);

            if params.do_timescale
                [fig_ts, ~] = local_run_timescale_postprocess( ...
                    EDMD_post, title_suffix, params);
                out.timescale_png = string(local_save_fig( ...
                    fig_ts, run_output_root, '02_timescale', ...
                    file_stub, params));
                local_close_fig(fig_ts, params);
            end

            if params.do_deconv
                [~, fig_deconv, ~] = local_run_deconv_postprocess( ...
                    EDMD_post, idx1, idx2, title_suffix, ...
                    params, params.deconv_plot_normalize_scope);
                out.deconv_png = string(local_save_fig( ...
                    fig_deconv, run_output_root, '03_deconv', ...
                    file_stub, params));
                local_close_fig(fig_deconv, params);
            end

            if params.do_deconv_window_norm
                [~, fig_deconv_local, ~] = local_run_deconv_postprocess( ...
                    EDMD_post, idx1, idx2, title_suffix, ...
                    params, params.deconv_window_plot_normalize_scope);
                out.deconv_localwin_png = string(local_save_fig( ...
                    fig_deconv_local, run_output_root, '04_deconv_localwin_norm', ...
                    file_stub, params));
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

    manifest_rows = [manifest_rows; run_table]; %#ok<AGROW>

    if params.do_spike_correlation
        spike_row = local_run_spike_residual_correlation_fulltime( ...
            cfg, run_info, run_output_root, params);
        spike_manifest_rows = [spike_manifest_rows; spike_row]; %#ok<AGROW>
    end
end

manifest = struct();
manifest.params = params;
manifest.runs = runs;
manifest.table = manifest_rows;
manifest.spike_table = spike_manifest_rows;
manifest.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

manifest.csv_file = fullfile(params.combined_manifest_root, ...
    sprintf('mlp_top_state_diversity_postprocessing_manifest_%s.csv', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
writetable(manifest_rows, manifest.csv_file);

manifest.mat_file = fullfile(params.combined_manifest_root, ...
    sprintf('mlp_top_state_diversity_postprocessing_manifest_%s.mat', ...
    char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));

if ~isempty(spike_manifest_rows)
    manifest.spike_csv_file = fullfile(params.combined_manifest_root, ...
        sprintf('mlp_fulltime_spike_residual_correlation_manifest_%s.csv', ...
        char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'))));
    writetable(spike_manifest_rows, manifest.spike_csv_file);
    fprintf('Full-time spike correlation manifest CSV:\n  %s\n', ...
        manifest.spike_csv_file);
else
    manifest.spike_csv_file = '';
end

save(manifest.mat_file, 'manifest', '-v7.3');

fprintf('\nCombined manifest CSV:\n  %s\n', manifest.csv_file);
fprintf('Combined manifest MAT:\n  %s\n', manifest.mat_file);
end


function params = local_apply_defaults(params)
if ~isfield(params, 'dataset_stems') || isempty(params.dataset_stems)
    params.dataset_stems = {'e10fV1', 'e10gb1', 'e10gh1'};
end
params.dataset_stems = cellstr(string(params.dataset_stems(:)).');

if ~isfield(params, 'exclude_dataset_stems')
    params.exclude_dataset_stems = {'f12m01'};
elseif isempty(params.exclude_dataset_stems)
    params.exclude_dataset_stems = {};
end
params.exclude_dataset_stems = lower(cellstr(string(params.exclude_dataset_stems(:)).'));

if ~isfield(params, 'autodl_root') || isempty(params.autodl_root)
    params.autodl_root = 'E:\autodl_results';
end

if ~isfield(params, 'processed_root') || isempty(params.processed_root)
    params.processed_root = io_project.get_project_processed_root();
end

if ~isfield(params, 'output_folder_name') || isempty(params.output_folder_name)
    params.output_folder_name = 'mlp_top_state_diversity_postprocessing';
end

if ~isfield(params, 'output_root') || isempty(params.output_root)
    params.output_root = '';
end

if ~isfield(params, 'combined_manifest_root') || isempty(params.combined_manifest_root)
    if isempty(params.output_root)
        params.combined_manifest_root = fullfile(params.processed_root, ...
            'postprocessing_manifests', params.output_folder_name);
    else
        params.combined_manifest_root = params.output_root;
    end
end

if ~isfield(params, 'filename_pattern') || isempty(params.filename_pattern)
    params.filename_pattern = '*_outputs_*.mat';
end
if ~isfield(params, 'variable_name') || isempty(params.variable_name)
    params.variable_name = 'EDMD_outputs';
end
if ~isfield(params, 'n_top_windows') || isempty(params.n_top_windows)
    params.n_top_windows = 30;
end
if ~isfield(params, 'window_length_samples') || isempty(params.window_length_samples)
    params.window_length_samples = 6000;
end
if ~isfield(params, 'window_mode') || isempty(params.window_mode)
    params.window_mode = 'global';
end

if ~isfield(params, 'abs_thresh') || isempty(params.abs_thresh)
    params.abs_thresh = 0.01;
end
if ~isfield(params, 'sort_by') || isempty(params.sort_by)
    params.sort_by = 'modulus';
end
if ~isfield(params, 'sort_dir') || isempty(params.sort_dir)
    params.sort_dir = 'descend';
end
if ~isfield(params, 'max_basis') || isempty(params.max_basis)
    params.max_basis = 30;
end
if ~isfield(params, 'timescale_max_modes_all') || isempty(params.timescale_max_modes_all)
    params.timescale_max_modes_all = 30;
end
if ~isfield(params, 'timescale_max_modes_sel') || isempty(params.timescale_max_modes_sel)
    params.timescale_max_modes_sel = 20;
end
if ~isfield(params, 'timescale_max_lag') || isempty(params.timescale_max_lag)
    params.timescale_max_lag = 200;
end

if ~isfield(params, 'do_main_plot') || isempty(params.do_main_plot)
    params.do_main_plot = true;
end
if ~isfield(params, 'do_timescale') || isempty(params.do_timescale)
    params.do_timescale = true;
end
if ~isfield(params, 'do_deconv') || isempty(params.do_deconv)
    params.do_deconv = true;
end
if ~isfield(params, 'deconv_method') || isempty(params.deconv_method)
    params.deconv_method = 'koopman_residual';
end
if ~isfield(params, 'deconv_lambda_source') || isempty(params.deconv_lambda_source)
    params.deconv_lambda_source = 'edmd';
end
if ~isfield(params, 'deconv_plot_normalize_scope') || isempty(params.deconv_plot_normalize_scope)
    params.deconv_plot_normalize_scope = 'global';
end
if ~isfield(params, 'do_deconv_window_norm') || isempty(params.do_deconv_window_norm)
    params.do_deconv_window_norm = true;
end
if ~isfield(params, 'deconv_window_plot_normalize_scope') || isempty(params.deconv_window_plot_normalize_scope)
    params.deconv_window_plot_normalize_scope = 'window';
end
if ~isfield(params, 'deconv_max_modes_all') || isempty(params.deconv_max_modes_all)
    params.deconv_max_modes_all = 30;
end
if ~isfield(params, 'deconv_max_modes_sel') || isempty(params.deconv_max_modes_sel)
    params.deconv_max_modes_sel = 20;
end
if ~isfield(params, 'deconv_normalize_exclude_idx') || isempty(params.deconv_normalize_exclude_idx)
    params.deconv_normalize_exclude_idx = 1;
end

if ~isfield(params, 'do_spike_correlation') || isempty(params.do_spike_correlation)
    params.do_spike_correlation = false;
end
if ~isfield(params, 'spike_corr_channels') || isempty(params.spike_corr_channels)
    params.spike_corr_channels = 'all';
end
if ~isfield(params, 'spike_corr_max_modes') || isempty(params.spike_corr_max_modes)
    params.spike_corr_max_modes = 20;
end
if ~isfield(params, 'spike_corr_top_n_rows') || isempty(params.spike_corr_top_n_rows)
    params.spike_corr_top_n_rows = 100;
end
if ~isfield(params, 'spike_corr_progress_every') || isempty(params.spike_corr_progress_every)
    params.spike_corr_progress_every = 50;
end
if ~isfield(params, 'spike_corr_skip_existing') || isempty(params.spike_corr_skip_existing)
    params.spike_corr_skip_existing = true;
end
if ~isfield(params, 'spike_corr_save_overview') || isempty(params.spike_corr_save_overview)
    params.spike_corr_save_overview = true;
end
if ~isfield(params, 'spike_corr_overview_top_n') || isempty(params.spike_corr_overview_top_n)
    params.spike_corr_overview_top_n = 20;
end
if ~isfield(params, 'spike_corr_overview_features') || isempty(params.spike_corr_overview_features)
    params.spike_corr_overview_features = {'abs_mean', 'abs_rms'};
end
params.spike_corr_overview_features = cellstr(string(params.spike_corr_overview_features(:)).');
if ~isfield(params, 'spike_corr_overview_resolution') || isempty(params.spike_corr_overview_resolution)
    params.spike_corr_overview_resolution = 220;
end
if ~isfield(params, 'spike_corr_verbose') || isempty(params.spike_corr_verbose)
    params.spike_corr_verbose = true;
end
if ~isfield(params, 'spike_corr_save_under_run') || isempty(params.spike_corr_save_under_run)
    params.spike_corr_save_under_run = false;
end
if ~isfield(params, 'spike_corr_save_dir_name') || isempty(params.spike_corr_save_dir_name)
    params.spike_corr_save_dir_name = 'mlp_fulltime_spike_residual_correlation';
end

if ~isfield(params, 'save_png') || isempty(params.save_png)
    params.save_png = true;
end
if ~isfield(params, 'save_fig') || isempty(params.save_fig)
    params.save_fig = false;
end
if ~isfield(params, 'skip_existing') || isempty(params.skip_existing)
    params.skip_existing = true;
end
if ~isfield(params, 'close_figures') || isempty(params.close_figures)
    params.close_figures = true;
end
if ~isfield(params, 'headless') || isempty(params.headless)
    params.headless = true;
end
if ~isfield(params, 'resolution') || isempty(params.resolution)
    params.resolution = 180;
end
end


function dataset_output_root = local_dataset_output_root(params, dataset)
if isempty(params.output_root)
    dataset_output_root = fullfile(params.processed_root, dataset, ...
        'postprocessing', params.output_folder_name);
else
    dataset_output_root = fullfile(params.output_root, dataset);
end
end


function runs = local_discover_completed_runs(params)
runs = struct('dataset', {}, 'observable_mode', {}, 'residual_form', {}, ...
    'run_name', {}, 'output_dir', {}, 'summary_count', {}, 'chunk_count', {}, ...
    'last_write_time', {});

for i = 1:numel(params.dataset_stems)
    dataset = params.dataset_stems{i};
    if any(strcmpi(dataset, params.exclude_dataset_stems))
        continue;
    end

    output_parent = fullfile(params.autodl_root, dataset, 'mlp', 'outputs');
    if exist(output_parent, 'dir') ~= 7
        continue;
    end

    L = dir(output_parent);
    L = L([L.isdir]);
    L = L(~ismember({L.name}, {'.', '..'}));

    for k = 1:numel(L)
        run_name = L(k).name;
        if contains(lower(run_name), 'smoke')
            continue;
        end

        output_dir = fullfile(L(k).folder, run_name);
        summary_files = dir(fullfile(output_dir, '*_summary.mat'));
        chunk_files = local_collect_chunk_files(output_dir, params.filename_pattern);
        if isempty(summary_files) || isempty(chunk_files)
            continue;
        end

        [observable_mode, residual_form] = local_infer_condition_from_run_name(run_name);
        file_times = [summary_files.datenum, chunk_files.datenum];
        last_write = max(file_times);

        run = struct();
        run.dataset = dataset;
        run.observable_mode = observable_mode;
        run.residual_form = residual_form;
        run.run_name = run_name;
        run.output_dir = output_dir;
        run.summary_count = numel(summary_files);
        run.chunk_count = numel(chunk_files);
        run.last_write_time = last_write;
        runs(end+1) = run; %#ok<AGROW>
    end
end

if isempty(runs)
    return;
end

% Deduplicate repeated local runs of the same condition, keeping newest.
keys = strings(numel(runs), 1);
for i = 1:numel(runs)
    keys(i) = string(lower(strjoin({runs(i).dataset, ...
        runs(i).observable_mode, runs(i).residual_form}, '|')));
end

keep = false(numel(runs), 1);
u = unique(keys, 'stable');
for i = 1:numel(u)
    idx = find(keys == u(i));
    [~, best_local] = max([runs(idx).last_write_time]);
    keep(idx(best_local)) = true;
end
runs = runs(keep);

runs = local_filter_discovered_runs(runs, params);
if isempty(runs)
    return;
end

[~, order] = sort(string({runs.dataset}) + "|" + string({runs.observable_mode}) + "|" + string({runs.residual_form}));
runs = runs(order);
end


function runs = local_filter_discovered_runs(runs, params)
if isempty(runs)
    return;
end

if isfield(params, 'run_name_filter') && ~isempty(params.run_name_filter)
    keep = false(numel(runs), 1);
    for i = 1:numel(runs)
        keep(i) = local_matches_filter(runs(i).run_name, params.run_name_filter);
    end
    runs = runs(keep);
end

if isempty(runs)
    return;
end

if isfield(params, 'condition_key_filter') && ~isempty(params.condition_key_filter)
    keep = false(numel(runs), 1);
    for i = 1:numel(runs)
        keep(i) = local_matches_filter(local_condition_key(runs(i)), ...
            params.condition_key_filter);
    end
    runs = runs(keep);
end
end


function tf = local_matches_filter(value, filter_values)
tf = any(strcmpi(string(value), string(filter_values)));
end


function key = local_condition_key(run_info)
key = lower(strjoin({run_info.dataset, run_info.observable_mode, ...
    run_info.residual_form}, '|'));
end


function [observable_mode, residual_form] = local_infer_condition_from_run_name(run_name)
name_lower = lower(run_name);
if contains(name_lower, 'complex_split') || contains(name_lower, 'complexsplit')
    observable_mode = 'complex_split';
elseif contains(name_lower, '_abs') || endsWith(name_lower, 'abs')
    observable_mode = 'abs';
else
    observable_mode = 'unknown';
end

if contains(name_lower, 'projected_vlambda')
    residual_form = 'projected_vlambda';
elseif contains(name_lower, 'projected_kv')
    residual_form = 'projected_kv';
else
    residual_form = 'unknown';
end
end


function cfg = local_dataset_cfg(dataset)
switch lower(dataset)
    case 'e10fv1'
        cfg = cfg_E10fV1();
    case 'e10gb1'
        cfg = cfg_E10gb1();
    case 'e10gh1'
        cfg = cfg_E10gH1();
    case 'f12m01'
        cfg = cfg_F12m01();
    otherwise
        error('No dataset config mapping is defined for %s.', dataset);
end

if ~isfield(cfg, 'spectrogram')
    cfg.spectrogram = struct();
end
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
end


function cache = local_load_top_windows_for_dataset(cfg, params)
window_file = fullfile(params.processed_root, cfg.file_stem, ...
    'consensus_state_diversity_windows', ...
    sprintf('%s_consensus_state_diversity_windows_%dsamp_globalwin.mat', ...
    cfg.file_stem, params.window_length_samples));

if exist(window_file, 'file') ~= 2
    fprintf('State-diversity window file missing; computing it:\n  %s\n', window_file);
    loader_cfg = struct('file_stem', cfg.file_stem);
    [C, source_consensus_file] = io_results.load_consensus_state_results(loader_cfg, params.processed_root, []);
    win_params = struct();
    win_params.window_length_samples = params.window_length_samples;
    win_params.window_mode = params.window_mode;
    win_params.keep_partial_window = false;
    win_params.top_k = params.n_top_windows;
    win_params.save_csv = true;
    win_params.force_recompute = false;
    W = analyze_blp_consensus_state_diversity_windows( ...
        cfg, params.processed_root, C, win_params, source_consensus_file);
    window_file = W.save_file;
else
    S = load(window_file, 'W');
    if ~isfield(S, 'W')
        error('State-diversity file does not contain variable W: %s', window_file);
    end
    W = S.W;
end

local_validate_top_window_result(W, window_file, params);
cache = struct('W', W, 'window_file', window_file);
end


function local_validate_top_window_result(W, window_file, params)
if ~isfield(W, 'top_windows_table') || isempty(W.top_windows_table)
    error('Top-window table is missing or empty: %s', window_file);
end
if ~isfield(W, 'window_mode') || ~strcmpi(char(string(W.window_mode)), params.window_mode)
    error('Top-window file is not a %s-window result: %s', params.window_mode, window_file);
end
if ~isfield(W, 'window_length_samples') || ...
        double(W.window_length_samples) ~= double(params.window_length_samples)
    error('Top-window file does not use %d-sample windows: %s', ...
        params.window_length_samples, window_file);
end

required_vars = {'state_diversity_rank', 'global_start_idx', 'global_end_idx'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table is missing required column "%s": %s', ...
            required_vars{i}, window_file);
    end
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


function EDMD_window = local_load_edmd_window(chunk_index, idx1, idx2, variable_name)
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


function [fig_ts, timescale_info] = local_run_timescale_postprocess(EDMD_post, title_suffix, params)
cfg_ts = struct();
cfg_ts.t_plot = 1:size(EDMD_post.efuns, 1);
cfg_ts.max_modes_all = min(params.timescale_max_modes_all, ...
    size(EDMD_post.original_sorted.efuns, 2));
cfg_ts.max_modes_sel = min(params.timescale_max_modes_sel, ...
    numel(EDMD_post.evalues));
cfg_ts.maxLag = min(params.timescale_max_lag, size(EDMD_post.efuns, 1) - 1);
cfg_ts.title_prefix = ['EDMD timescale diagnostics | ', title_suffix];

[fig_ts, timescale_info] = postprocess_EDMD_outputs_timescale(EDMD_post, cfg_ts);
if ~isempty(fig_ts) && isgraphics(fig_ts)
    set(fig_ts, 'Name', ['EDMD timescale | ', title_suffix], ...
        'NumberTitle', 'off');
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
cfg_deconv.max_modes_all = min(params.deconv_max_modes_all, numel(EDMD_post.evalues));
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
    set(fig_deconv, 'Name', ['EDMD deconv | ', title_suffix], ...
        'NumberTitle', 'off');
end
end


function spike_row = local_run_spike_residual_correlation_fulltime( ...
    cfg, run_info, run_output_root, params)

out = struct();
out.dataset = string(run_info.dataset);
out.observable_mode = string(run_info.observable_mode);
out.residual_form = string(run_info.residual_form);
out.run_name = string(run_info.run_name);
out.run_output_dir = string(run_info.output_dir);
out.spike_corr_dir = "";
out.spike_result_mat = "";
out.pooled_corr_csv = "";
out.pooled_top_csv = "";
out.session_corr_csv = "";
out.mode_csv = "";
out.session_summary_csv = "";
out.overview_png = "";
out.overview_fig = "";
out.n_pooled_corr_rows = NaN;
out.n_pooled_top_rows = NaN;
out.top_abs_corr = NaN;
out.top_corr = NaN;
out.top_channel_label = "";
out.top_residual_feature = "";
out.top_mode_rank = NaN;
out.status = "ok";
out.error_message = "";

save_dir = local_spike_corr_output_root(params, run_info, run_output_root);
out.spike_corr_dir = string(save_dir);

try
    if exist(save_dir, 'dir') ~= 7
        mkdir(save_dir);
    end

    fprintf('  Full-time spike-residual correlation:\n    %s\n', save_dir);

    existing_mat = '';
    if params.spike_corr_skip_existing
        existing_mat = local_find_existing_spike_result(save_dir, cfg, run_info, params);
    end

    if ~isempty(existing_mat)
        fprintf('    Using existing result:\n      %s\n', existing_mat);
        try
            S = load(existing_mat, 'result');
            if ~isfield(S, 'result')
                error('Existing spike-correlation MAT lacks variable result: %s', existing_mat);
            end
            result = local_ensure_spike_save_paths(S.result, existing_mat);
        catch ME_load
            warning('Ignoring unreadable existing spike-correlation MAT and recomputing: %s', ...
                ME_load.message);
            existing_mat = '';
        end
    end

    if isempty(existing_mat)
        cmp_params = struct();
        cmp_params.source_cfg = struct();
        cmp_params.source_cfg.mode = 'chunk_dir';
        cmp_params.source_cfg.data_dir = run_info.output_dir;
        cmp_params.source = struct();
        cmp_params.source.base_dir = fileparts(run_info.output_dir);
        cmp_params.source.name_contains = run_info.run_name;
        cmp_params.source.prefer_non_smoke = true;
        cmp_params.filename_pattern = params.filename_pattern;
        cmp_params.variable_name = params.variable_name;
        cmp_params.post = struct();
        cmp_params.post.abs_thresh = params.abs_thresh;
        cmp_params.post.sort_by = params.sort_by;
        cmp_params.post.sort_dir = params.sort_dir;
        cmp_params.post.max_basis = params.max_basis;
        cmp_params.residual = struct();
        cmp_params.residual.method = params.deconv_method;
        cmp_params.residual.lambda_source = params.deconv_lambda_source;
        cmp_params.residual.lambdaType = 'discrete';
        cmp_params.residual.first_u_mode = 'phi1';
        cmp_params.residual.max_modes = params.spike_corr_max_modes;
        cmp_params.spike_channels = params.spike_corr_channels;
        cmp_params.output_root = params.processed_root;
        cmp_params.save_dir = save_dir;
        cmp_params.save_results = true;
        cmp_params.verbose = params.spike_corr_verbose;
        cmp_params.progress_every = params.spike_corr_progress_every;
        cmp_params.top_n_rows = params.spike_corr_top_n_rows;

        result = compute_spike_residual_comparison(cfg, cmp_params);
        result = local_ensure_spike_save_paths(result, result.save_paths.main_mat);
    end

    out.spike_result_mat = string(result.save_paths.main_mat);
    out.pooled_corr_csv = string(result.save_paths.pooled_corr_csv);
    out.pooled_top_csv = string(result.save_paths.pooled_top_csv);
    out.session_corr_csv = string(result.save_paths.session_corr_csv);
    out.mode_csv = string(result.save_paths.mode_csv);
    out.session_summary_csv = string(result.save_paths.session_summary_csv);
    out.n_pooled_corr_rows = height(result.pooled_corr_table);
    out.n_pooled_top_rows = height(result.pooled_top_corr);

    if ~isempty(result.pooled_top_corr)
        top_row = result.pooled_top_corr(1, :);
        out.top_abs_corr = double(top_row.abs_corr);
        out.top_corr = double(top_row.corr);
        out.top_channel_label = string(top_row.channel_label);
        out.top_residual_feature = string(top_row.residual_feature);
        out.top_mode_rank = double(top_row.mode_rank);
    end

    if params.spike_corr_save_overview
        [overview_png, overview_fig] = local_save_spike_corr_overview( ...
            result, cfg, run_info, save_dir, params);
        out.overview_png = string(overview_png);
        out.overview_fig = string(overview_fig);
    end
catch ME
    out.status = "failed";
    out.error_message = string(local_single_line_error(ME));
    warning('Full-time spike-residual correlation failed for %s: %s', ...
        run_info.run_name, out.error_message);
    close all force;
end

spike_row = struct2table(out, 'AsArray', true);
end


function save_dir = local_spike_corr_output_root(params, run_info, run_output_root)
if params.spike_corr_save_under_run
    save_dir = fullfile(run_output_root, params.spike_corr_save_dir_name);
else
    save_dir = fullfile(params.processed_root, run_info.dataset, ...
        'postprocessing', params.spike_corr_save_dir_name);
end
end


function result = local_ensure_spike_save_paths(result, main_mat)
if ~isfield(result, 'save_paths') || isempty(result.save_paths) || ...
        ~isfield(result.save_paths, 'main_mat') || isempty(result.save_paths.main_mat)
    [save_dir, stem] = fileparts(main_mat);
    result.save_paths = struct();
    result.save_paths.main_mat = main_mat;
    result.save_paths.session_corr_csv = fullfile(save_dir, [stem, '_session_corr.csv']);
    result.save_paths.pooled_corr_csv = fullfile(save_dir, [stem, '_pooled_corr.csv']);
    result.save_paths.pooled_top_csv = fullfile(save_dir, [stem, '_pooled_top.csv']);
    result.save_paths.mode_csv = fullfile(save_dir, [stem, '_modes.csv']);
    result.save_paths.session_summary_csv = fullfile(save_dir, [stem, '_session_summary.csv']);
end
end


function existing_mat = local_find_existing_spike_result(save_dir, cfg, run_info, params)
existing_mat = '';
run_label = local_filename_safe(run_info.run_name);
spike_label = local_spike_channel_label(params.spike_corr_channels);
tag = sprintf('%s_%s_%s_koopman_residual_edmd_top%d', ...
    cfg.file_stem, run_label, spike_label, params.spike_corr_max_modes);

candidate = fullfile(save_dir, [tag, '.mat']);
if exist(candidate, 'file') == 2
    D = dir(candidate);
    if ~isempty(D) && D.bytes > 0
        existing_mat = candidate;
        return;
    end
end

L = dir(fullfile(save_dir, sprintf('%s_%s*_koopman_residual_edmd_top*.mat', ...
    cfg.file_stem, run_label)));
if isempty(L)
    return;
end

names = {L.name};
keep = ~contains(names, '_lagged', 'IgnoreCase', true);
keep = keep & [L.bytes] > 0;
L = L(keep);
if isempty(L)
    return;
end

[~, order] = sort([L.datenum], 'descend');
L = L(order);
existing_mat = fullfile(L(1).folder, L(1).name);
end


function label = local_spike_channel_label(spike_channel_spec)
if ischar(spike_channel_spec) || isstring(spike_channel_spec)
    spec = lower(strtrim(char(string(spike_channel_spec))));
    switch spec
        case 'selected'
            label = 'spk_selected';
            return;
        case 'all'
            label = 'spk_all';
            return;
    end
end

if isnumeric(spike_channel_spec)
    label = sprintf('spk_custom_%dch', numel(spike_channel_spec));
else
    label = 'spk_custom';
end
end


function [overview_png, overview_fig] = local_save_spike_corr_overview( ...
    result, cfg, run_info, save_dir, params)

[~, result_stem] = fileparts(result.save_paths.main_mat);
overview_png = fullfile(save_dir, [result_stem, '_overview.png']);
overview_fig = fullfile(save_dir, [result_stem, '_overview.fig']);

if params.spike_corr_skip_existing && exist(overview_png, 'file') == 2
    return;
end

fig = local_plot_spike_corr_overview(result, cfg, run_info, params);
drawnow;
exportgraphics(fig, overview_png, 'Resolution', params.spike_corr_overview_resolution);
if params.save_fig
    savefig(fig, overview_fig);
else
    overview_fig = '';
end
local_close_fig(fig, params);
end


function fig = local_plot_spike_corr_overview(result, cfg, run_info, params)
T = result.pooled_corr_table;
if isempty(T)
    error('pooled_corr_table is empty for %s.', run_info.run_name);
end

features_to_plot = params.spike_corr_overview_features;
if numel(features_to_plot) < 2
    features_to_plot = [features_to_plot, features_to_plot(1)];
end
features_to_plot = features_to_plot(1:2);

channel_labels = cellstr(string(result.channel_table.channel_label));
mode_labels = compose('m%d', double(result.mode_table.mode_rank));
n_channels = height(result.channel_table);
n_modes = height(result.mode_table);

plot_mats = cell(numel(features_to_plot), 1);
for i = 1:numel(features_to_plot)
    plot_mats{i} = local_build_spike_corr_matrix( ...
        T, features_to_plot{i}, result.channel_table, result.mode_table);
end

all_vals = [];
for i = 1:numel(plot_mats)
    vals = plot_mats{i};
    all_vals = [all_vals; vals(isfinite(vals))]; %#ok<AGROW>
end

if isempty(all_vals)
    clim_use = [-1, 1];
else
    cmax = max(abs(all_vals));
    if cmax == 0
        cmax = 1;
    end
    clim_use = [-cmax, cmax];
end

top_n = min(params.spike_corr_overview_top_n, height(T));
top_tbl = T(1:top_n, :);
top_labels = strings(top_n, 1);
for i = 1:top_n
    top_labels(i) = sprintf('%s | %s | m%d', ...
        char(top_tbl.channel_label(i)), ...
        char(top_tbl.residual_feature(i)), ...
        double(top_tbl.mode_rank(i)));
end

fig = figure('Color', 'w', 'Position', [100, 100, 1500, 900]);
tiledlayout(2, 2, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:numel(features_to_plot)
    nexttile(i);
    imagesc(plot_mats{i});
    axis tight;
    set(gca, 'YTick', 1:n_channels, 'YTickLabel', channel_labels, ...
        'XTick', 1:n_modes, 'XTickLabel', mode_labels, ...
        'TickLabelInterpreter', 'none');
    xtickangle(0);
    xlabel('Residual Mode Rank');
    ylabel('Spike Channel');
    title(sprintf('Full-time pooled corr: %s', ...
        strrep(features_to_plot{i}, '_', '\_')));
    colormap(gca, parula(256));
    colorbar;
    caxis(clim_use);
end

nexttile([1, 2]);
if top_n > 0
    barh(1:top_n, top_tbl.corr, ...
        'FaceColor', [0.18 0.46 0.71], 'EdgeColor', 'none');
    set(gca, 'YTick', 1:top_n, 'YTickLabel', top_labels, ...
        'TickLabelInterpreter', 'none', 'YDir', 'reverse');
else
    text(0.5, 0.5, 'No finite pooled correlations found', ...
        'Units', 'normalized', 'HorizontalAlignment', 'center');
end
xlabel('Correlation');
ylabel('Top Pairs');
title(sprintf('Top %d full-time pooled spike-residual correlations', top_n));
grid on;
xline(0, 'k-');

sgtitle(sprintf('%s | %s | full-time spike vs deconv residual correlation', ...
    cfg.file_stem, run_info.run_name), 'FontWeight', 'bold', 'Interpreter', 'none');
end


function M = local_build_spike_corr_matrix(T, feature_name, channel_table, mode_table)
n_channels = height(channel_table);
n_modes = height(mode_table);
M = nan(n_channels, n_modes);

mask = strcmp(string(T.residual_feature), string(feature_name));
Tf = T(mask, :);

for i = 1:height(Tf)
    ch_idx = find(channel_table.channel_index == Tf.channel_index(i), 1, 'first');
    mode_idx = find(mode_table.mode_rank == Tf.mode_rank(i), 1, 'first');
    if isempty(ch_idx) || isempty(mode_idx)
        continue;
    end
    M(ch_idx, mode_idx) = Tf.corr(i);
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


function out = local_filename_safe(name_in)
out = regexprep(char(string(name_in)), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled';
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
if params.do_timescale
    required_files(end+1, 1) = out.timescale_png;
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
