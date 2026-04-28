% Plot e10gb1 eigenfunction/SVD results in the top consensus-state-diversity windows.
%
% This script intentionally reuses the existing eigenfunction plotting
% functions. It only handles loading window slices from the compact result
% and the original EDMD chunk files.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
    % Older MATLAB releases do not expose this graphics default.
end
set(groot, 'defaultFigureVisible', 'off');

if ~exist('cfg', 'var') || ~isstruct(cfg) || ~isfield(cfg, 'source')
    cfg = cfg_eigenfunction_reduction_minimal();
end
output_root = cfg.dataset.processed_root;

params = struct();
params.n_top_windows = 30;
params.make_overview = true;                 % Eigenfunction heatmap + temporal components.
params.make_state_space = false;             % Optional trajectory colored by dim/time.
params.make_consensus_state_space = true;    % Same trajectory colored by consensus state.
params.close_figures = true;
params.force_recompute_norm_cache = false;
params.chunk_scan_progress_every = 100;
params.baseline_code = cfg.viz.state_space_consensus.baseline_code;
params.state_start_idx = cfg.viz.state_space_consensus.state_start_idx;
params.save_dir = fullfile(cfg.output.root, 'top30');
params.norm_cache_file = fullfile(params.save_dir, 'norm.mat');

if exist('top_window_params', 'var') && isstruct(top_window_params)
    params = local_merge_struct(params, top_window_params);
    if ~isfield(top_window_params, 'baseline_code') || isempty(top_window_params.baseline_code)
        params.baseline_code = cfg.viz.state_space_consensus.baseline_code;
    end
    if ~isfield(top_window_params, 'state_start_idx') || isempty(top_window_params.state_start_idx)
        params.state_start_idx = cfg.viz.state_space_consensus.state_start_idx;
    end
    if ~isfield(top_window_params, 'save_dir') || isempty(top_window_params.save_dir)
        params.save_dir = fullfile(cfg.output.root, 'top30');
    end
    if ~isfield(top_window_params, 'norm_cache_file') || isempty(top_window_params.norm_cache_file)
        params.norm_cache_file = fullfile(params.save_dir, 'norm.mat');
    end
end

if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

if ~exist('result_file', 'var') || isempty(result_file)
    result_file = local_find_latest_result_file(cfg.save.dir);
end
window_result_file = fullfile(output_root, cfg.dataset.name, ...
    'consensus_state_diversity_windows', ...
    sprintf('%s_consensus_state_diversity_windows_6000samp_globalwin.mat', ...
    cfg.dataset.name));

fprintf('Loading compact eigenfunction reduction result:\n  %s\n', result_file);
S = load(result_file, 'result');
if ~isfield(S, 'result')
    error('Result file %s does not contain variable result.', result_file);
end
result = S.result;

fprintf('Loading top consensus-state-diversity windows:\n  %s\n', window_result_file);
S = load(window_result_file, 'W');
if ~isfield(S, 'W')
    error('Top-window file %s does not contain variable W.', window_result_file);
end
W = S.W;
local_validate_state_diversity_result(W, window_result_file);

top_windows = W.top_windows_table;
n_windows = min(params.n_top_windows, height(top_windows));
top_windows = top_windows(1:n_windows, :);

consensus_loader_cfg = struct();
consensus_loader_cfg.file_stem = cfg.dataset.name;
[C_consensus, source_consensus_file] = io_results.load_consensus_state_results( ...
    consensus_loader_cfg, output_root, []);
fprintf('Loading consensus states:\n  %s\n', source_consensus_file);

selected_mode_idx = local_get_selected_mode_idx(result);
chunk_files = local_collect_edmd_chunk_files( ...
    cfg.source.data_dir, cfg.source.concat.filename_pattern);
[norm_denominator, chunk_table] = local_get_or_compute_norm_denominator( ...
    chunk_files, selected_mode_idx, cfg, params);

T_result = local_get_result_length(result);
T_chunks = chunk_table.chunk_end_idx(end);
if T_result ~= T_chunks
    warning('Result length (%d) differs from EDMD chunk length (%d).', ...
        T_result, T_chunks);
end

overview_plot = strings(n_windows, 1);
state_space_plot = strings(n_windows, 1);
consensus_state_space_plot = strings(n_windows, 1);

fprintf('Saving top-window eigenfunction figures to:\n  %s\n', params.save_dir);

for i = 1:n_windows
    idx1 = double(top_windows.global_start_idx(i));
    idx2 = double(top_windows.global_end_idx(i));
    local_validate_window_span(idx1, idx2, T_result, i);

    file_stub = local_build_file_stub(top_windows, i);
    title_text = local_build_window_title(top_windows, i, idx1, idx2);
    window_idx = (idx1:idx2).';

    fprintf('[%02d/%02d] %s | samples [%d, %d]\n', ...
        i, n_windows, file_stub, idx1, idx2);

    if params.make_overview
        efun_feature_window = local_load_efun_feature_window( ...
            chunk_files, chunk_table, selected_mode_idx, norm_denominator, ...
            idx1, idx2, result.meta.feature_variant);
        result_window = local_build_window_result( ...
            result, efun_feature_window, idx1, idx2);

        overview_cfg = cfg.viz.overview;
        overview_cfg.window_idx = [];
        overview_cfg.save_figure = false;
        overview_cfg.figure_visible = 'off';
        overview_cfg.save_dir = params.save_dir;
        overview_cfg.save_tag = [file_stub, '_ov'];
        overview_cfg.title = ['Eigenfunctions And Temporal Components | ', title_text];
        overview_cfg.event_windows = local_build_consensus_event_windows( ...
            C_consensus, result_window.input.time_axis, idx1, idx2, params);
        overview_cfg.event_colors = local_build_consensus_event_colors(C_consensus, cfg);

        [fig_overview, overview_info] = ...
            plot_eigenfunction_component_overview(result_window, overview_cfg);
        overview_info.save_path = local_export_existing_plot( ...
            fig_overview, params.save_dir, 'efun_comp', ...
            result, overview_cfg.save_tag, 200);
        overview_plot(i) = string(overview_info.save_path);
        local_close_if_requested(fig_overview, params);
    end

    if params.make_state_space
        state_cfg = cfg.viz.state_space;
        state_cfg.window_idx = window_idx;
        state_cfg.save_figure = false;
        state_cfg.figure_visible = 'off';
        state_cfg.save_dir = params.save_dir;
        state_cfg.save_tag = [file_stub, '_ss'];
        state_cfg.title = ['State-Space Trajectory | ', title_text];

        [fig_state_space, state_info] = ...
            plot_eigenfunction_state_space_trajectory(result, state_cfg);
        state_info.save_path = local_export_existing_plot( ...
            fig_state_space, params.save_dir, 'efun_ss', ...
            result, state_cfg.save_tag, 220);
        state_space_plot(i) = string(state_info.save_path);
        local_close_if_requested(fig_state_space, params);
    end

    if params.make_consensus_state_space
        consensus_cfg = cfg.viz.state_space_consensus;
        consensus_cfg.window_idx = window_idx;
        consensus_cfg.save_figure = false;
        consensus_cfg.figure_visible = 'off';
        consensus_cfg.save_dir = params.save_dir;
        consensus_cfg.save_tag = [file_stub, '_ssc'];
        consensus_cfg.title = ['State-Space Trajectory By Consensus State | ', title_text];

        [fig_consensus_state_space, consensus_state_info] = ...
            plot_eigenfunction_state_space_consensus_trajectory( ...
            result, C_consensus, consensus_cfg);
        consensus_state_info.save_path = local_export_existing_plot( ...
            fig_consensus_state_space, params.save_dir, ...
            'efun_ssc', result, ...
            consensus_cfg.save_tag, 240);
        consensus_state_space_plot(i) = string(consensus_state_info.save_path);
        local_close_if_requested(fig_consensus_state_space, params);
    end
end

manifest_table = top_windows;
manifest_table.source_result_file = repmat(string(result_file), n_windows, 1);
manifest_table.source_state_diversity_file = repmat(string(window_result_file), n_windows, 1);
manifest_table.source_consensus_file = repmat(string(source_consensus_file), n_windows, 1);
manifest_table.eigenfunction_overview_png = overview_plot;
if params.make_state_space
    manifest_table.state_space_png = state_space_plot;
end
manifest_table.consensus_state_space_png = consensus_state_space_plot;

manifest_file = fullfile(params.save_dir, 'manifest.csv');
writetable(manifest_table, manifest_file);

fprintf('\nSaved %d top-window plot sets.\n', n_windows);
fprintf('Manifest:\n  %s\n', manifest_file);


function result_file = local_find_latest_result_file(result_dir)
if exist(result_dir, 'dir') ~= 7
    error(['Result directory does not exist:\n  %s\n', ...
        'Run script_run_eigenfunction_reduction_minimal.m first.'], result_dir);
end

L = dir(fullfile(result_dir, '*.mat'));
if isempty(L)
    error(['No eigenfunction reduction MAT files were found in:\n  %s\n', ...
        'Run script_run_eigenfunction_reduction_minimal.m first.'], result_dir);
end

[~, idx] = max([L.datenum]);
result_file = fullfile(L(idx).folder, L(idx).name);
end


function params = local_merge_struct(params, override)
names = fieldnames(override);
for i = 1:numel(names)
    params.(names{i}) = override.(names{i});
end
end


function local_validate_state_diversity_result(W, window_result_file)
if ~isfield(W, 'top_windows_table') || isempty(W.top_windows_table)
    error('Top-window table is empty or missing in %s.', window_result_file);
end

if ~isfield(W, 'window_mode') || ~strcmpi(char(string(W.window_mode)), 'global')
    error('%s is not a global-window state-diversity result.', window_result_file);
end

if ~isfield(W, 'window_length_samples') || double(W.window_length_samples) ~= 6000
    error('%s does not use 6000-sample windows.', window_result_file);
end

required_vars = {'state_diversity_rank', 'global_window_idx', ...
    'global_start_idx', 'global_end_idx'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table in %s is missing required column "%s".', ...
            window_result_file, required_vars{i});
    end
end
end


function selected_mode_idx = local_get_selected_mode_idx(result)
if ~isfield(result, 'input') || ...
        ~isfield(result.input, 'selected_mode_idx_in_original') || ...
        isempty(result.input.selected_mode_idx_in_original)
    error('result.input.selected_mode_idx_in_original is required.');
end

selected_mode_idx = double(result.input.selected_mode_idx_in_original(:)).';
end


function files = local_collect_edmd_chunk_files(data_dir, filename_pattern)
if exist(data_dir, 'dir') ~= 7
    error('EDMD chunk directory does not exist:\n  %s', data_dir);
end

L = dir(fullfile(data_dir, filename_pattern));
files = repmat(struct('name', '', 'fullpath', '', 'chunk_id', NaN), numel(L), 1);
n_keep = 0;

for i = 1:numel(L)
    tokens = regexp(L(i).name, '^(.*)_outputs_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tokens)
        continue;
    end

    n_keep = n_keep + 1;
    files(n_keep).name = L(i).name;
    files(n_keep).fullpath = fullfile(L(i).folder, L(i).name);
    files(n_keep).chunk_id = str2double(tokens{2});
end

files = files(1:n_keep);
if isempty(files)
    error('No EDMD chunk files matching %s were found in %s.', ...
        filename_pattern, data_dir);
end

[~, order] = sort([files.chunk_id]);
files = files(order);
end


function [denominator, chunk_table] = local_get_or_compute_norm_denominator( ...
    files, selected_mode_idx, cfg, params)
cache_file = params.norm_cache_file;

if ~params.force_recompute_norm_cache && exist(cache_file, 'file') == 2
    S = load(cache_file, 'norm_cache');
    if isfield(S, 'norm_cache') && ...
            local_norm_cache_matches(S.norm_cache, files, selected_mode_idx, cfg)
        fprintf('Using cached eigenfunction normalization denominator:\n  %s\n', ...
            cache_file);
        denominator = double(S.norm_cache.denominator_max_abs_per_mode(:)).';
        chunk_table = S.norm_cache.chunk_table;
        return;
    end
end

fprintf('Scanning %d EDMD chunks for selected-mode normalization.\n', numel(files));

n_files = numel(files);
n_modes = numel(selected_mode_idx);
denominator = zeros(1, n_modes);
chunk_id = zeros(n_files, 1);
chunk_start_idx = zeros(n_files, 1);
chunk_end_idx = zeros(n_files, 1);
chunk_length = zeros(n_files, 1);
chunk_file = strings(n_files, 1);
cursor = 1;

for k = 1:n_files
    if k == 1 || mod(k, params.chunk_scan_progress_every) == 0 || k == n_files
        fprintf('[norm scan %d/%d] %s\n', k, n_files, files(k).name);
    end

    S = load(files(k).fullpath, cfg.source.concat.variable_name);
    if ~isfield(S, cfg.source.concat.variable_name)
        error('File %s does not contain variable %s.', ...
            files(k).fullpath, cfg.source.concat.variable_name);
    end
    EDMD_outputs = S.(cfg.source.concat.variable_name);
    if ~isfield(EDMD_outputs, 'efuns')
        error('EDMD_outputs.efuns is missing from %s.', files(k).fullpath);
    end
    if max(selected_mode_idx) > size(EDMD_outputs.efuns, 2)
        error('Selected mode index exceeds efuns columns in %s.', files(k).fullpath);
    end

    efuns_selected = EDMD_outputs.efuns(:, selected_mode_idx);
    denominator = max(denominator, double(max(abs(efuns_selected), [], 1)));

    chunk_id(k) = files(k).chunk_id;
    chunk_start_idx(k) = cursor;
    chunk_length(k) = size(EDMD_outputs.efuns, 1);
    chunk_end_idx(k) = cursor + chunk_length(k) - 1;
    chunk_file(k) = string(files(k).fullpath);
    cursor = chunk_end_idx(k) + 1;
end

chunk_table = table(chunk_id, chunk_start_idx, chunk_end_idx, ...
    chunk_length, chunk_file);

norm_cache = struct();
norm_cache.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
norm_cache.source_data_dir = cfg.source.data_dir;
norm_cache.filename_pattern = cfg.source.concat.filename_pattern;
norm_cache.n_chunks = numel(files);
norm_cache.first_chunk_file = files(1).name;
norm_cache.last_chunk_file = files(end).name;
norm_cache.selected_mode_idx_in_original = selected_mode_idx(:);
norm_cache.denominator_max_abs_per_mode = denominator(:);
norm_cache.chunk_table = chunk_table;

save(cache_file, 'norm_cache', '-v7.3');
fprintf('Saved eigenfunction normalization denominator cache:\n  %s\n', cache_file);
end


function tf = local_norm_cache_matches(norm_cache, files, selected_mode_idx, cfg)
tf = false;

required_fields = {'source_data_dir', 'filename_pattern', 'n_chunks', ...
    'first_chunk_file', 'last_chunk_file', 'selected_mode_idx_in_original', ...
    'denominator_max_abs_per_mode', 'chunk_table'};
for i = 1:numel(required_fields)
    if ~isfield(norm_cache, required_fields{i})
        return;
    end
end

tf = strcmp(char(norm_cache.source_data_dir), char(cfg.source.data_dir)) && ...
    strcmp(char(norm_cache.filename_pattern), char(cfg.source.concat.filename_pattern)) && ...
    double(norm_cache.n_chunks) == numel(files) && ...
    strcmp(char(norm_cache.first_chunk_file), files(1).name) && ...
    strcmp(char(norm_cache.last_chunk_file), files(end).name) && ...
    isequal(double(norm_cache.selected_mode_idx_in_original(:)).', ...
    double(selected_mode_idx(:)).');
end


function T = local_get_result_length(result)
if ~isfield(result, 'core') || ...
        ~isfield(result.core, 'temporal_components_time_by_comp')
    error('result.core.temporal_components_time_by_comp is required.');
end

T = size(result.core.temporal_components_time_by_comp, 1);
end


function local_validate_window_span(idx1, idx2, T, row_idx)
if idx1 < 1 || idx2 > T || idx2 < idx1
    error('Invalid global sample range [%d, %d] for top-window row %d.', ...
        idx1, idx2, row_idx);
end
end


function X_feature = local_load_efun_feature_window(files, chunk_table, ...
    selected_mode_idx, denominator, idx1, idx2, feature_variant)
window_len = idx2 - idx1 + 1;
n_modes = numel(selected_mode_idx);
X_feature = zeros(window_len, n_modes);

overlap = find(chunk_table.chunk_end_idx >= idx1 & chunk_table.chunk_start_idx <= idx2);
if isempty(overlap)
    error('No EDMD chunks overlap sample range [%d, %d].', idx1, idx2);
end

denom_safe = double(denominator(:)).';
denom_safe(~isfinite(denom_safe) | denom_safe <= 0) = 1;

for ii = 1:numel(overlap)
    k = overlap(ii);
    chunk_start = double(chunk_table.chunk_start_idx(k));
    chunk_end = double(chunk_table.chunk_end_idx(k));

    global_read_start = max(idx1, chunk_start);
    global_read_end = min(idx2, chunk_end);
    chunk_read_start = global_read_start - chunk_start + 1;
    chunk_read_end = global_read_end - chunk_start + 1;
    dest_start = global_read_start - idx1 + 1;
    dest_end = global_read_end - idx1 + 1;

    S = load(files(k).fullpath, 'EDMD_outputs');
    if ~isfield(S, 'EDMD_outputs') || ~isfield(S.EDMD_outputs, 'efuns')
        error('File %s does not contain EDMD_outputs.efuns.', files(k).fullpath);
    end

    raw_slice = S.EDMD_outputs.efuns( ...
        chunk_read_start:chunk_read_end, selected_mode_idx);
    switch lower(char(feature_variant))
        case 'real'
            feature_slice = real(raw_slice);
        case 'abs'
            feature_slice = abs(raw_slice);
        otherwise
            error('Unsupported feature variant: %s.', char(feature_variant));
    end

    X_feature(dest_start:dest_end, :) = bsxfun(@rdivide, ...
        double(feature_slice), denom_safe);
end
end


function result_window = local_build_window_result(result, X_feature, idx1, idx2)
result_window = struct();
result_window.meta = result.meta;
result_window.cfg = result.cfg;
result_window.input = result.input;
result_window.input.time_axis = local_slice_time_axis(result, idx1, idx2);

result_window.data = struct();
if isfield(result, 'data')
    data_fields = {'evalues_discrete', 'evalues_bilinear', ...
        'omitted_compact_fields'};
    for i = 1:numel(data_fields)
        field_name = data_fields{i};
        if isfield(result.data, field_name)
            result_window.data.(field_name) = result.data.(field_name);
        end
    end
end
result_window.data.efun_feature_time_by_mode = X_feature;

if isfield(result, 'feature')
    result_window.feature = result.feature;
else
    result_window.feature = struct('variant', result.meta.feature_variant);
end

result_window.core = struct();
result_window.core.temporal_components_time_by_comp = ...
    result.core.temporal_components_time_by_comp(idx1:idx2, :);

core_fields = {'mode_weights_mode_by_comp', 'singular_values', ...
    'explained_variance', 'explained_variance_ratio'};
for i = 1:numel(core_fields)
    field_name = core_fields{i};
    if isfield(result.core, field_name)
        result_window.core.(field_name) = result.core.(field_name);
    end
end

result_window.summary = struct();
if isfield(result, 'summary') && ...
        isfield(result.summary, 'temporal_components_smooth_time_by_comp') && ...
        ~isempty(result.summary.temporal_components_smooth_time_by_comp)
    result_window.summary.temporal_components_smooth_time_by_comp = ...
        result.summary.temporal_components_smooth_time_by_comp(idx1:idx2, :);
else
    result_window.summary.temporal_components_smooth_time_by_comp = [];
end
end


function time_axis = local_slice_time_axis(result, idx1, idx2)
if isfield(result, 'input') && isfield(result.input, 'time_axis') && ...
        numel(result.input.time_axis) >= idx2
    time_axis = result.input.time_axis(idx1:idx2);
else
    time_axis = (idx1:idx2).';
end
time_axis = time_axis(:);
end


function file_stub = local_build_file_stub(top_windows, row_idx)
rank_value = local_table_numeric_value(top_windows, row_idx, ...
    'state_diversity_rank', row_idx);
global_window_idx = local_table_numeric_value(top_windows, row_idx, ...
    'global_window_idx', row_idx);

file_stub = sprintf('r%02d_w%03d', ...
    round(rank_value), round(global_window_idx));
end


function title_text = local_build_window_title(top_windows, row_idx, idx1, idx2)
rank_value = local_table_numeric_value(top_windows, row_idx, ...
    'state_diversity_rank', row_idx);
global_window_idx = local_table_numeric_value(top_windows, row_idx, ...
    'global_window_idx', row_idx);
richness = local_table_numeric_value(top_windows, row_idx, ...
    'active_state_richness', NaN);
total_state_count = local_table_numeric_value(top_windows, row_idx, ...
    'total_state_window_count', NaN);
h_norm = local_table_numeric_value(top_windows, row_idx, ...
    'normalized_state_entropy', NaN);
dominant_state = local_table_text_value(top_windows, row_idx, ...
    'dominant_state', 'NA');

title_text = sprintf(['Rank %d | globalwin %d | samples [%d,%d] | ', ...
    'richness=%g | total states=%g | dominant=%s | Hnorm=%.3f'], ...
    round(rank_value), round(global_window_idx), idx1, idx2, ...
    richness, total_state_count, dominant_state, h_norm);
end


function value = local_table_numeric_value(T, row_idx, name, default_value)
if ~ismember(name, T.Properties.VariableNames)
    value = default_value;
    return;
end

x = T.(name)(row_idx);
if iscell(x)
    x = x{1};
end
value = double(x);
if isempty(value) || ~isfinite(value)
    value = default_value;
end
end


function value = local_table_text_value(T, row_idx, name, default_value)
if ~ismember(name, T.Properties.VariableNames)
    value = default_value;
    return;
end

x = T.(name)(row_idx);
if iscell(x)
    x = x{1};
end
value = char(string(x));
if isempty(value)
    value = default_value;
end
end


function event_windows = local_build_consensus_event_windows(C, x_axis, idx1, idx2, params)
event_windows = zeros(0, 3);

if ~isfield(C, 'state_code_by_time') || isempty(C.state_code_by_time)
    return;
end

state_idx1 = double(params.state_start_idx) + idx1 - 1;
state_idx2 = double(params.state_start_idx) + idx2 - 1;
if state_idx1 < 1 || state_idx2 > numel(C.state_code_by_time)
    warning('Consensus state vector does not cover sample range [%d, %d].', ...
        idx1, idx2);
    return;
end

state_codes = double(C.state_code_by_time(state_idx1:state_idx2));
state_codes = state_codes(:);
x_axis = x_axis(:);

if numel(x_axis) ~= numel(state_codes)
    error('x_axis and consensus state window lengths do not match.');
end

run_bounds = local_find_constant_runs(state_codes);
for r = 1:size(run_bounds, 1)
    s = run_bounds(r, 1);
    e = run_bounds(r, 2);
    code = state_codes(s);
    if code == double(params.baseline_code)
        continue;
    end
    event_windows(end + 1, :) = [x_axis(s), x_axis(e), code]; %#ok<AGROW>
end
end


function run_bounds = local_find_constant_runs(x)
x = x(:);
if isempty(x)
    run_bounds = zeros(0, 2);
    return;
end

starts = [1; find(diff(x) ~= 0) + 1];
ends = [starts(2:end) - 1; numel(x)];
run_bounds = [starts, ends];
end


function colors = local_build_consensus_event_colors(C, cfg)
if isfield(cfg.viz, 'state_space_consensus') && ...
        isfield(cfg.viz.state_space_consensus, 'state_colormap') && ...
        ~isempty(cfg.viz.state_space_consensus.state_colormap)
    state_cmap = cfg.viz.state_space_consensus.state_colormap;
else
    state_cmap = [ ...
        0.34, 0.36, 0.38; ...
        0.86, 0.74, 0.42; ...
        0.45, 0.64, 0.90; ...
        0.86, 0.55, 0.75; ...
        0.63, 0.53, 0.88; ...
        0.78, 0.45, 0.60];
end

max_code = 5;
if isfield(C, 'state_catalog') && ~isempty(C.state_catalog) && ...
        isfield(C.state_catalog, 'code')
    max_code = max(max_code, max(double([C.state_catalog.code])));
end

needed_rows = max_code + 1;
if size(state_cmap, 1) < needed_rows
    state_cmap = [state_cmap; lines(needed_rows - size(state_cmap, 1))];
end

colors = state_cmap(2:needed_rows, :);
end


function local_close_if_requested(fig, params)
if params.close_figures && ~isempty(fig) && isvalid(fig)
    close(fig);
end
end


function save_path = local_export_existing_plot(fig, save_dir, file_prefix, ...
    result, save_tag, resolution)
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

local_disable_axes_toolbars(fig);
save_path = local_build_eigenfunction_png_path( ...
    save_dir, file_prefix, result, save_tag);
exportgraphics(fig, save_path, 'Resolution', resolution);
end


function local_disable_axes_toolbars(fig)
axes_list = findall(fig, 'Type', 'axes');
for i = 1:numel(axes_list)
    ax = axes_list(i);
    try
        ax.Toolbar.Visible = 'off';
    catch
    end

    try
        axtoolbar(ax, {});
    catch
    end
end
drawnow;
end


function save_path = local_build_eigenfunction_png_path( ...
    save_dir, file_prefix, result, save_tag)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {file_prefix, lower(result.meta.path_kind), ...
    lower(result.meta.feature_variant)};

if nargin >= 4 && ~isempty(save_tag)
    pieces{end + 1} = save_tag;
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(save_dir, sprintf('%s__%s.png', filename, timestamp));
end
