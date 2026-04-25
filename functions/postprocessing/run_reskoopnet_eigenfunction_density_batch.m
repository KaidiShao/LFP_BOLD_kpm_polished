function manifest = run_reskoopnet_eigenfunction_density_batch(params)
%RUN_RESKOOPNET_EIGENFUNCTION_DENSITY_BATCH
% Compute thresholded raw-eigenfunction density for completed ResKoopNet runs.
%
% This is intentionally density-only: it does not run eigenfunction
% dimension reduction and does not compute density on reduced components.

if nargin < 1 || isempty(params)
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
    error('No completed ResKoopNet output runs were found under %s.', ...
        params.autodl_root);
end

created_tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest_csv_file = fullfile(params.combined_manifest_root, ...
    sprintf('reskoopnet_efun_density_scan_manifest_%s.csv', created_tag));
manifest_mat_file = fullfile(params.combined_manifest_root, ...
    sprintf('reskoopnet_efun_density_scan_manifest_%s.mat', created_tag));

fprintf('Discovered %d completed ResKoopNet runs.\n', numel(runs));
fprintf('Threshold ratios: %s\n', mat2str(params.threshold_ratios));
fprintf('Combined manifest:\n  %s\n\n', manifest_csv_file);

max_rows = max(1, numel(runs) * numel(params.threshold_ratios));
rows = repmat(local_empty_manifest_row(), max_rows, 1);
row_count = 0;

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n============================================================\n');
    fprintf('[RUN %d/%d] %s | %s\n', ...
        i_run, numel(runs), run_info.dataset, run_info.run_name);
    fprintf('Source EDMD chunks:\n  %s\n', run_info.output_dir);

    run_output_root = local_run_output_root(params, run_info);
    mat_dir = fullfile(run_output_root, 'mat');
    fig_dir = fullfile(run_output_root, 'fig');
    observable_file_guess = local_observable_file_for_run(params, run_info);

    ratios = params.threshold_ratios(:).';
    needs_run = true(size(ratios));

    if params.skip_existing
        for i_thr = 1:numel(ratios)
            [existing_mat, existing_fig] = local_find_existing_density_result( ...
                mat_dir, fig_dir, run_info, params, ratios(i_thr));
            if ~isempty(existing_mat)
                row_count = row_count + 1;
                rows = local_ensure_row_capacity(rows, row_count);
                rows(row_count) = local_manifest_row( ...
                    run_info, params, ratios(i_thr), 'skipped_existing', ...
                    'Existing density MAT found.', existing_mat, existing_fig, ...
                    '', run_output_root, observable_file_guess, NaN, '', ...
                    run_info.chunk_count, NaN, NaN, NaN, NaN, 0);
                needs_run(i_thr) = false;
            end
        end

        if ~any(needs_run)
            fprintf('All requested thresholds already exist; skipping load.\n');
            local_write_manifest(rows(1:row_count), manifest_csv_file);
            continue;
        end
    end

    t_run = tic;
    try
        source_cfg = local_source_cfg_for_run(run_info, params);
        [EDMD_outputs, concat_info] = load_edmd_source(source_cfg);

        observable_file = local_observable_file_for_run( ...
            params, run_info, EDMD_outputs);
        prep = local_prepare_eigenfunction_density_input( ...
            EDMD_outputs, params, observable_file);

        fprintf('Loaded T=%d samples, selected modes=%d, dt=%.10g sec.\n', ...
            size(prep.phi, 1), size(prep.phi, 2), prep.dt);
        if ~isempty(observable_file)
            fprintf('Observable metadata:\n  %s\n', observable_file);
        end

        for i_thr = 1:numel(ratios)
            ratio = ratios(i_thr);
            if ~needs_run(i_thr)
                continue;
            end

            fprintf('  threshold %g (%d/%d)...\n', ...
                ratio, i_thr, numel(ratios));

            density_params = local_density_params( ...
                params, run_info, prep, ratio, mat_dir, fig_dir, ...
                observable_file);

            t_thr = tic;
            [D, fig] = get_thresholded_density(prep.phi, prep.dt, density_params);
            if params.close_figures && ~isempty(fig) && isvalid(fig)
                close(fig);
            end
            runtime_sec = toc(t_thr);

            row_count = row_count + 1;
            rows = local_ensure_row_capacity(rows, row_count);
            rows(row_count) = local_manifest_row( ...
                run_info, params, ratio, 'ok', '', ...
                D.artifacts.mat_file, D.artifacts.figure_file, ...
                D.artifacts.file_stem, run_output_root, observable_file, ...
                prep.dt, prep.dt_source_text, concat_info.n_chunks, ...
                size(prep.phi, 1), size(prep.phi, 2), ...
                numel(D.t_centers), runtime_sec, toc(t_run));

            fprintf('    saved MAT:\n      %s\n', D.artifacts.mat_file);
            local_write_manifest(rows(1:row_count), manifest_csv_file);
        end

        clear EDMD_outputs prep D fig
    catch ME
        message = getReport(ME, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s failed:\n%s\n', run_info.run_name, message);

        for i_thr = 1:numel(ratios)
            if ~needs_run(i_thr)
                continue;
            end
            row_count = row_count + 1;
            rows = local_ensure_row_capacity(rows, row_count);
            rows(row_count) = local_manifest_row( ...
                run_info, params, ratios(i_thr), 'error', message, ...
                '', '', '', run_output_root, observable_file_guess, ...
                NaN, '', run_info.chunk_count, NaN, NaN, NaN, NaN, ...
                toc(t_run));
        end

        local_write_manifest(rows(1:row_count), manifest_csv_file);
        close all force;
        if ~params.continue_on_error
            rethrow(ME);
        end
    end
end

if row_count == 0
    out_rows = local_empty_manifest_row();
    out_rows = out_rows([]);
else
    out_rows = rows(1:row_count);
end

manifest = struct();
manifest.params = params;
manifest.runs = runs;
manifest.rows = out_rows;
manifest.table = local_rows_to_table(out_rows);
manifest.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
manifest.csv_file = manifest_csv_file;
manifest.mat_file = manifest_mat_file;

if ~isempty(out_rows)
    writetable(manifest.table, manifest.csv_file);
end
save(manifest.mat_file, 'manifest', '-v7.3');

fprintf('\nEigenfunction density batch finished.\n');
fprintf('Rows: %d\n', height(manifest.table));
fprintf('CSV manifest:\n  %s\n', manifest.csv_file);
fprintf('MAT manifest:\n  %s\n', manifest.mat_file);
end


function params = local_apply_defaults(params)
if ~isfield(params, 'autodl_root') || isempty(params.autodl_root)
    params.autodl_root = 'E:\autodl_results';
end
if ~isfield(params, 'processed_root') || isempty(params.processed_root)
    params.processed_root = get_project_processed_root();
end
if ~isfield(params, 'dataset_stems')
    params.dataset_stems = {};
end
params.dataset_stems = local_normalize_cellstr(params.dataset_stems);

if ~isfield(params, 'exclude_dataset_stems')
    params.exclude_dataset_stems = {};
end
params.exclude_dataset_stems = lower(local_normalize_cellstr( ...
    params.exclude_dataset_stems));

if ~isfield(params, 'output_folder_name') || isempty(params.output_folder_name)
    params.output_folder_name = 'efun_den';
end
if ~isfield(params, 'legacy_output_folder_names')
    params.legacy_output_folder_names = {'eigenfunction_density_scan'};
end
params.legacy_output_folder_names = local_normalize_cellstr( ...
    params.legacy_output_folder_names);
if ~isfield(params, 'combined_manifest_root') || ...
        isempty(params.combined_manifest_root)
    params.combined_manifest_root = fullfile(params.processed_root, ...
        'postprocessing_manifests', params.output_folder_name);
end

if ~isfield(params, 'filename_pattern') || isempty(params.filename_pattern)
    params.filename_pattern = '*_outputs_*.mat';
end
if ~isfield(params, 'variable_name') || isempty(params.variable_name)
    params.variable_name = 'EDMD_outputs';
end
if ~isfield(params, 'min_chunk_files') || isempty(params.min_chunk_files)
    params.min_chunk_files = 1;
end
if ~isfield(params, 'deduplicate_conditions') || ...
        isempty(params.deduplicate_conditions)
    params.deduplicate_conditions = false;
end

if ~isfield(params, 'threshold_ratios') || isempty(params.threshold_ratios)
    params.threshold_ratios = [0.1:0.1:0.9, 0.99];
end
params.threshold_ratios = double(params.threshold_ratios(:)).';
if any(~isfinite(params.threshold_ratios))
    error('params.threshold_ratios must contain finite numeric values.');
end

if ~isfield(params, 'threshold_mode') || isempty(params.threshold_mode)
    params.threshold_mode = 'quantile';
end
if ~isfield(params, 'window_sec') || isempty(params.window_sec)
    params.window_sec = 2;
end
if ~isfield(params, 'step_sec') || isempty(params.step_sec)
    params.step_sec = params.window_sec;
end
if ~isfield(params, 'value_transform') || isempty(params.value_transform)
    params.value_transform = 'none';
end
if ~isfield(params, 'feature_variant') || isempty(params.feature_variant)
    params.feature_variant = 'abs';
end
if ~isfield(params, 'density_denominator') || ...
        isempty(params.density_denominator)
    params.density_denominator = 'window_samples';
end
if ~isfield(params, 'smooth_density') || isempty(params.smooth_density)
    params.smooth_density = false;
end
if ~isfield(params, 'smooth_window_sec') || isempty(params.smooth_window_sec)
    params.smooth_window_sec = 0;
end
if ~isfield(params, 'output_class') || isempty(params.output_class)
    params.output_class = 'single';
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
if ~isfield(params, 'max_modes') || isempty(params.max_modes)
    params.max_modes = Inf;
end
if ~isfield(params, 'dt') || isempty(params.dt)
    params.dt = [];
end
if ~isfield(params, 'require_session_metadata') || ...
        isempty(params.require_session_metadata)
    params.require_session_metadata = true;
end

if ~isfield(params, 'save_results') || isempty(params.save_results)
    params.save_results = true;
end
if ~isfield(params, 'make_figure') || isempty(params.make_figure)
    params.make_figure = true;
end
if ~isfield(params, 'save_figure') || isempty(params.save_figure)
    params.save_figure = true;
end
if ~isfield(params, 'save_v7_3') || isempty(params.save_v7_3)
    params.save_v7_3 = true;
end
if ~isfield(params, 'figure_resolution') || isempty(params.figure_resolution)
    params.figure_resolution = 200;
end
if ~isfield(params, 'max_plot_modes') || isempty(params.max_plot_modes)
    params.max_plot_modes = 100;
end

if ~isfield(params, 'skip_existing') || isempty(params.skip_existing)
    params.skip_existing = true;
end
if ~isfield(params, 'continue_on_error') || isempty(params.continue_on_error)
    params.continue_on_error = true;
end
if ~isfield(params, 'headless') || isempty(params.headless)
    params.headless = true;
end
if ~isfield(params, 'close_figures') || isempty(params.close_figures)
    params.close_figures = true;
end
if ~isfield(params, 'verbose') || isempty(params.verbose)
    params.verbose = true;
end
if ~isfield(params, 'progress_every') || isempty(params.progress_every)
    params.progress_every = 10;
end
if ~isfield(params, 'progress_timing') || isempty(params.progress_timing)
    params.progress_timing = true;
end
end


function runs = local_discover_completed_runs(params)
runs = struct('dataset', {}, 'observable_mode', {}, 'residual_form', {}, ...
    'run_name', {}, 'output_dir', {}, 'summary_count', {}, ...
    'chunk_count', {}, 'last_write_time', {});

datasets = params.dataset_stems;
if isempty(datasets)
    L = dir(params.autodl_root);
    L = L([L.isdir]);
    L = L(~ismember({L.name}, {'.', '..'}));
    datasets = {L.name};
end

for i = 1:numel(datasets)
    dataset = char(datasets{i});
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
        if isempty(summary_files) || numel(chunk_files) < params.min_chunk_files
            continue;
        end

        [observable_mode, residual_form] = ...
            local_infer_condition_from_run_name(run_name);
        file_times = [summary_files.datenum, chunk_files.datenum];

        run = struct();
        run.dataset = dataset;
        run.observable_mode = observable_mode;
        run.residual_form = residual_form;
        run.run_name = run_name;
        run.output_dir = output_dir;
        run.summary_count = numel(summary_files);
        run.chunk_count = numel(chunk_files);
        run.last_write_time = max(file_times);
        runs(end+1) = run; %#ok<AGROW>
    end
end

if isempty(runs)
    return;
end

if params.deduplicate_conditions
    runs = local_deduplicate_runs_by_condition(runs);
end

[~, order] = sort(string({runs.dataset}) + "|" + ...
    string({runs.observable_mode}) + "|" + ...
    string({runs.residual_form}) + "|" + string({runs.run_name}));
runs = runs(order);
end


function runs = local_deduplicate_runs_by_condition(runs)
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


function source_cfg = local_source_cfg_for_run(run_info, params)
source_cfg = struct();
source_cfg.mode = 'chunk_dir';
source_cfg.data_dir = run_info.output_dir;
source_cfg.concat = struct();
source_cfg.concat.filename_pattern = params.filename_pattern;
source_cfg.concat.variable_name = params.variable_name;
source_cfg.concat.concat_fields = {'efuns'};
source_cfg.concat.concat_dim = 1;
source_cfg.concat.required_equal_fields = { ...
    'evalues', 'N_dict', 'residual_form', ...
    'observable_tag', 'observable_mode', ...
    'dt', 'dx', 'sampling_period', 'sample_period', ...
    'fs', 'sampling_frequency', ...
    'session_dx', 'session_fs', 'session_ids', 'session_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
source_cfg.concat.allow_missing_chunks = false;
source_cfg.concat.scan_mode = 'uniform_except_last';
source_cfg.concat.verbose = params.verbose;
source_cfg.concat.progress_every = params.progress_every;
source_cfg.concat.progress_timing = params.progress_timing;
end


function prep = local_prepare_eigenfunction_density_input( ...
        EDMD_outputs, params, observable_file)
required_fields = {'evalues', 'efuns'};
for i = 1:numel(required_fields)
    if ~isfield(EDMD_outputs, required_fields{i})
        error('EDMD_outputs.%s is required.', required_fields{i});
    end
end

evalues0 = EDMD_outputs.evalues(:);
efuns0 = EDMD_outputs.efuns;
if size(efuns0, 2) ~= numel(evalues0)
    error('size(EDMD_outputs.efuns, 2) must equal numel(EDMD_outputs.evalues).');
end

ord = local_sort_eigenvalues(evalues0, params.sort_by, params.sort_dir);
evalues_sorted = evalues0(ord);
mask_sorted = abs(evalues_sorted) > params.abs_thresh;
idx_sorted_selected = ord(mask_sorted);

if isempty(idx_sorted_selected)
    error('No modes remain after applying abs_thresh=%g.', params.abs_thresh);
end

if isfinite(params.max_modes)
    max_modes = min(numel(idx_sorted_selected), params.max_modes);
    idx_sorted_selected = idx_sorted_selected(1:max_modes);
end

efun_raw = efuns0(:, idx_sorted_selected);
phi = normalize_efun(efun_raw, params.feature_variant);
evalues_selected = evalues0(idx_sorted_selected);

[dt, dt_source_text] = local_resolve_dt(params, EDMD_outputs, observable_file);
if isempty(dt)
    error('Could not resolve dt for density windows. %s', dt_source_text);
end

prep = struct();
prep.phi = phi;
prep.dt = dt;
prep.dt_source_text = dt_source_text;
prep.time_axis = (0:size(phi, 1)-1).' * dt;
prep.evalues_discrete = evalues_selected(:);
prep.mode_index = (1:numel(evalues_selected)).';
prep.selected_mode_idx_in_original = idx_sorted_selected(:);
prep.session = local_session_params_from_edmd(EDMD_outputs);
end


function ord = local_sort_eigenvalues(evalues, sort_by, sort_dir)
switch lower(sort_by)
    case {'modulus', 'abs'}
        key = abs(evalues);
    case {'real', 'realpart'}
        key = real(evalues);
    otherwise
        error('Unknown sort_by = %s.', sort_by);
end
[~, ord] = sort(key, sort_dir);
end


function session_params = local_session_params_from_edmd(EDMD_outputs)
session_params = struct();
names = {'session_start_idx', 'session_end_idx', 'session_lengths', ...
    'session_ids', 'session_dx'};
for i = 1:numel(names)
    name = names{i};
    if isfield(EDMD_outputs, name)
        session_params.(name) = EDMD_outputs.(name);
    end
end
end


function density_params = local_density_params( ...
        params, run_info, prep, ratio, mat_dir, fig_dir, observable_file)
density_params = struct();
density_params.window_sec = params.window_sec;
density_params.step_sec = params.step_sec;
density_params.threshold_ratio = ratio;
density_params.threshold_mode = params.threshold_mode;
density_params.value_transform = params.value_transform;
density_params.density_denominator = params.density_denominator;
density_params.smooth_density = params.smooth_density;
density_params.smooth_window_sec = params.smooth_window_sec;
density_params.output_class = params.output_class;
density_params.time_axis = prep.time_axis;
density_params.mode_index = prep.mode_index;
density_params.selected_mode_idx_in_original = ...
    prep.selected_mode_idx_in_original;
density_params.observable_file = observable_file;
density_params.require_session_metadata = params.require_session_metadata;
density_params.save_results = params.save_results;
density_params.save_dir = mat_dir;
density_params.save_stem = local_density_save_stem(run_info);
density_params.save_tag = local_run_save_tag(run_info);
density_params.save_v7_3 = params.save_v7_3;
density_params.make_figure = params.make_figure;
density_params.save_figure = params.save_figure;
density_params.figure_dir = fig_dir;
density_params.figure_visible = 'off';
density_params.figure_resolution = params.figure_resolution;
density_params.max_plot_modes = params.max_plot_modes;
density_params.title = sprintf('%s | %s | Eigenfunction Density', ...
    run_info.dataset, run_info.run_name);

session_fields = fieldnames(prep.session);
for i = 1:numel(session_fields)
    density_params.(session_fields{i}) = prep.session.(session_fields{i});
end
end


function [dt, source_text] = local_resolve_dt(params, EDMD_outputs, observable_file)
dt = []; %#ok<NASGU>
source_text = ''; %#ok<NASGU>

if ~isempty(params.dt)
    dt = local_positive_scalar(params.dt);
    if isempty(dt)
        error('params.dt must be empty or a positive finite scalar.');
    end
    source_text = 'params.dt';
    return;
end

[dt, field_name] = local_read_first_positive_scalar( ...
    EDMD_outputs, {'dt', 'dx', 'sampling_period', 'sample_period'});
if ~isempty(dt)
    source_text = sprintf('EDMD_outputs.%s', field_name);
    return;
end

if isempty(observable_file)
    source_text = ['No dt/dx was found in EDMD_outputs, and no ', ...
        'observable dictionary file could be inferred.'];
    return;
end
if exist(observable_file, 'file') ~= 2
    source_text = sprintf(['No dt/dx was found in EDMD_outputs. ', ...
        'Observable dictionary file does not exist: %s'], observable_file);
    return;
end

[dt, field_name, detail_msg] = local_read_dt_from_observable_file(observable_file);
if ~isempty(dt)
    source_text = sprintf('observable_file.%s', field_name);
else
    source_text = sprintf(['No dt/dx was found in EDMD_outputs, and ', ...
        'observable dictionary %s did not contain usable time metadata. %s'], ...
        observable_file, detail_msg);
end
end


function [dt, field_name, detail_msg] = local_read_dt_from_observable_file(observable_file)
dt = [];
field_name = '';
detail_msg = '';

try
    vars = who('-file', observable_file);
catch ME
    detail_msg = sprintf('Could not inspect observable file: %s', ME.message);
    return;
end

fields_to_load = intersect( ...
    {'dx', 'dt', 'sampling_period', 'sample_period', ...
    'fs', 'sampling_frequency', 'session_dx', 'session_fs'}, ...
    vars, 'stable');
if isempty(fields_to_load)
    detail_msg = 'No candidate time-metadata variables were present.';
    return;
end

try
    S = load(observable_file, fields_to_load{:});
catch ME
    detail_msg = sprintf('Could not load observable time metadata: %s', ME.message);
    return;
end

[dt, field_name] = local_read_first_positive_scalar( ...
    S, {'dx', 'dt', 'sampling_period', 'sample_period'});
if ~isempty(dt)
    return;
end

if isfield(S, 'session_dx')
    dt = local_uniform_positive_scalar(S.session_dx);
    if ~isempty(dt)
        field_name = 'session_dx';
        return;
    end
end

[fs, fs_field] = local_read_first_positive_scalar( ...
    S, {'fs', 'sampling_frequency'});
if isempty(fs) && isfield(S, 'session_fs')
    fs = local_uniform_positive_scalar(S.session_fs);
    if ~isempty(fs)
        fs_field = 'session_fs';
    end
end
if ~isempty(fs)
    dt = 1 / fs;
    field_name = fs_field;
    return;
end

detail_msg = 'Candidate variables were present but did not resolve to a positive scalar dt.';
end


function [value, field_name] = local_read_first_positive_scalar(S, field_names)
value = [];
field_name = '';
for i = 1:numel(field_names)
    name = field_names{i};
    if ~isfield(S, name)
        continue;
    end
    value = local_positive_scalar(S.(name));
    if ~isempty(value)
        field_name = name;
        return;
    end
end
end


function value = local_positive_scalar(value_in)
value = [];
if ~isnumeric(value_in) || isempty(value_in)
    return;
end
if ~isreal(value_in)
    value_in = real(value_in);
end
value_in = double(value_in);
if ~isscalar(value_in) || ~isfinite(value_in) || value_in <= 0
    return;
end
value = value_in;
end


function value = local_uniform_positive_scalar(value_in)
value = [];
if ~isnumeric(value_in) || isempty(value_in)
    return;
end
if ~isreal(value_in)
    value_in = real(value_in);
end
values = double(value_in(:));
values = values(isfinite(values) & values > 0);
if isempty(values)
    return;
end
value0 = median(values);
tol = max(1e-12, abs(value0) * 1e-5);
if all(abs(values - value0) <= tol)
    value = value0;
end
end


function observable_file = local_observable_file_for_run(params, run_info, EDMD_outputs)
if nargin >= 3 && isfield(EDMD_outputs, 'observable_file') && ...
        ~isempty(EDMD_outputs.observable_file)
    observable_file = char(string(EDMD_outputs.observable_file));
    if exist(observable_file, 'file') == 2
        return;
    end
end
if nargin >= 3 && isfield(EDMD_outputs, 'data_full_path') && ...
        ~isempty(EDMD_outputs.data_full_path)
    observable_file = char(string(EDMD_outputs.data_full_path));
    if exist(observable_file, 'file') == 2
        return;
    end
end

tag = run_info.observable_mode;
if strcmpi(tag, 'unknown')
    tag = 'abs';
end
observable_file = fullfile(params.processed_root, run_info.dataset, ...
    'reskoopnet_dictionary', sprintf('%s_low50_high250_g2_%s_single.mat', ...
    run_info.dataset, tag));
if exist(observable_file, 'file') ~= 2
    observable_file = '';
end
end


function output_root = local_run_output_root(params, run_info)
output_root = fullfile(params.processed_root, run_info.dataset, ...
    'postprocessing', params.output_folder_name, ...
    local_short_run_label(run_info));
end


function save_stem = local_density_save_stem(run_info)
save_stem = sprintf('%s_efd', local_filename_safe(run_info.dataset));
end


function save_tag = local_run_save_tag(run_info)
save_tag = local_short_run_label(run_info);
end


function [mat_file, fig_file] = local_find_existing_density_result( ...
        mat_dir, fig_dir, run_info, params, ratio)
[mat_file, fig_file] = local_find_existing_density_in_dir( ...
    mat_dir, fig_dir, local_density_file_base(run_info, params, ratio));
if ~isempty(mat_file)
    return;
end

legacy_roots = local_legacy_output_roots(params, run_info);
legacy_base = local_legacy_density_file_base(run_info, params, ratio);
for i = 1:numel(legacy_roots)
    legacy_mat_dir = fullfile(legacy_roots{i}, 'mat');
    legacy_fig_dir = fullfile(legacy_roots{i}, 'fig');
    [mat_file, fig_file] = local_find_existing_density_in_dir( ...
        legacy_mat_dir, legacy_fig_dir, legacy_base);
    if ~isempty(mat_file)
        return;
    end
end
end


function [mat_file, fig_file] = local_find_existing_density_in_dir( ...
        mat_dir, fig_dir, base)
mat_file = '';
fig_file = '';
L = dir(fullfile(mat_dir, [base, '__*.mat']));
if isempty(L)
    return;
end

[~, order] = sort([L.datenum], 'descend');
L = L(order);
for i = 1:numel(L)
    candidate = fullfile(L(i).folder, L(i).name);
    if ~local_is_valid_existing_mat(candidate)
        continue;
    end

    mat_file = candidate;
    [~, stem] = fileparts(mat_file);
    candidate_fig = fullfile(fig_dir, [stem, '.png']);
    if exist(candidate_fig, 'file') == 2
        fig_file = candidate_fig;
    end
    return;
end
end


function tf = local_is_valid_existing_mat(mat_file)
tf = false;
if exist(mat_file, 'file') ~= 2
    return;
end
info = dir(mat_file);
if isempty(info) || info.bytes <= 0
    return;
end
tf = true;
end


function roots = local_legacy_output_roots(params, run_info)
folder_names = params.legacy_output_folder_names;
if isempty(folder_names)
    roots = {};
    return;
end

roots = cell(numel(folder_names), 1);
n = 0;
for i = 1:numel(folder_names)
    folder_name = char(folder_names{i});
    if isempty(folder_name)
        continue;
    end
    n = n + 1;
    roots{n} = fullfile(params.processed_root, run_info.dataset, ...
        'postprocessing', folder_name, local_filename_safe(run_info.run_name));
end
roots = roots(1:n);
end


function base = local_density_file_base(run_info, params, ratio)
pieces = {local_density_save_stem(run_info), char(params.threshold_mode), ...
    local_ratio_tag(ratio), local_run_save_tag(run_info)};
base = strjoin(pieces, '__');
base = regexprep(base, '[^\w\-]+', '_');
end


function base = local_legacy_density_file_base(run_info, params, ratio)
pieces = {sprintf('%s_efun_density', local_filename_safe(run_info.dataset)), ...
    char(params.threshold_mode), local_ratio_tag(ratio), ...
    local_filename_safe(run_info.run_name)};
base = strjoin(pieces, '__');
base = regexprep(base, '[^\w\-]+', '_');
end


function tag = local_ratio_tag(ratio)
tag = sprintf('ratio_%03d', round(double(ratio) * 100));
end


function row = local_manifest_row(run_info, params, ratio, status, message, ...
        mat_file, fig_file, file_stem, output_root, observable_file, ...
        dt, dt_source, n_chunks, n_samples, n_modes, n_windows, ...
        runtime_sec, run_elapsed_sec)
row = struct();
row.created_at = string(char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
row.dataset = string(run_info.dataset);
row.run_name = string(run_info.run_name);
row.observable_mode = string(run_info.observable_mode);
row.residual_form = string(run_info.residual_form);
row.source_output_dir = string(run_info.output_dir);
row.output_root = string(output_root);
row.observable_file = string(observable_file);
row.threshold_mode = string(params.threshold_mode);
row.threshold_ratio = double(ratio);
row.window_sec = double(params.window_sec);
row.step_sec = double(params.step_sec);
row.feature_variant = string(params.feature_variant);
row.abs_thresh = double(params.abs_thresh);
row.sort_by = string(params.sort_by);
row.sort_dir = string(params.sort_dir);
row.status = string(status);
row.message = string(message);
row.mat_file = string(mat_file);
row.figure_file = string(fig_file);
row.file_stem = string(file_stem);
row.dt = double(dt);
row.dt_source = string(dt_source);
row.n_chunks = double(n_chunks);
row.n_samples = double(n_samples);
row.n_modes = double(n_modes);
row.n_windows = double(n_windows);
row.runtime_sec = double(runtime_sec);
row.run_elapsed_sec = double(run_elapsed_sec);
end


function row = local_empty_manifest_row()
empty_run = struct('dataset', '', 'observable_mode', '', ...
    'residual_form', '', 'run_name', '', 'output_dir', '');
params = struct('threshold_mode', '', 'window_sec', NaN, ...
    'step_sec', NaN, 'feature_variant', '', 'abs_thresh', NaN, ...
    'sort_by', '', 'sort_dir', '');
row = local_manifest_row(empty_run, params, NaN, '', '', '', '', '', ...
    '', '', NaN, '', NaN, NaN, NaN, NaN, NaN, NaN);
end


function rows = local_ensure_row_capacity(rows, row_count)
if row_count <= numel(rows)
    return;
end
rows(end+max(100, numel(rows))) = local_empty_manifest_row();
end


function T = local_rows_to_table(rows)
if isempty(rows)
    T = table();
else
    T = struct2table(rows(:));
end
end


function local_write_manifest(rows, manifest_file)
if isempty(rows)
    return;
end
manifest_dir = fileparts(manifest_file);
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end
writetable(local_rows_to_table(rows), manifest_file);
end


function out = local_normalize_cellstr(value)
if isempty(value)
    out = {};
elseif ischar(value)
    out = cellstr(string(value));
elseif isstring(value)
    out = cellstr(value(:).');
elseif iscell(value)
    out = cellstr(string(value(:)).');
else
    error('Expected a char, string, or cellstr value.');
end
end


function safe = local_filename_safe(value)
safe = char(string(value));
safe = regexprep(safe, '[^\w\-]+', '_');
safe = regexprep(safe, '_+', '_');
safe = regexprep(safe, '^_|_$', '');
if isempty(safe)
    safe = 'unnamed';
end
end


function label = local_short_run_label(run_info)
residual_tag = local_residual_tag(run_info.residual_form);
observable_tag = local_observable_tag(run_info.observable_mode);
hash_tag = local_text_hash(run_info.run_name);
label = sprintf('%s_%s_%s', residual_tag, observable_tag, hash_tag);
label = local_filename_safe(label);
end


function tag = local_residual_tag(residual_form)
switch lower(char(string(residual_form)))
    case 'projected_kv'
        tag = 'kv';
    case 'projected_vlambda'
        tag = 'vl';
    otherwise
        tag = 'res';
end
end


function tag = local_observable_tag(observable_mode)
switch lower(char(string(observable_mode)))
    case 'abs'
        tag = 'abs';
    case 'complex_split'
        tag = 'cs';
    otherwise
        tag = 'obs';
end
end


function hash = local_text_hash(value)
txt = char(string(value));
h = uint32(5381);
bytes = uint8(txt);
for i = 1:numel(bytes)
    h = uint32(mod(double(h) * 33 + double(bytes(i)), 2^32));
end
hex = lower(dec2hex(double(h), 8));
hash = hex(max(1, end-5):end);
end
