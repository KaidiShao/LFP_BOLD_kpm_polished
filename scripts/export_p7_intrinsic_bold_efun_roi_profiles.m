function result = export_p7_intrinsic_bold_efun_roi_profiles(varargin)
%EXPORT_P7_INTRINSIC_BOLD_EFUN_ROI_PROFILES Export all P7 BOLD efun ROI profiles.
%
% This exporter is independent of P8/P10 xcorr and BLP density.  It maps the
% original P7 BOLD Koopman modes into ROI-average vectors for every retained
% sorted BOLD efun mode in the current best-checkpoint BOLD_POST runs.

p = inputParser;
p.addParameter('processed_root', io_project.get_project_processed_root(), @(x) ischar(x) || isstring(x));
p.addParameter('datapons_root', 'E:\DataPons', @(x) ischar(x) || isstring(x));
p.addParameter('autodl_roots', {'E:\autodl_results_local\bold_wsl'});
p.addParameter('dataset_stems', {'e10gb1', 'e10fV1', 'e10gh1', 'e10gw1', ...
    'f12m01', 'k13m17', 'k13m23'});
p.addParameter('exclude_dataset_stems', {});
p.addParameter('observable_modes', {'global_svd100', 'gsvd100_ds', 'HP_svd100', 'roi_mean'});
p.addParameter('residual_forms', {'projected_vlambda'});
p.addParameter('max_modes', Inf, @(x) isnumeric(x) && isscalar(x));
p.addParameter('include_processed_p7_fallback', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('feature_reduce', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('roi_value_mode', 'mean_abs', @(x) ischar(x) || isstring(x));
p.addParameter('output_dir', fullfile('results', 'pipeline_roi_profile_consistency_current'), ...
    @(x) ischar(x) || isstring(x));
p.parse(varargin{:});
opts = p.Results;

opts.processed_root = char(string(opts.processed_root));
opts.datapons_root = char(string(opts.datapons_root));
opts.output_dir = char(string(opts.output_dir));
opts.autodl_roots = cellstr(string(opts.autodl_roots(:)).');
opts.dataset_stems = cellstr(string(opts.dataset_stems(:)).');
opts.exclude_dataset_stems = cellstr(string(opts.exclude_dataset_stems(:)).');
opts.observable_modes = cellstr(string(opts.observable_modes(:)).');
opts.residual_forms = cellstr(string(opts.residual_forms(:)).');
opts.feature_reduce = char(string(opts.feature_reduce));
opts.roi_value_mode = char(string(opts.roi_value_mode));

if exist(opts.output_dir, 'dir') ~= 7
    mkdir(opts.output_dir);
end

profile_file = fullfile(opts.output_dir, 'p7_intrinsic_bold_efun_roi_profiles_long.csv');
missing_file = fullfile(opts.output_dir, 'p7_intrinsic_bold_efun_roi_profiles_missing.csv');

p7 = build_bold_reskoopnet_postprocessing_params();
p7.processed_root = opts.processed_root;
p7.autodl_roots = opts.autodl_roots;
p7.dataset_stems = opts.dataset_stems;
p7.exclude_dataset_stems = opts.exclude_dataset_stems;
p7.observable_modes = opts.observable_modes;
p7.residual_forms = opts.residual_forms;
p7.dedupe_by_condition = true;
p7.run_selection_mode = 'best_val_metric';
p7.include_smoke = false;
p7.require_summary_file = true;

runs = discover_completed_bold_reskoopnet_runs(p7);
if logical(opts.include_processed_p7_fallback)
    runs = local_append_processed_p7_fallback_runs(runs, opts);
end

fid = fopen(profile_file, 'w');
if fid < 0
    error('Could not open profile CSV for writing: %s', profile_file);
end
cleanup_profile = onCleanup(@() fclose(fid));
local_write_csv_row(fid, local_profile_header());

fid_missing = fopen(missing_file, 'w');
if fid_missing < 0
    error('Could not open missing CSV for writing: %s', missing_file);
end
cleanup_missing = onCleanup(@() fclose(fid_missing));
local_write_csv_row(fid_missing, local_missing_header());

n_runs = numel(runs);
n_profile_rows = 0;
n_missing = 0;

fprintf('Discovered %d current P7 BOLD run(s).\n', n_runs);
for i_run = 1:n_runs
    run = runs(i_run);
    fprintf('[%d/%d] %s | %s | %s\n', i_run, n_runs, ...
        run.dataset_stem, run.observable_mode, run.run_name);
    try
        bold_post_file = local_find_bold_post_file(opts.processed_root, run);
        if isempty(bold_post_file)
            n_missing = n_missing + 1;
            local_write_csv_row(fid_missing, local_missing_row(run, ...
                'missing_bold_post', 'No current P7 BOLD_POST MAT was found.'));
            continue;
        end

        [plot_ctx, B] = build_bold_activation_plot_context(bold_post_file, struct( ...
            'datapons_root', opts.datapons_root, ...
            'feature_reduce', opts.feature_reduce));
        [raw_indices, mode_labels] = local_resolve_sorted_raw_indices( ...
            plot_ctx.evalues, opts.max_modes);
        if isempty(raw_indices)
            n_missing = n_missing + 1;
            local_write_csv_row(fid_missing, local_missing_row(run, ...
                'no_modes', 'No sorted BOLD efun modes were available.'));
            continue;
        end

        [roi_values, map_infos] = local_compute_roi_values(plot_ctx, raw_indices, opts);
        roi_labels = cellstr(string(plot_ctx.region_labels(:)));
        n_roi = numel(roi_labels);
        for i_mode = 1:numel(raw_indices)
            raw_idx = raw_indices(i_mode);
            lambda = plot_ctx.evalues(raw_idx);
            map_info = map_infos{i_mode};
            for i_roi = 1:n_roi
                value = roi_values(i_mode, i_roi);
                if ~isfinite(value)
                    continue;
                end
                n_profile_rows = n_profile_rows + 1;
                local_write_csv_row(fid, { ...
                    run.dataset_stem, run.dataset_id, run.observable_mode, ...
                    run.residual_form, run.run_name, run.output_dir, ...
                    run.autodl_root, run.best_val_metric, run.best_outer_epoch, ...
                    run.completed_outer_epochs, run.training_state_path, ...
                    run.selection_metric_source, bold_post_file, ...
                    local_get_field(B, 'observable_file', ''), ...
                    local_get_field(plot_ctx, 'roi_ts_file', ''), ...
                    i_mode, mode_labels{i_mode}, raw_idx, ...
                    real(lambda), imag(lambda), abs(lambda), ...
                    numel(raw_indices), n_roi, i_roi, roi_labels{i_roi}, ...
                    value, opts.roi_value_mode, opts.feature_reduce, ...
                    local_get_field(plot_ctx, 'source_space_kind', ''), ...
                    local_get_map_info_field(map_info, 'coverage_fraction'), ...
                    local_get_map_info_field(map_info, 'n_covered_voxels'), ...
                    local_get_map_info_field(map_info, 'n_uncovered_voxels')});
            end
        end
    catch ME
        n_missing = n_missing + 1;
        local_write_csv_row(fid_missing, local_missing_row(run, ...
            'error', ME.message));
        fprintf('  ERROR: %s\n', ME.message);
    end
end

clear cleanup_profile cleanup_missing

result = struct();
result.profile_file = profile_file;
result.missing_file = missing_file;
result.n_runs = n_runs;
result.n_profile_rows = n_profile_rows;
result.n_missing = n_missing;

fprintf('Profile rows: %d\n', n_profile_rows);
fprintf('Missing/error rows: %d\n', n_missing);
fprintf('Profile CSV:\n  %s\n', profile_file);
fprintf('Missing CSV:\n  %s\n', missing_file);
end


function file = local_find_bold_post_file(processed_root, run)
post_root = fullfile(processed_root, run.dataset_stem, ...
    io_project.get_pipeline_stage_name(7, 'bold_postprocessing'), ...
    run.run_name, 'mat');
file = fullfile(post_root, [run.run_name, '_bold_post.mat']);
if exist(file, 'file') == 2
    return;
end
L = dir(fullfile(post_root, '*_bold_post.mat'));
if isempty(L)
    file = '';
    return;
end
[~, order] = sort([L.datenum], 'descend');
L = L(order);
file = fullfile(L(1).folder, L(1).name);
end


function runs = local_append_processed_p7_fallback_runs(runs, opts)
existing = containers.Map('KeyType', 'char', 'ValueType', 'logical');
for i = 1:numel(runs)
    existing(local_condition_key(runs(i).dataset_stem, runs(i).observable_mode, ...
        runs(i).residual_form)) = true;
end

for i_ds = 1:numel(opts.dataset_stems)
    dataset = opts.dataset_stems{i_ds};
    if any(strcmpi(dataset, opts.exclude_dataset_stems))
        continue;
    end
    post_root = fullfile(opts.processed_root, dataset, ...
        io_project.get_pipeline_stage_name(7, 'bold_postprocessing'));
    if exist(post_root, 'dir') ~= 7
        continue;
    end
    run_dirs = dir(post_root);
    run_dirs = run_dirs([run_dirs.isdir]);
    run_dirs = run_dirs(~ismember({run_dirs.name}, {'.', '..'}));
    for i_obs = 1:numel(opts.observable_modes)
        obs = opts.observable_modes{i_obs};
        for i_form = 1:numel(opts.residual_forms)
            form = opts.residual_forms{i_form};
            key = local_condition_key(dataset, obs, form);
            if isKey(existing, key)
                continue;
            end
            candidates = repmat(local_empty_run(), 0, 1);
            for i_run = 1:numel(run_dirs)
                run_name = run_dirs(i_run).name;
                if ~contains(run_name, obs, 'IgnoreCase', true) || ...
                        ~contains(run_name, form, 'IgnoreCase', true)
                    continue;
                end
                mat_dir = fullfile(run_dirs(i_run).folder, run_name, 'mat');
                L = dir(fullfile(mat_dir, '*_bold_post.mat'));
                if isempty(L)
                    continue;
                end
                one = local_empty_run();
                one.dataset_stem = dataset;
                one.dataset_id = io_project.get_dataset_id_from_stem(dataset);
                one.run_name = run_name;
                one.output_dir = '';
                one.autodl_root = 'processed_p7_fallback';
                one.observable_mode = obs;
                one.residual_form = form;
                one.n_output_files = 0;
                one.n_summary_files = 0;
                one.last_write_time = max([L.datenum]);
                one.selection_metric_source = 'processed_p7_fallback';
                candidates(end + 1, 1) = one; %#ok<AGROW>
            end
            if isempty(candidates)
                continue;
            end
            [~, newest] = max([candidates.last_write_time]);
            runs(end + 1, 1) = candidates(newest); %#ok<AGROW>
            existing(key) = true;
        end
    end
end

if isempty(runs)
    return;
end
[~, order] = sort(string({runs.dataset_stem}) + "|" + ...
    string({runs.observable_mode}) + "|" + string({runs.residual_form}) + "|" + ...
    string({runs.run_name}));
runs = runs(order);
end


function key = local_condition_key(dataset, observable, residual_form)
key = char(lower(strjoin(string({dataset, observable, residual_form}), '|')));
end


function [raw_indices, labels] = local_resolve_sorted_raw_indices(evalues, max_modes)
evalues = evalues(:);
[~, order] = sort(abs(evalues), 'descend');
if isfinite(max_modes)
    n_keep = min(numel(order), max(0, floor(max_modes)));
else
    n_keep = numel(order);
end
raw_indices = order(1:n_keep).';
labels = arrayfun(@(k) sprintf('sorted%03d', k), 1:n_keep, 'UniformOutput', false);
end


function [roi_values, map_infos] = local_compute_roi_values(plot_ctx, raw_indices, opts)
n_mode = numel(raw_indices);
n_roi = numel(plot_ctx.region_labels);
roi_values = nan(n_mode, n_roi);
map_infos = cell(n_mode, 1);
for i_mode = 1:n_mode
    raw_idx = raw_indices(i_mode);
    if raw_idx < 1 || raw_idx > size(plot_ctx.kpm_modes, 1)
        continue;
    end
    mode_obs = double(plot_ctx.kpm_modes(raw_idx, :));
    if ~plot_ctx.direct_voxel_mode
        source_values = mode_obs * plot_ctx.coeff_t;
    else
        source_values = mode_obs;
    end
    [voxel_values, map_info] = map_bold_source_values_to_voxels( ...
        source_values(:), plot_ctx, struct('feature_reduce', opts.feature_reduce));
    map_infos{i_mode} = map_info;
    roi_values(i_mode, :) = local_reduce_voxels_to_regions( ...
        voxel_values(:), plot_ctx.voxel_region_idx(:), n_roi, opts.roi_value_mode);
end
end


function roi_values = local_reduce_voxels_to_regions(voxel_values, voxel_region_idx, n_roi, roi_value_mode)
roi_values = nan(1, n_roi);
for i_roi = 1:n_roi
    vals = voxel_values(voxel_region_idx == i_roi);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    switch lower(char(string(roi_value_mode)))
        case 'mean_abs'
            roi_values(i_roi) = mean(abs(vals), 'omitnan');
        case 'abs_mean'
            roi_values(i_roi) = abs(mean(vals, 'omitnan'));
        case 'real_mean'
            roi_values(i_roi) = mean(real(vals), 'omitnan');
        case 'imag_mean'
            roi_values(i_roi) = mean(imag(vals), 'omitnan');
        otherwise
            error('Unsupported roi_value_mode: %s', roi_value_mode);
    end
end
end


function header = local_profile_header()
header = { ...
    'dataset', 'dataset_id', 'observable', 'residual_form', ...
    'run_name', 'output_dir', 'autodl_root', ...
    'best_val_metric', 'best_outer_epoch', 'completed_outer_epochs', ...
    'training_state_path', 'selection_metric_source', ...
    'bold_post_file', 'observable_file', 'roi_ts_file', ...
    'sorted_position', 'mode_label', 'raw_index', ...
    'eigenvalue_real', 'eigenvalue_imag', 'eigenvalue_abs', ...
    'n_export_modes', 'n_rois', 'roi_index', 'roi_label', ...
    'roi_value', 'roi_value_mode', 'feature_reduce', ...
    'source_space_kind', 'map_coverage_fraction', ...
    'map_n_covered_voxels', 'map_n_uncovered_voxels'};
end


function header = local_missing_header()
header = { ...
    'dataset', 'dataset_id', 'observable', 'residual_form', ...
    'run_name', 'output_dir', 'autodl_root', ...
    'status', 'message'};
end


function row = local_missing_row(run, status, message)
row = {run.dataset_stem, run.dataset_id, run.observable_mode, ...
    run.residual_form, run.run_name, run.output_dir, run.autodl_root, ...
    status, message};
end


function run_info = local_empty_run()
run_info = struct('dataset_stem', '', 'dataset_id', '', 'run_name', '', ...
    'output_dir', '', 'autodl_root', '', 'observable_mode', '', ...
    'residual_form', '', 'n_output_files', 0, 'n_summary_files', 0, ...
    'last_write_time', NaN, 'best_val_metric', NaN, ...
    'best_outer_epoch', NaN, 'completed_outer_epochs', NaN, ...
    'training_state_path', '', 'selection_metric_source', 'metric_missing');
end


function value = local_get_map_info_field(map_info, name)
if isstruct(map_info) && isfield(map_info, name) && ~isempty(map_info.(name))
    value = map_info.(name);
else
    value = NaN;
end
end


function value = local_get_field(S, name, default_value)
if nargin < 3
    default_value = '';
end
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function local_write_csv_row(fid, row)
for i = 1:numel(row)
    if i > 1
        fprintf(fid, ',');
    end
    fprintf(fid, '%s', local_csv_cell(row{i}));
end
fprintf(fid, '\n');
end


function text = local_csv_cell(value)
if isnumeric(value)
    if isempty(value) || ~isscalar(value) || ~isfinite(value)
        text = '';
    else
        text = sprintf('%.12g', value);
    end
elseif isstring(value) || ischar(value)
    text = char(string(value));
elseif islogical(value)
    text = char(string(value));
else
    try
        text = char(string(value));
    catch
        text = '';
    end
end
text = strrep(text, '"', '""');
if contains(text, {',', '"', newline, char(13)})
    text = ['"', text, '"'];
end
end
