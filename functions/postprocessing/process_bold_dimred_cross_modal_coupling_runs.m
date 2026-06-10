function manifest = process_bold_dimred_cross_modal_coupling_runs(params)
%PROCESS_BOLD_DIMRED_CROSS_MODAL_COUPLING_RUNS Batch pipeline 10.

if nargin < 1 || isempty(params)
    params = build_bold_dimred_cross_modal_coupling_params();
else
    params = local_apply_defaults(params);
end

if params.headless
    set(groot, 'defaultFigureVisible', 'off');
end
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
end

candidates = discover_completed_bold_dimred_cross_modal_coupling_runs(params);
if isempty(candidates)
    if params.continue_on_error
        warning('No pipeline 10 candidate P9 result files were found. Skipping pipeline 10.');
        manifest = local_empty_manifest(params, candidates);
        return;
    end
    error('No pipeline 10 candidate P9 result files were found.');
end
if exist(params.manifest_dir, 'dir') ~= 7
    mkdir(params.manifest_dir);
end

rows = repmat(local_empty_row(), numel(candidates), 1);
row_count = 0;
fprintf('Discovered %d pipeline 10 candidate run(s).\n', numel(candidates));

for i_run = 1:numel(candidates)
    candidate = candidates(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s | %s | %s | %s\n', ...
        i_run, numel(candidates), candidate.dataset_stem, ...
        candidate.observable_mode, candidate.residual_form, ...
        candidate.feature_name, candidate.method_tag);
    fprintf('P9 result:\n  %s\n', candidate.dimred_result_file);

    try
        one_result = run_one_bold_dimred_cross_modal_coupling_core(candidate, params);
        fprintf('XCORR status: %s\n', one_result.xcorr_status);
        fprintf('Activation status: %s\n', one_result.activation_status);
        fprintf('ROI summary status: %s\n', one_result.roi_summary_status);
        fprintf('XCORR MAT:\n  %s\n', one_result.xcorr_mat_file);
        row_count = row_count + 1;
        rows(row_count) = local_make_row(one_result);
    catch ME
        message = getReport(ME, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s failed:\n%s\n', candidate.p10_tag, message);
        row_count = row_count + 1;
        rows(row_count) = local_error_row(candidate, message);
        if ~params.continue_on_error
            rethrow(ME);
        end
    end

    if params.close_figures
        close all force;
        drawnow limitrate;
    end
    clear one_result
end

rows = rows(1:row_count);
T = struct2table(rows, 'AsArray', true);
tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest = struct();
manifest.params = params;
manifest.candidates = candidates;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(params.manifest_dir, ...
    ['pipeline10_dimred_xcorr_manifest_', tag, '.csv']);
manifest.mat_file = fullfile(params.manifest_dir, ...
    ['pipeline10_dimred_xcorr_manifest_', tag, '.mat']);

writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');
fprintf('\nPipeline 10 dimred cross-modal coupling finished.\n');
fprintf('Manifest:\n  %s\n', manifest.csv_file);
end


function manifest = local_empty_manifest(params, candidates)
if exist(params.manifest_dir, 'dir') ~= 7
    mkdir(params.manifest_dir);
end
rows = repmat(local_empty_row(), 0, 1);
T = struct2table(rows);
tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest = struct();
manifest.params = params;
manifest.candidates = candidates;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(params.manifest_dir, ...
    ['pipeline10_dimred_xcorr_manifest_empty_', tag, '.csv']);
manifest.mat_file = fullfile(params.manifest_dir, ...
    ['pipeline10_dimred_xcorr_manifest_empty_', tag, '.mat']);
writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');
fprintf('Pipeline 10 skipped because no P9 dimred candidates were found.\n');
fprintf('Empty manifest:\n  %s\n', manifest.csv_file);
end


function params = local_apply_defaults(params)
defaults = build_bold_dimred_cross_modal_coupling_params();
params = local_merge_defaults(defaults, params);
end


function out = local_merge_defaults(defaults, overrides)
out = defaults;
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    value = overrides.(name);
    if isstruct(value) && isfield(defaults, name) && isstruct(defaults.(name)) && ...
            isscalar(value) && isscalar(defaults.(name))
        out.(name) = local_merge_defaults(defaults.(name), value);
    elseif ~isempty(value)
        out.(name) = value;
    else
        out.(name) = defaults.(name);
    end
end
end


function row = local_empty_row()
row = struct('dataset_stem', '', 'dataset_id', '', 'run_name', '', ...
    'run_tag', '', 'p10_tag', '', 'observable_mode', '', 'residual_form', '', ...
    'feature_name', '', 'method_tag', '', 'path_kind', '', 'n_components', NaN, ...
    'status', '', 'xcorr_status', '', 'activation_status', '', ...
    'roi_summary_status', '', 'message', '', 'bold_post_file', '', ...
    'dimred_result_file', '', 'density_source_names', '', 'n_density_sources', 0, ...
    'xcorr_mat_file', '', 'xcorr_peak_csv_file', '', 'xcorr_top_csv_file', '', ...
    'xcorr_by_density_dir', '', 'activation_dir', '', 'activation_group_dirs', '', ...
    'n_activation_groups', 0, 'roi_summary_dir', '', 'roi_summary_group_dirs', '', ...
    'n_roi_summary_groups', 0, 'n_requested_maps', 0, 'n_saved_maps', 0, ...
    'n_selected_roi_summary_modes', 0, 'n_saved_roi_summary_figures', 0, ...
    'runtime_sec', NaN);
end


function row = local_make_row(one_result)
row = local_empty_row();
names = fieldnames(row);
for i = 1:numel(names)
    name = names{i};
    if isfield(one_result, name) && ~isempty(one_result.(name))
        row.(name) = one_result.(name);
    end
end
end


function row = local_error_row(candidate, message)
row = local_empty_row();
row.dataset_stem = candidate.dataset_stem;
row.dataset_id = candidate.dataset_id;
row.run_name = candidate.run_name;
row.run_tag = candidate.run_tag;
row.p10_tag = candidate.p10_tag;
row.observable_mode = candidate.observable_mode;
row.residual_form = candidate.residual_form;
row.feature_name = candidate.feature_name;
row.method_tag = candidate.method_tag;
row.path_kind = candidate.path_kind;
row.n_components = candidate.n_components;
row.status = 'error';
row.xcorr_status = 'error';
row.activation_status = 'error';
row.roi_summary_status = 'error';
row.message = message;
row.bold_post_file = candidate.bold_post_file;
row.dimred_result_file = candidate.dimred_result_file;
row.xcorr_mat_file = candidate.xcorr_file;
row.xcorr_by_density_dir = fullfile(candidate.xcorr_dir, 'density');
end
