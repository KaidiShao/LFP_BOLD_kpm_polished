function manifest = process_bold_cross_modal_coupling_runs(params)
%PROCESS_BOLD_CROSS_MODAL_COUPLING_RUNS Batch pipeline 8 cross-modal coupling.

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

candidates = discover_completed_bold_cross_modal_coupling_runs(params);
if isempty(candidates)
    error('No pipeline 8 candidate BOLD_POST files were found.');
end

if ~isempty(params.max_runs)
    candidates = candidates(1:min(numel(candidates), params.max_runs));
end

if exist(params.manifest_dir, 'dir') ~= 7
    mkdir(params.manifest_dir);
end

rows = repmat(local_empty_row(), numel(candidates), 1);
row_count = 0;

fprintf('Discovered %d pipeline 8 candidate run(s).\n', numel(candidates));

for i_run = 1:numel(candidates)
    candidate = candidates(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s | %s\n', ...
        i_run, numel(candidates), candidate.dataset_stem, ...
        candidate.observable_mode, candidate.residual_form);
    fprintf('Run name:\n  %s\n', candidate.run_name);
    fprintf('BOLD_POST:\n  %s\n', candidate.bold_post_file);

    try
        one_result = run_one_bold_cross_modal_coupling_core(candidate, params);
        fprintf('XCORR status: %s\n', one_result.xcorr_status);
        fprintf('Activation status: %s\n', one_result.activation_status);
        fprintf('ROI summary status: %s\n', one_result.roi_summary_status);
        fprintf('XCORR MAT:\n  %s\n', one_result.xcorr_mat_file);
        if ~isempty(one_result.activation_dir)
            fprintf('Activation maps:\n  %s\n', one_result.activation_dir);
        end
        if ~isempty(one_result.roi_summary_dir)
            fprintf('ROI summaries:\n  %s\n', one_result.roi_summary_dir);
        end

        row_count = row_count + 1;
        rows(row_count) = local_make_row(one_result);
    catch ME
        message = getReport(ME, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s failed:\n%s\n', candidate.run_name, message);
        row_count = row_count + 1;
        rows(row_count) = local_error_row(candidate, message);
        if ~params.continue_on_error
            rethrow(ME);
        end
    end

    if params.close_figures
        close all force;
    end
end

rows = rows(1:row_count);
T = struct2table(rows);
tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest = struct();
manifest.params = params;
manifest.candidates = candidates;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(params.manifest_dir, ...
    ['pipeline8_xcorr_manifest_', tag, '.csv']);
manifest.mat_file = fullfile(params.manifest_dir, ...
    ['pipeline8_xcorr_manifest_', tag, '.mat']);

writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');
fprintf('\nPipeline 8 cross-modal coupling finished.\n');
fprintf('Manifest:\n  %s\n', manifest.csv_file);
end


function params = local_apply_defaults(params)
defaults = build_bold_cross_modal_coupling_params();
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
row = struct('dataset_stem', '', 'run_name', '', 'run_tag', '', 'observable_mode', '', ...
    'residual_form', '', 'status', '', 'xcorr_status', '', ...
    'activation_status', '', 'roi_summary_status', '', 'message', '', 'bold_post_file', '', ...
    'density_source_names', '', 'n_density_sources', 0, ...
    'xcorr_mat_file', '', 'xcorr_peak_csv_file', '', ...
    'xcorr_top_csv_file', '', 'xcorr_by_density_dir', '', ...
    'summary_sync_status', '', 'summary_dir', '', 'n_summary_figures', 0, ...
    'activation_dir', '', 'activation_group_dirs', '', ...
    'n_activation_groups', 0, 'roi_summary_dir', '', 'roi_summary_group_dirs', '', ...
    'n_roi_summary_groups', 0, ...
    'n_requested_maps', 0, 'n_saved_maps', 0, ...
    'n_selected_roi_summary_modes', 0, 'n_saved_roi_summary_figures', 0, ...
    'runtime_sec', NaN);
end


function row = local_make_row(one_result)
row = local_empty_row();
row.dataset_stem = one_result.dataset_stem;
row.run_name = one_result.run_name;
row.run_tag = local_get_field(one_result, 'run_tag', '');
row.observable_mode = one_result.observable_mode;
row.residual_form = one_result.residual_form;
row.status = one_result.status;
row.xcorr_status = one_result.xcorr_status;
row.activation_status = one_result.activation_status;
row.roi_summary_status = one_result.roi_summary_status;
row.message = one_result.message;
row.bold_post_file = one_result.bold_post_file;
row.density_source_names = one_result.density_source_names;
row.n_density_sources = one_result.n_density_sources;
row.xcorr_mat_file = one_result.xcorr_mat_file;
row.xcorr_peak_csv_file = one_result.xcorr_peak_csv_file;
row.xcorr_top_csv_file = one_result.xcorr_top_csv_file;
row.xcorr_by_density_dir = one_result.xcorr_by_density_dir;
row.summary_sync_status = local_get_field(one_result, 'summary_sync_status', '');
row.summary_dir = local_get_field(one_result, 'summary_dir', '');
row.n_summary_figures = local_get_field(one_result, 'n_summary_figures', 0);
row.activation_dir = one_result.activation_dir;
row.activation_group_dirs = one_result.activation_group_dirs;
row.n_activation_groups = one_result.n_activation_groups;
row.roi_summary_dir = local_get_field(one_result, 'roi_summary_dir', '');
row.roi_summary_group_dirs = local_get_field(one_result, 'roi_summary_group_dirs', '');
row.n_roi_summary_groups = local_get_field(one_result, 'n_roi_summary_groups', 0);
row.n_requested_maps = one_result.n_requested_maps;
row.n_saved_maps = one_result.n_saved_maps;
row.n_selected_roi_summary_modes = local_get_field(one_result, 'n_selected_roi_summary_modes', 0);
row.n_saved_roi_summary_figures = local_get_field(one_result, 'n_saved_roi_summary_figures', 0);
row.runtime_sec = one_result.runtime_sec;
end


function row = local_error_row(candidate, message)
row = local_empty_row();
row.dataset_stem = candidate.dataset_stem;
row.run_name = candidate.run_name;
row.run_tag = local_get_field(candidate, 'run_tag', '');
row.observable_mode = candidate.observable_mode;
row.residual_form = candidate.residual_form;
row.status = 'error';
row.xcorr_status = 'error';
row.activation_status = 'error';
row.roi_summary_status = 'error';
row.message = message;
row.bold_post_file = candidate.bold_post_file;
row.density_source_names = '';
row.n_density_sources = 0;
row.xcorr_mat_file = candidate.xcorr_file;
row.xcorr_by_density_dir = fullfile(candidate.xcorr_dir, 'density');
row.summary_sync_status = '';
row.summary_dir = '';
row.n_summary_figures = 0;
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
