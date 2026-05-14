function manifest = postprocess_bold_reskoopnet_results(params)
%POSTPROCESS_BOLD_RESKOOPNET_RESULTS Batch pipeline 7 postprocessing for completed BOLD runs.

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

runs = discover_completed_bold_reskoopnet_runs(params);
if isempty(runs)
    error('No completed BOLD ResKoopNet output folders were found.');
end

if exist(params.manifest_dir, 'dir') ~= 7
    mkdir(params.manifest_dir);
end

rows = repmat(local_empty_row(), numel(runs), 1);
row_count = 0;

fprintf('Discovered %d BOLD ResKoopNet runs.\n', numel(runs));

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s\n', i_run, numel(runs), ...
        run_info.dataset_stem, run_info.run_name);
    fprintf('Source:\n  %s\n', run_info.output_dir);

    t_run = tic;
    try
        one_result = run_one_bold_reskoopnet_postprocessing_core(run_info, params);
        if strcmp(one_result.status, 'skipped_existing')
            fprintf('Existing postprocessing found; skipping.\n');
        else
            fprintf('Saved BOLD postprocessing:\n  %s\n', one_result.post_file);
        end

        sync_one_pipeline7_postprocessing_summary_figures( ...
            run_info, one_result, params.processed_root);

        row_count = row_count + 1;
        rows(row_count) = local_make_row(run_info, one_result.status, ...
            one_result.message, one_result.post_file, ...
            one_result.main_png, one_result.deconv_png, ...
            one_result.timescale_png, ...
            local_get_field(one_result, 'intrinsic_activation_dir', ''), ...
            local_get_field(one_result, 'n_intrinsic_maps', 0), toc(t_run));
    catch ME
        message = getReport(ME, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s failed:\n%s\n', run_info.run_name, message);
        row_count = row_count + 1;
        rows(row_count) = local_make_row(run_info, 'error', message, ...
            '', '', '', '', '', 0, toc(t_run));
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
manifest.runs = runs;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(params.manifest_dir, ...
    ['pipeline7_bold_reskoopnet_postprocessing_manifest_', tag, '.csv']);
manifest.mat_file = fullfile(params.manifest_dir, ...
    ['pipeline7_bold_reskoopnet_postprocessing_manifest_', tag, '.mat']);

writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');
fprintf('\nBOLD ResKoopNet postprocessing finished.\n');
fprintf('Manifest:\n  %s\n', manifest.csv_file);
end


function params = local_apply_defaults(params)
defaults = build_bold_reskoopnet_postprocessing_params();
params = local_merge_defaults(defaults, params);
end


function out = local_merge_defaults(defaults, overrides)
out = defaults;
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    value = overrides.(name);
    if isstruct(value) && isfield(defaults, name) && isstruct(defaults.(name))
        out.(name) = local_merge_defaults(defaults.(name), value);
    elseif ~isempty(value)
        out.(name) = value;
    else
        out.(name) = defaults.(name);
    end
end
end


function row = local_empty_row()
row = struct('dataset_stem', '', 'run_name', '', 'observable_mode', '', ...
    'residual_form', '', 'status', '', 'message', '', 'post_file', '', ...
    'main_png', '', 'deconv_png', '', 'timescale_png', '', ...
    'intrinsic_activation_dir', '', 'n_intrinsic_maps', 0, 'runtime_sec', NaN);
end


function row = local_make_row(run_info, status, message, post_file, ...
        main_png, deconv_png, timescale_png, intrinsic_activation_dir, ...
        n_intrinsic_maps, runtime_sec)
row = local_empty_row();
row.dataset_stem = run_info.dataset_stem;
row.run_name = run_info.run_name;
row.observable_mode = local_get_field(run_info, 'observable_mode', '');
row.residual_form = local_get_field(run_info, 'residual_form', '');
row.status = status;
row.message = message;
row.post_file = post_file;
row.main_png = main_png;
row.deconv_png = deconv_png;
row.timescale_png = timescale_png;
row.intrinsic_activation_dir = intrinsic_activation_dir;
row.n_intrinsic_maps = n_intrinsic_maps;
row.runtime_sec = runtime_sec;
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
