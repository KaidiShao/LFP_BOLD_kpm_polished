% Canonical interactive BOLD eigenfunction/deconvolved-eigenfunction DR entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_bold_eigenfunction_reduction.m');
%
% Optional controls, set before run(...) if needed:
%   observable_modes = {'global_svd100'};
%   residual_forms = {'projected_vlambda'};
%   run_name_filter = {'mlp_obs_bold_...'};
%   run_name_contains = {'global_svd100'};
%   current_best_p7_only = true;
%   current_p7_autodl_roots = {'E:\autodl_results_local\bold_wsl'};
%   feature_names = {'efun_abs','efun_real','deconv_abs','deconv_real'};
%   method_filter = {'svd','logsvd','nmf','mds','umap'};
%   component_count_sweep = 3:20;
%   time_component_count_sweep = 3:20;
%   spectrum_component_count_sweep = 3:20;
%   max_modes = Inf;
%   force_recompute = false;
%   make_summary_plot = true;
%   max_runs = [];
%   continue_on_error = true;
%
% Role:
%   - canonical single-dataset entry for pipeline 9
%   - discovers BOLD_POST artifacts from pipeline 7
%   - runs P5-style DR methods on BOLD efun/deconv efun features
%   - saves one result per feature/method/component count
%   - writes a three-panel BOLD observable / BOLD efun / reduced component plot

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

cfg_name = char(string(cfg_name));
cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
params = build_bold_eigenfunction_reduction_params();
params.dataset_stems = {cfg.file_stem};

if exist('observable_modes', 'var') && ~isempty(observable_modes)
    params.observable_modes = cellstr(string(observable_modes(:)).');
end
if exist('residual_forms', 'var') && ~isempty(residual_forms)
    params.residual_forms = cellstr(string(residual_forms(:)).');
end
if exist('run_name_filter', 'var') && ~isempty(run_name_filter)
    params.run_name_filter = cellstr(string(run_name_filter(:)).');
end
if exist('run_name_contains', 'var') && ~isempty(run_name_contains)
    params.run_name_contains = cellstr(string(run_name_contains(:)).');
end
if exist('current_best_p7_only', 'var') && ~isempty(current_best_p7_only)
    params.current_best_p7_only = logical(current_best_p7_only);
end
if exist('current_p7_autodl_roots', 'var') && ~isempty(current_p7_autodl_roots)
    params.current_p7_autodl_roots = cellstr(string(current_p7_autodl_roots(:)).');
end
if exist('current_p7_allow_summary_only_outputs', 'var') && ~isempty(current_p7_allow_summary_only_outputs)
    params.current_p7_allow_summary_only_outputs = logical(current_p7_allow_summary_only_outputs);
end
if exist('feature_names', 'var') && ~isempty(feature_names)
    params.feature_names = cellstr(string(feature_names(:)).');
end
if exist('method_filter', 'var') && ~isempty(method_filter)
    params.method_filter = cellstr(string(method_filter(:)).');
end
if exist('component_count_sweep', 'var') && ~isempty(component_count_sweep)
    params.component_count_sweep = double(component_count_sweep(:)).';
end
if exist('time_component_count_sweep', 'var') && ~isempty(time_component_count_sweep)
    params.time_component_count_sweep = double(time_component_count_sweep(:)).';
end
if exist('spectrum_component_count_sweep', 'var') && ~isempty(spectrum_component_count_sweep)
    params.spectrum_component_count_sweep = double(spectrum_component_count_sweep(:)).';
end
if exist('max_modes', 'var') && ~isempty(max_modes)
    params.selection.max_modes = double(max_modes);
end
if exist('abs_thresh', 'var') && ~isempty(abs_thresh)
    params.selection.abs_thresh = double(abs_thresh);
end
if exist('sort_by', 'var') && ~isempty(sort_by)
    params.selection.sort_by = char(string(sort_by));
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing = ~logical(force_recompute);
end
if exist('make_summary_plot', 'var') && ~isempty(make_summary_plot)
    params.make_summary_plot = logical(make_summary_plot);
end
if exist('max_runs', 'var') && ~isempty(max_runs)
    params.max_runs = double(max_runs);
end
if exist('continue_on_error', 'var') && ~isempty(continue_on_error)
    params.continue_on_error = logical(continue_on_error);
end
if exist('plot_window_idx', 'var') && ~isempty(plot_window_idx)
    params.plot.window_idx = double(plot_window_idx(:));
end
if exist('plot_window_start', 'var') && ~isempty(plot_window_start)
    params.plot.window_start = double(plot_window_start);
end
if exist('plot_max_samples', 'var') && ~isempty(plot_max_samples)
    params.plot.max_samples = double(plot_max_samples);
end

fprintf('Running pipeline 9 BOLD efun/deconv DR for %s (%s)\n', ...
    cfg_name, cfg.dataset_id);
fprintf('Processed root:\n  %s\n', params.processed_root);
fprintf('Features: %s\n', strjoin(params.feature_names, ', '));
fprintf('Component sweep: %s\n', mat2str(params.component_count_sweep));

candidates = discover_completed_bold_eigenfunction_reduction_runs(params);
if isempty(candidates)
    if params.continue_on_error
        warning('No pipeline 7 BOLD_POST candidates were found for cfg %s. Skipping pipeline 9 for this cfg.', cfg_name);
        return;
    end
    error('No pipeline 7 BOLD_POST candidates were found for cfg %s.', cfg_name);
end

method_specs = build_bold_eigenfunction_reduction_method_specs(params);
if isempty(method_specs)
    error('No P9 reduction methods remain after applying method_filter.');
end

fprintf('Discovered %d BOLD_POST candidate run(s).\n', numel(candidates));
fprintf('Methods per feature: %d\n', numel(method_specs));

stage_root = io_project.get_pipeline_stage_dir( ...
    params.processed_root, cfg.file_stem, 9, 'bold_eigenfunction_reduction');
manifest_file = fullfile(stage_root, ...
    sprintf('%s_pipeline9_bold_eigenfunction_reduction_manifest.csv', cfg.file_stem));

rows = repmat(local_empty_row(), 0, 1);

for i_run = 1:numel(candidates)
    candidate = candidates(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s | %s\n', ...
        i_run, numel(candidates), candidate.dataset_stem, ...
        candidate.observable_mode, candidate.residual_form);
    fprintf('BOLD_POST:\n  %s\n', candidate.bold_post_file);

    for i_feature = 1:numel(params.feature_names)
        feature_name_i = params.feature_names{i_feature};

        for i_method = 1:numel(method_specs)
            spec = method_specs(i_method);
            if ~spec.enabled
                continue;
            end

            [method_cfg, method_tag, feature_tag] = ...
                build_bold_eigenfunction_reduction_cfg( ...
                candidate, feature_name_i, spec, params, repo_root);
            result_file = local_result_file(method_cfg);

            fprintf('  [%s | %s] -> %s\n', feature_tag, method_tag, method_cfg.output.root);
            t_start = tic;
            try
                if params.skip_existing && exist(result_file, 'file') == 2
                    row_i = local_make_row(candidate, feature_name_i, method_tag, ...
                        spec.kind, local_spec_components(spec), 'skipped_existing', ...
                        '', result_file, '', toc(t_start));
                else
                    [result, ~, ~] = run_bold_eigenfunction_reduction_pipeline(method_cfg);
                    row_i = local_make_row(candidate, feature_name_i, method_tag, ...
                        spec.kind, local_spec_components(spec), 'ok', '', ...
                        result.artifacts.result_mat_file, ...
                        result.artifacts.summary_png_file, toc(t_start));
                end
            catch ME
                message = getReport(ME, 'basic', 'hyperlinks', 'off');
                fprintf(2, '  [ERROR] %s | %s failed:\n%s\n', ...
                    feature_name_i, method_tag, message);
                row_i = local_make_row(candidate, feature_name_i, method_tag, ...
                    spec.kind, local_spec_components(spec), 'error', message, ...
                    result_file, '', toc(t_start));
                if ~params.continue_on_error
                    rethrow(ME);
                end
            end

            rows(end + 1, 1) = row_i; %#ok<SAGROW>
            if params.write_manifest
                local_write_manifest(rows, manifest_file);
            end
            if params.close_figures
                close all force;
            end
        end
    end
end

manifest = struct();
manifest.params = params;
manifest.candidates = candidates;
manifest.method_specs = method_specs;
manifest.file = manifest_file;
manifest.table = struct2table(rows);

if params.write_manifest
    local_write_manifest(rows, manifest_file);
end

fprintf('\nFinished pipeline 9 for %s.\n', cfg.dataset_id);
fprintf('Manifest:\n  %s\n', manifest_file);
fprintf('Rows: %d\n', height(manifest.table));


function row = local_empty_row()
row = struct();
row.dataset_stem = "";
row.dataset_id = "";
row.run_name = "";
row.run_tag = "";
row.observable_mode = "";
row.residual_form = "";
row.feature_name = "";
row.method_tag = "";
row.path_kind = "";
row.n_components = NaN;
row.status = "";
row.message = "";
row.bold_post_file = "";
row.result_mat_file = "";
row.summary_png_file = "";
row.runtime_sec = NaN;
end


function row = local_make_row(candidate, feature_name, method_tag, path_kind, ...
        n_components, status, message, result_file, summary_png, runtime_sec)
row = local_empty_row();
row.dataset_stem = string(candidate.dataset_stem);
row.dataset_id = string(candidate.dataset_id);
row.run_name = string(candidate.run_name);
row.run_tag = string(candidate.run_tag);
row.observable_mode = string(candidate.observable_mode);
row.residual_form = string(candidate.residual_form);
row.feature_name = string(feature_name);
row.method_tag = string(method_tag);
row.path_kind = string(path_kind);
row.n_components = double(n_components);
row.status = string(status);
row.message = string(message);
row.bold_post_file = string(candidate.bold_post_file);
row.result_mat_file = string(result_file);
row.summary_png_file = string(summary_png);
row.runtime_sec = double(runtime_sec);
end


function local_write_manifest(rows, manifest_file)
manifest_dir = fileparts(manifest_file);
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end
T = struct2table(rows);
writetable(T, manifest_file);
end


function result_file = local_result_file(method_cfg)
if isempty(method_cfg.save.tag)
    file_name = sprintf('%s.mat', method_cfg.save.file_stem);
else
    file_name = sprintf('%s_%s.mat', method_cfg.save.file_stem, method_cfg.save.tag);
end
result_file = fullfile(method_cfg.save.dir, file_name);
end


function k = local_spec_components(spec)
k = NaN;
if isfield(spec, 'options') && isstruct(spec.options) && ...
        isfield(spec.options, 'n_components') && ~isempty(spec.options.n_components)
    k = double(spec.options.n_components);
end
end
