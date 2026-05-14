% Canonical interactive BLP eigenfunction reduction entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_blp_eigenfunction_reduction.m');
%
% Optional controls, set before run(...) if needed:
%   condition_key_filter = {'e10gh1|abs|projected_kv'};
%   run_name_filter = {'mlp_obs_e10fV1_local4cond_run01_projected_kv_abs'};
%   autodl_root = 'E:\autodl_results';
%   condition_output_tag = 'best20260507_abs_projected_vlambda';
%   condition_tag_mode = 'condition'; % 'condition', 'run_name', or 'condition_run_name'
%   method_filter = {'svd','umap'};
%   component_count_sweep = [2 4 8];
%   time_component_count_sweep = [2 4 8];
%   spectrum_component_count_sweep = [2 3 4];
%   make_overview_plot = false; % run-level overview is deprecated; top30 overview remains enabled
%   make_state_space_plot = true;
%   make_consensus_state_space_plot = true;
%   make_spectrum_diagnostics = true;
%   make_top30_window_plots = true;
%   top30_n_windows = 30;
%   top30_force_recompute_norm_cache = false;
%   make_thresholded_density = true;
%   make_thresholded_events = false;
%   make_dimred_thresholded_density = false;
%   make_dimred_thresholded_events = false;
%   threshold_mode = 'quantile';
%   threshold_ratio = 0.7;
%   threshold_ratio_sweep = [0.5 0.6 0.7 0.8 0.9];
%   continue_on_error = true;
%
% Role:
%   - canonical single-dataset entry for pipeline 5
%   - discovers completed runs for one cfg
%   - loads EDMD_outputs once per run
%   - reuses the preloaded workspace source across all reduction methods
%   - saves reduction results and figures under pipeline5_eigenfunction_reduction

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ', ...
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
params = build_blp_eigenfunction_reduction_params();
params.dataset_stems = {cfg.file_stem};
params.exclude_dataset_stems = {};

if exist('condition_key_filter', 'var') && ~isempty(condition_key_filter)
    params.condition_key_filter = cellstr(string(condition_key_filter(:)).');
end
if exist('run_name_filter', 'var') && ~isempty(run_name_filter)
    params.run_name_filter = cellstr(string(run_name_filter(:)).');
end
if exist('autodl_root', 'var') && ~isempty(autodl_root)
    params.autodl_root = char(string(autodl_root));
end
if exist('condition_output_tag', 'var') && ~isempty(condition_output_tag)
    params.condition_output_tag = char(string(condition_output_tag));
end
if exist('condition_tag_mode', 'var') && ~isempty(condition_tag_mode)
    params.condition_tag_mode = char(string(condition_tag_mode));
end
if exist('method_filter', 'var') && ~isempty(method_filter)
    params.method_filter = cellstr(string(method_filter(:)).');
end
if exist('component_count_sweep', 'var') && ~isempty(component_count_sweep)
    params.component_count_sweep = double(component_count_sweep(:).');
end
if exist('time_component_count_sweep', 'var') && ~isempty(time_component_count_sweep)
    params.time_component_count_sweep = double(time_component_count_sweep(:).');
end
if exist('spectrum_component_count_sweep', 'var') && ~isempty(spectrum_component_count_sweep)
    params.spectrum_component_count_sweep = double(spectrum_component_count_sweep(:).');
end
if exist('make_overview_plot', 'var') && ~isempty(make_overview_plot)
    params.make_overview_plot = logical(make_overview_plot);
end
if exist('make_state_space_plot', 'var') && ~isempty(make_state_space_plot)
    params.make_state_space_plot = logical(make_state_space_plot);
end
if exist('make_consensus_state_space_plot', 'var') && ~isempty(make_consensus_state_space_plot)
    params.make_consensus_state_space_plot = logical(make_consensus_state_space_plot);
end
if exist('make_spectrum_diagnostics', 'var') && ~isempty(make_spectrum_diagnostics)
    params.make_spectrum_diagnostics = logical(make_spectrum_diagnostics);
end
if exist('make_top30_window_plots', 'var') && ~isempty(make_top30_window_plots)
    params.make_top30_window_plots = logical(make_top30_window_plots);
end
if exist('top30_n_windows', 'var') && ~isempty(top30_n_windows)
    params.top30_n_windows = double(top30_n_windows);
end
if exist('top30_force_recompute_norm_cache', 'var') && ~isempty(top30_force_recompute_norm_cache)
    params.top30_force_recompute_norm_cache = logical(top30_force_recompute_norm_cache);
end
if exist('make_thresholded_density', 'var') && ~isempty(make_thresholded_density)
    params.do_thresholded_density = logical(make_thresholded_density);
end
if exist('make_thresholded_events', 'var') && ~isempty(make_thresholded_events)
    params.do_thresholded_events = logical(make_thresholded_events);
end
if exist('make_dimred_thresholded_density', 'var') && ~isempty(make_dimred_thresholded_density)
    params.do_dimred_thresholded_density = logical(make_dimred_thresholded_density);
end
if exist('make_dimred_thresholded_events', 'var') && ~isempty(make_dimred_thresholded_events)
    params.do_dimred_thresholded_events = logical(make_dimred_thresholded_events);
end
if exist('threshold_mode', 'var') && ~isempty(threshold_mode)
    params.threshold_mode = char(string(threshold_mode));
end
if exist('threshold_ratio', 'var') && ~isempty(threshold_ratio)
    params.threshold_ratio = double(threshold_ratio);
end
if exist('threshold_ratio_sweep', 'var') && ~isempty(threshold_ratio_sweep)
    params.threshold_ratio_sweep = double(threshold_ratio_sweep(:).');
end
if exist('save_thresholded_density_results', 'var') && ~isempty(save_thresholded_density_results)
    params.save_thresholded_density_results = logical(save_thresholded_density_results);
end
if exist('save_dimred_thresholded_density_results', 'var') && ~isempty(save_dimred_thresholded_density_results)
    params.save_dimred_thresholded_density_results = logical(save_dimred_thresholded_density_results);
end
if exist('continue_on_error', 'var') && ~isempty(continue_on_error)
    params.continue_on_error = logical(continue_on_error);
end

fprintf('Running BLP eigenfunction reduction for %s (%s)\n', cfg_name, cfg.dataset_id);
fprintf('AutoDL results root:\n  %s\n', params.autodl_root);
fprintf('Processed root:\n  %s\n', params.processed_root);

runs = discover_completed_blp_mlp_runs(params);
if isempty(runs)
    error('No completed MLP output runs were found under %s.', params.autodl_root);
end

method_specs = build_blp_eigenfunction_reduction_method_specs(params);
if isempty(method_specs)
    error('No eigenfunction reduction methods remain after applying method_filter.');
end

fprintf('Discovered %d completed unique MLP conditions.\n', numel(runs));
fprintf('Methods per run: %d\n', numel(method_specs));

C_consensus = [];
source_consensus_file = '';
if params.make_consensus_state_space_plot
    consensus_loader_cfg = struct();
    consensus_loader_cfg.file_stem = cfg.file_stem;
    [C_consensus, source_consensus_file] = io_results.load_consensus_state_results( ...
        consensus_loader_cfg, params.processed_root, []);
    fprintf('Consensus states:\n  %s\n', source_consensus_file);
end

manifest_rows = struct([]);
manifest_count = 0;
stage_root = io_project.get_pipeline_stage_dir( ...
    params.processed_root, cfg, 5, 'eigenfunction_reduction');
manifest_file = fullfile(stage_root, ...
    sprintf('%s_pipeline5_eigenfunction_reduction_manifest.csv', cfg.file_stem));

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n[%d/%d] %s | %s | %s\n', ...
        i_run, numel(runs), run_info.dataset, run_info.observable_mode, run_info.residual_form);
    fprintf('Source EDMD chunks:\n  %s\n', run_info.output_dir);

    fprintf('  Loading EDMD_outputs once for this run.\n');
    [EDMD_outputs, concat_info, source_info] = ...
        load_blp_eigenfunction_reduction_source(run_info, params);

    raw_cfg = build_blp_raw_eigenfunction_threshold_cfg( ...
        cfg, run_info, params, repo_root, ...
        EDMD_outputs, concat_info, source_info);
    raw_prep = prepare_blp_eigenfunction_reduction_inputs(EDMD_outputs, raw_cfg);
    raw_result = build_blp_raw_eigenfunction_threshold_result(raw_prep, run_info);

    if params.do_thresholded_density
        fprintf('  [raw] thresholded density -> %s\n', raw_cfg.thresholded_density.save_dir);
        run_blp_thresholded_density_stage(raw_prep, raw_result, raw_cfg.thresholded_density);
    end

    if params.do_thresholded_events
        fprintf('  [raw] thresholded events -> %s\n', raw_cfg.thresholded_events.save_dir);
        run_blp_thresholded_events_stage(raw_prep, raw_result, raw_cfg.thresholded_events);
    end

    for i_method = 1:numel(method_specs)
        spec = method_specs(i_method);
        if ~spec.enabled
            continue;
        end

        [method_cfg, method_tag] = build_blp_eigenfunction_reduction_cfg( ...
            cfg, run_info, spec, params, repo_root, ...
            EDMD_outputs, concat_info, source_info);

        fprintf('  [%d/%d] %s -> %s\n', ...
            i_method, numel(method_specs), method_tag, method_cfg.output.root);

        t_start = tic;
        try
            [result, ~, ~, ~] = run_eigenfunction_reduction_pipeline(method_cfg);
            plot_paths = plot_blp_eigenfunction_reduction_outputs( ...
                result, method_cfg, params, C_consensus);
            sync_pipeline5_eigenfunction_reduction_ssc_summary( ...
                method_cfg, run_info, method_tag, plot_paths, params.processed_root);

            if params.make_top30_window_plots
                plot_paths.top30_window_manifest_file = string( ...
                    run_blp_eigenfunction_dimred_top30_plots( ...
                    method_cfg, char(string(result.artifacts.result_mat_file)), ...
                    run_info, method_tag, params, repo_root));
            else
                plot_paths.top30_window_manifest_file = "";
            end

            row_i = build_blp_eigenfunction_reduction_manifest_row( ...
                run_info, method_tag, 'ok', '', result, plot_paths, toc(t_start));
            manifest_rows = local_append_manifest_row(manifest_rows, row_i);
            manifest_count = numel(manifest_rows);
        catch ME
            row_i = build_blp_eigenfunction_reduction_manifest_row( ...
                run_info, method_tag, 'error', ME.message, [], struct(), toc(t_start));
            manifest_rows = local_append_manifest_row(manifest_rows, row_i);
            manifest_count = numel(manifest_rows);
            fprintf(2, '  [ERROR] %s failed:\n    %s\n', method_tag, ME.message);
            if ~params.continue_on_error
                rethrow(ME);
            end
        end

        if params.write_manifest
            local_write_manifest(manifest_rows, manifest_file);
        end
    end
end

manifest = struct();
manifest.params = params;
manifest.runs = runs;
manifest.method_specs = method_specs;
manifest.source_consensus_file = source_consensus_file;
manifest.file = manifest_file;
manifest.table = local_rows_to_table(manifest_rows);

if params.write_manifest
    local_write_manifest(manifest_rows, manifest_file);
end

fprintf('\nFinished BLP eigenfunction reduction for %s.\n', cfg.dataset_id);
fprintf('Manifest:\n  %s\n', manifest_file);
fprintf('Rows: %d\n', height(manifest.table));


function local_write_manifest(rows, manifest_file)
T = local_rows_to_table(rows);
manifest_dir = fileparts(manifest_file);
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end
writetable(T, manifest_file);
end


function rows = local_append_manifest_row(rows, row)
if isempty(rows)
    rows = row;
else
    rows(end + 1) = row;
end
end


function T = local_rows_to_table(rows)
if isempty(rows)
    T = table();
    return;
end
T = struct2table(rows(:), 'AsArray', true);
end
