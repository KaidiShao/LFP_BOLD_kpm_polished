% Canonical interactive BLP eigenfunction postprocessing entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_blp_eigenfunction_postprocessing.m');
%
% Optional controls, set before run(...) if needed:
%   n_top_windows = 30;
%   window_length_samples = 6000;
%   condition_key_filter = {'e10gh1|abs|projected_kv'};
%   run_name_filter = {'mlp_obs_e10fV1_local4cond_run01_projected_kv_abs'};
%   autodl_root = 'E:\autodl_results';
%   output_root = [];
%   force_recompute = false;
%   skip_existing = true;
%   max_basis = 30;
%   timescale_max_modes_all = Inf;
%   timescale_max_modes_sel = 30;
%   timescale_max_lag = [];
%   timescale_xlim_time = 20;
%   timescale_match_empirical_scale = true;
%   timescale_similarity_max_samples = 100000;
%   deconv_max_modes_all = Inf;
%   deconv_max_modes_sel = 30;
%   spkt_cross_skip_existing = true;
%   mua_cross_skip_existing = true;
%   spkt_cross_save_under_run = false;
%   mua_cross_save_under_run = false;
%   make_main_plot = true;
%   make_timescale_plot = true;
%   make_deconv_plot = true;
%   make_deconv_window_norm_plot = true;
%   make_spkt_residual_cross_correlation = true;
%   make_mua_residual_cross_correlation = true;
%
% Role:
%   - canonical single-dataset entry for pipeline 6
%   - discovers completed runs for one cfg
%   - exports one run-level timescale diagnostic figure per run
%   - exports the three window-level top-window postprocessing figure types
%   - optionally runs SPKT and MUA residual cross-correlation
%   - leaves cfg / params / runs / manifest in the workspace
%
% Deferred follow-up:
%   - lagged SPKT remains optional and is not yet part of this canonical
%     runner; see docs/pipeline6_lagged_spkt_followup_2026-04-29.md

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
params = build_blp_eigenfunction_postprocessing_params();
params.dataset_stems = {cfg.file_stem};
params.exclude_dataset_stems = {};

if exist('n_top_windows', 'var') && ~isempty(n_top_windows)
    params.n_top_windows = n_top_windows;
end
if exist('window_length_samples', 'var') && ~isempty(window_length_samples)
    params.window_length_samples = window_length_samples;
end
if exist('condition_key_filter', 'var') && ~isempty(condition_key_filter)
    params.condition_key_filter = cellstr(string(condition_key_filter(:)).');
end
if exist('run_name_filter', 'var') && ~isempty(run_name_filter)
    params.run_name_filter = cellstr(string(run_name_filter(:)).');
end
if exist('autodl_root', 'var') && ~isempty(autodl_root)
    params.autodl_root = char(string(autodl_root));
end
if exist('output_root', 'var') && ~isempty(output_root)
    params.output_root = char(string(output_root));
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing = ~logical(force_recompute);
    params.spkt_cross_skip_existing = ~logical(force_recompute);
    params.mua_cross_skip_existing = ~logical(force_recompute);
end
if exist('skip_existing', 'var') && ~isempty(skip_existing)
    params.skip_existing = logical(skip_existing);
end
if exist('max_basis', 'var') && ~isempty(max_basis)
    params.max_basis = double(max_basis);
end
if exist('timescale_max_modes_all', 'var') && ~isempty(timescale_max_modes_all)
    params.timescale_max_modes_all = double(timescale_max_modes_all);
end
if exist('timescale_max_modes_sel', 'var') && ~isempty(timescale_max_modes_sel)
    params.timescale_max_modes_sel = double(timescale_max_modes_sel);
end
if exist('timescale_max_lag', 'var')
    params.timescale_max_lag = double(timescale_max_lag);
end
if exist('timescale_xlim_time', 'var') && ~isempty(timescale_xlim_time)
    params.timescale_xlim_time = double(timescale_xlim_time);
end
if exist('timescale_match_empirical_scale', 'var') && ~isempty(timescale_match_empirical_scale)
    params.timescale_match_empirical_scale = logical(timescale_match_empirical_scale);
end
if exist('timescale_similarity_max_samples', 'var') && ~isempty(timescale_similarity_max_samples)
    params.timescale_similarity_max_samples = double(timescale_similarity_max_samples);
end
if exist('deconv_max_modes_all', 'var') && ~isempty(deconv_max_modes_all)
    params.deconv_max_modes_all = double(deconv_max_modes_all);
end
if exist('deconv_max_modes_sel', 'var') && ~isempty(deconv_max_modes_sel)
    params.deconv_max_modes_sel = double(deconv_max_modes_sel);
end
if exist('spkt_cross_skip_existing', 'var') && ~isempty(spkt_cross_skip_existing)
    params.spkt_cross_skip_existing = logical(spkt_cross_skip_existing);
end
if exist('mua_cross_skip_existing', 'var') && ~isempty(mua_cross_skip_existing)
    params.mua_cross_skip_existing = logical(mua_cross_skip_existing);
end
if exist('spkt_cross_save_under_run', 'var') && ~isempty(spkt_cross_save_under_run)
    params.spkt_cross_save_under_run = logical(spkt_cross_save_under_run);
end
if exist('mua_cross_save_under_run', 'var') && ~isempty(mua_cross_save_under_run)
    params.mua_cross_save_under_run = logical(mua_cross_save_under_run);
end
if exist('make_main_plot', 'var') && ~isempty(make_main_plot)
    params.do_main_plot = logical(make_main_plot);
end
if exist('make_timescale_plot', 'var') && ~isempty(make_timescale_plot)
    params.do_timescale = logical(make_timescale_plot);
end
if exist('make_deconv_plot', 'var') && ~isempty(make_deconv_plot)
    params.do_deconv = logical(make_deconv_plot);
end
if exist('make_deconv_window_norm_plot', 'var') && ~isempty(make_deconv_window_norm_plot)
    params.do_deconv_window_norm = logical(make_deconv_window_norm_plot);
end
if exist('make_spkt_residual_cross_correlation', 'var') && ~isempty(make_spkt_residual_cross_correlation)
    params.do_spkt_cross_correlation = logical(make_spkt_residual_cross_correlation);
end
if exist('make_mua_residual_cross_correlation', 'var') && ~isempty(make_mua_residual_cross_correlation)
    params.do_mua_cross_correlation = logical(make_mua_residual_cross_correlation);
end

fprintf('Running BLP eigenfunction postprocessing for %s (%s)\n', cfg_name, cfg.dataset_id);
fprintf('AutoDL results root:\n  %s\n', params.autodl_root);
fprintf('Processed root:\n  %s\n', params.processed_root);
fprintf('Top-window output dir:\n  %s\n', ...
    io_project.get_pipeline_stage_dir(params.processed_root, cfg, 6, 'top_state_diversity_postprocessing'));

runs = discover_completed_blp_mlp_runs(params);
if isempty(runs)
    error('No completed MLP output runs were found under %s.', params.autodl_root);
end

fprintf('Discovered %d completed unique MLP conditions.\n', numel(runs));
for i = 1:numel(runs)
    fprintf('  %s | %s | %s | %s\n', ...
        runs(i).dataset, runs(i).observable_mode, runs(i).residual_form, runs(i).run_name);
end

window_rows = table();
spkt_rows = table();
mua_rows = table();
dataset_cache = struct();
need_window_outputs = params.do_main_plot || params.do_deconv || ...
    params.do_deconv_window_norm;
need_full_edmd = params.do_timescale || ...
    params.do_spkt_cross_correlation || params.do_mua_cross_correlation;

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n[%d/%d] Processing %s | %s | %s\n', ...
        i_run, numel(runs), run_info.dataset, run_info.observable_mode, run_info.residual_form);
    fprintf('Source EDMD chunks:\n  %s\n', run_info.output_dir);

    cfg_run = load_blp_pipeline_cfg_by_stem(run_info.dataset);
    dataset_key = matlab.lang.makeValidName(run_info.dataset);
    if need_window_outputs && ~isfield(dataset_cache, dataset_key)
        dataset_cache.(dataset_key) = load_blp_state_diversity_top_windows(cfg_run, params);
    end

    EDMD_full = [];
    residual_bundle = [];
    if need_full_edmd
        fprintf('  Loading full EDMD outputs into workspace.\n');
        EDMD_full = load_blp_run_edmd_outputs(run_info, params);
    end
    if params.do_spkt_cross_correlation || params.do_mua_cross_correlation
        fprintf('  Computing reusable residual bundle in workspace.\n');
        residual_bundle = prepare_blp_run_residual_workspace(EDMD_full, params);
    end

    timescale_png = export_blp_run_timescale_diagnostics(run_info, params, EDMD_full);
    if strlength(timescale_png) > 0
        fprintf('  Saved run-level timescale diagnostics:\n    %s\n', timescale_png);
    end

    run_table = table();
    run_manifest_file = "";
    if need_window_outputs
        [run_table, run_manifest_file] = export_blp_top_window_postprocess_bundle( ...
            run_info, dataset_cache.(dataset_key), params, EDMD_full);
        window_rows = [window_rows; run_table]; %#ok<AGROW>
        fprintf('  Saved top-window manifest:\n    %s\n', run_manifest_file);
        sync_pipeline6_dataset_window_figures( ...
            cfg_run, dataset_cache.(dataset_key), run_table, params.processed_root);
    else
        fprintf('  Skipping top-window postprocessing outputs; all window plot flags are false.\n');
    end

    run_output_root = fullfile( ...
        io_project.get_pipeline_stage_dir(params.processed_root, run_info.dataset, 6, 'top_state_diversity_postprocessing'), ...
        run_info.run_name);

    if params.do_spkt_cross_correlation
        spkt_save_dir = resolve_pipeline6_spkt_cross_save_dir(params, run_info, run_output_root);
        if exist(spkt_save_dir, 'dir') ~= 7
            mkdir(spkt_save_dir);
        end
        fprintf('  Full-time SPKT-residual cross-correlation:\n    %s\n', spkt_save_dir);

        existing_spkt_mat = '';
        if params.spkt_cross_skip_existing
            existing_spkt_mat = find_existing_pipeline6_cross_correlation_result( ...
                spkt_save_dir, cfg_run, run_info);
        end

        if ~isempty(existing_spkt_mat)
            fprintf('    Using existing result:\n      %s\n', existing_spkt_mat);
            spkt_result = load_saved_pipeline6_cross_correlation_result(existing_spkt_mat);
        else
            spkt_cmp_params = build_spkt_residual_cross_correlation_params( ...
                run_info, residual_bundle, spkt_save_dir, params);
            spkt_result = compute_spkt_residual_cross_correlation(cfg_run, spkt_cmp_params);
            spkt_result = save_spkt_residual_cross_correlation_result(spkt_result, spkt_save_dir);
        end

        spkt_overview_png = "";
        spkt_overview_fig = "";
        if params.spkt_cross_save_overview
            spkt_figure_dir = resolve_pipeline6_spkt_cross_figure_dir(params, run_info, run_output_root);
            [spkt_overview_png, spkt_overview_fig] = export_spkt_residual_cross_correlation_overview( ...
                spkt_result, cfg_run, run_info, spkt_figure_dir, params);
        end

        spkt_rows = [spkt_rows; build_spkt_residual_cross_correlation_row( ... %#ok<AGROW>
            run_info, spkt_save_dir, spkt_result, spkt_overview_png, spkt_overview_fig)];
    end

    if params.do_mua_cross_correlation
        if params.mua_cross_save_under_run
            mua_save_dir = fullfile(run_output_root, params.mua_cross_save_dir_name);
        else
            mua_save_dir = io_project.get_pipeline_stage_dir( ...
                params.processed_root, run_info.dataset, 6, 'mua_residual_cross_correlation');
        end
        if exist(mua_save_dir, 'dir') ~= 7
            mkdir(mua_save_dir);
        end
        fprintf('  Full-time MUA-proxy residual cross-correlation:\n    %s\n', mua_save_dir);

        existing_mua_mat = '';
        if params.mua_cross_skip_existing
            existing_mua_mat = find_existing_pipeline6_cross_correlation_result( ...
                mua_save_dir, cfg_run, run_info);
        end

        if ~isempty(existing_mua_mat)
            fprintf('    Using existing result:\n      %s\n', existing_mua_mat);
            mua_result = load_saved_pipeline6_cross_correlation_result(existing_mua_mat);
        else
            mua_cmp_params = build_mua_residual_cross_correlation_params( ...
                run_info, residual_bundle, mua_save_dir, params);
            mua_result = compute_mua_residual_cross_correlation(cfg_run, mua_cmp_params);
            mua_result = save_mua_residual_cross_correlation_result(mua_result, mua_save_dir);
        end

        mua_overview_png = "";
        mua_overview_fig = "";
        if params.mua_cross_save_overview
            if params.mua_cross_save_under_run
                mua_figure_dir = fullfile(run_output_root, ...
                    io_project.get_pipeline_stage_name(6, 'figures_mua_residual_cross_correlation'));
            else
                mua_figure_dir = io_project.get_pipeline_stage_dir( ...
                    params.processed_root, run_info.dataset, 6, 'figures_mua_residual_cross_correlation');
            end
            [mua_overview_png, mua_overview_fig] = export_mua_residual_cross_correlation_overview( ...
                mua_result, cfg_run, mua_figure_dir, params);
        end

        mua_rows = [mua_rows; build_mua_residual_cross_correlation_row( ... %#ok<AGROW>
            run_info, mua_save_dir, mua_result, mua_overview_png, mua_overview_fig)];
    end
end

manifest = struct();
manifest.params = params;
manifest.runs = runs;
manifest.table = window_rows;
manifest.spkt_table = spkt_rows;
manifest.mua_table = mua_rows;
manifest.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fprintf('\nFinished BLP eigenfunction postprocessing for %s.\n', cfg.dataset_id);
fprintf('Window rows: %d\n', height(manifest.table));
fprintf('SPKT rows: %d\n', height(manifest.spkt_table));
fprintf('MUA rows: %d\n', height(manifest.mua_table));
