% Canonical interactive BOLD ResKoopNet postprocessing entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m');
%
% Optional controls, set before run(...) if needed:
%   observable_modes = {'HP_svd100', 'global_svd100'};
%   residual_forms = {'projected_kv'};
%   run_name_filter = {'mlp_obs_bold_wsl_20260423_e10gb1_projected_kv_HP_svd100'};
%   run_name_contains = {'HP_svd100', 'projected_kv'};
%   autodl_roots = {'E:\autodl_results_local\bold_wsl'};
%   force_recompute = false;
%   timescale_only = false;
%   max_runs = [];
%   make_main_plot = true;
%   compute_deconv = true;
%   make_deconv_plot = true;
%   make_timescale_plot = true;
%   make_intrinsic_activation_maps = true;
%   make_intrinsic_roi_summary = true;
%   intrinsic_selection_mode = 'sorted';
%   intrinsic_basis_indices = 1:5;
%   intrinsic_value_mode = 'abs';
%   intrinsic_slice_list = 1:20;
%   intrinsic_tiles_per_row = 10;
%   intrinsic_feature_reduce = 'mean';
%   abs_thresh = 0.01;
%   max_basis = 80;
%   timescale_max_modes_all = Inf;
%   timescale_max_modes_sel = 80;
%   timescale_max_lag = [];
%   timescale_xlim_time = 240;
%   timescale_match_empirical_scale = true;
%   figure_export_method = 'print';  % 'print' or 'exportgraphics'
%   figure_renderer = 'opengl';
%   figure_resolution = 180;
%   close_figures_after_each_run = false;
%
% Role:
%   - canonical single-dataset entry for pipeline 7
%   - discovers completed BOLD ResKoopNet runs for one cfg
%   - expands the run-level steps in script form so intermediate variables
%     remain in the workspace
%   - exports BOLD_POST MAT files plus main / deconv / timescale figures
%   - leaves cfg / params / runs / manifest and the most recent run-level
%     intermediates in the workspace

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

%% -------------------------
%  Repo and plotting setup
%  -------------------------
this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

%% -------------------------
%  Config and user options
%  -------------------------
cfg_name = char(string(cfg_name));
cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
params = build_bold_reskoopnet_postprocessing_params();
params.dataset_stems = {cfg.file_stem};
params.exclude_dataset_stems = {};
params.headless = false;
params.close_figures = false;
params.continue_on_error = true;

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
if exist('autodl_roots', 'var') && ~isempty(autodl_roots)
    params.autodl_roots = cellstr(string(autodl_roots(:)).');
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing = ~logical(force_recompute);
end
if exist('make_main_plot', 'var') && ~isempty(make_main_plot)
    params.make_main_plot = logical(make_main_plot);
end
if exist('compute_deconv', 'var') && ~isempty(compute_deconv)
    params.compute_deconv = logical(compute_deconv);
end
if exist('make_deconv_plot', 'var') && ~isempty(make_deconv_plot)
    params.make_deconv_plot = logical(make_deconv_plot);
end
if exist('make_timescale_plot', 'var') && ~isempty(make_timescale_plot)
    params.make_timescale_plot = logical(make_timescale_plot);
else
    params.make_timescale_plot = true;
end
if exist('figure_export_method', 'var') && ~isempty(figure_export_method)
    params.figure_export_method = char(string(figure_export_method));
end
if exist('figure_renderer', 'var') && ~isempty(figure_renderer)
    params.figure_renderer = char(string(figure_renderer));
end
if exist('figure_resolution', 'var') && ~isempty(figure_resolution)
    params.resolution = double(figure_resolution);
end
if exist('make_intrinsic_activation_maps', 'var') && ~isempty(make_intrinsic_activation_maps)
    params.intrinsic_activation.enabled = logical(make_intrinsic_activation_maps);
end
if exist('make_intrinsic_roi_summary', 'var') && ~isempty(make_intrinsic_roi_summary)
    params.intrinsic_roi_summary.enabled = logical(make_intrinsic_roi_summary);
end
if exist('intrinsic_selection_mode', 'var') && ~isempty(intrinsic_selection_mode)
    params.intrinsic_activation.selection_mode = char(string(intrinsic_selection_mode));
end
if exist('intrinsic_basis_indices', 'var') && ~isempty(intrinsic_basis_indices)
    params.intrinsic_activation.basis_indices = intrinsic_basis_indices;
end
if exist('intrinsic_value_mode', 'var') && ~isempty(intrinsic_value_mode)
    params.intrinsic_activation.value_mode = char(string(intrinsic_value_mode));
end
if exist('intrinsic_slice_list', 'var') && ~isempty(intrinsic_slice_list)
    params.intrinsic_activation.slice_list = intrinsic_slice_list;
end
if exist('intrinsic_tiles_per_row', 'var') && ~isempty(intrinsic_tiles_per_row)
    params.intrinsic_activation.tiles_per_row = intrinsic_tiles_per_row;
end
if exist('intrinsic_feature_reduce', 'var') && ~isempty(intrinsic_feature_reduce)
    params.intrinsic_activation.feature_reduce = char(string(intrinsic_feature_reduce));
end
if exist('abs_thresh', 'var') && ~isempty(abs_thresh)
    params.post_opts.abs_thresh = abs_thresh;
end
if exist('max_basis', 'var') && ~isempty(max_basis)
    params.post_opts.max_basis = max_basis;
    params.timescale_max_modes_sel = max_basis;
end
if exist('timescale_max_modes_all', 'var') && ~isempty(timescale_max_modes_all)
    params.timescale_max_modes_all = double(timescale_max_modes_all);
end
if exist('timescale_max_modes_sel', 'var') && ~isempty(timescale_max_modes_sel)
    params.timescale_max_modes_sel = double(timescale_max_modes_sel);
end
if exist('timescale_max_lag', 'var')
    params.timescale_max_lag = double(timescale_max_lag);
elseif exist('timescale_maxLag', 'var')
    params.timescale_max_lag = double(timescale_maxLag);
end
if exist('timescale_xlim_time', 'var') && ~isempty(timescale_xlim_time)
    params.timescale_xlim_time = double(timescale_xlim_time);
end
if exist('timescale_match_empirical_scale', 'var') && ~isempty(timescale_match_empirical_scale)
    params.timescale_match_empirical_scale = logical(timescale_match_empirical_scale);
end
if exist('close_figures_after_each_run', 'var') && ~isempty(close_figures_after_each_run)
    params.close_figures = logical(close_figures_after_each_run);
end
if ~isfield(params, 'make_timescale_plot') || isempty(params.make_timescale_plot)
    params.make_timescale_plot = true;
end
if exist('timescale_only', 'var') && ~isempty(timescale_only) && logical(timescale_only)
    params.make_main_plot = false;
    params.compute_deconv = false;
    params.make_deconv_plot = false;
    params.make_timescale_plot = true;
    params.intrinsic_activation.enabled = false;
    params.intrinsic_roi_summary.enabled = false;
end
params = local_sync_timescale_struct(params);

fprintf('Running pipeline 7 for %s (%s)\n', cfg_name, cfg.dataset_id);
fprintf('AutoDL roots:\n');
for i_root = 1:numel(params.autodl_roots)
    fprintf('  %s\n', params.autodl_roots{i_root});
end
fprintf('Processed root:\n  %s\n', params.processed_root);
fprintf('Pipeline 7 output dir:\n  %s\n', ...
    io_project.get_pipeline_stage_dir(params.processed_root, cfg, 7, 'bold_postprocessing'));

%% -------------------------
%  Step 1. Discover completed runs for this cfg
%  -------------------------
runs = discover_completed_bold_reskoopnet_runs(params);
if isempty(runs)
    error('No completed BOLD ResKoopNet runs were found for cfg %s.', cfg_name);
end

if exist('max_runs', 'var') && ~isempty(max_runs)
    runs = runs(1:min(numel(runs), max_runs));
end

fprintf('Discovered %d completed BOLD run(s).\n', numel(runs));
for i = 1:numel(runs)
    fprintf('  %s | %s | %s | %s\n', ...
        runs(i).dataset_stem, runs(i).observable_mode, ...
        runs(i).residual_form, runs(i).run_name);
end

%% -------------------------
%  Step 2. Process each run
%  -------------------------
rows = repmat(local_empty_row(), numel(runs), 1);
row_count = 0;

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s | %s\n', ...
        i_run, numel(runs), run_info.dataset_stem, ...
        run_info.observable_mode, run_info.residual_form);
    fprintf('Run name:\n  %s\n', run_info.run_name);
    fprintf('Source output dir:\n  %s\n', run_info.output_dir);

    t_run = tic;
    try
        %% -------------------------
        %  Step 2.1. Load shared run inputs and resolve metadata
        %  -------------------------
        ctx = prepare_bold_reskoopnet_postprocessing_run(run_info, params);
        run_info = ctx.run_info;
        out_root = ctx.out_root;
        mat_dir = ctx.mat_dir;
        fig_dir = ctx.fig_dir;
        post_file = ctx.post_file;
        if ctx.skip_existing
            fprintf('Existing postprocessing found; refreshing cached figure outputs where needed.\n');
            [main_result, attempted_main] = backfill_bold_main_postprocessing_plot(post_file, params);
            main_png = '';
            if attempted_main && strcmp(local_get_field(main_result, 'status', ''), 'ok')
                main_png = local_get_field(main_result, 'main_png', '');
                if ~isempty(main_png)
                    fprintf('Refreshed main efun figure:\n  %s\n', main_png);
                end
            end
            [deconv_result, attempted_deconv] = backfill_bold_deconv_efun_plot(post_file, params);
            deconv_png = '';
            if attempted_deconv && strcmp(local_get_field(deconv_result, 'status', ''), 'ok')
                deconv_png = local_get_field(deconv_result, 'deconv_png', '');
                if ~isempty(deconv_png)
                    fprintf('Refreshed deconv efun figure:\n  %s\n', deconv_png);
                end
            end
            [timescale_result, attempted_timescale] = backfill_bold_timescale_plot(post_file, params);
            timescale_png = '';
            if attempted_timescale && strcmp(local_get_field(timescale_result, 'status', ''), 'ok')
                timescale_png = local_get_field(timescale_result, 'timescale_png', '');
                if ~isempty(timescale_png)
                    fprintf('Refreshed timescale figure:\n  %s\n', timescale_png);
                end
            end
            [act_result, attempted_backfill] = backfill_bold_intrinsic_activation_maps(post_file, params);
            intrinsic_activation_dir = '';
            n_intrinsic_maps = 0;
            row_status = 'skipped_existing';
            row_message = '';
            if attempted_deconv && ~isempty(deconv_png)
                row_status = 'ok_deconv_backfill';
            end
            if attempted_timescale && ~isempty(timescale_png)
                row_status = 'ok_timescale_backfill';
            end
            if attempted_backfill
                intrinsic_activation_dir = local_get_field(act_result, 'activation_dir', '');
                n_intrinsic_maps = local_get_field(act_result, 'n_requested_maps', 0);
                row_message = local_get_field(act_result, 'message', '');
                if strcmp(local_get_field(act_result, 'status', ''), 'ok')
                    row_status = 'ok_intrinsic_only';
                    fprintf('Backfilled intrinsic activation maps:\n  %s\n', intrinsic_activation_dir);
                else
                    fprintf('Intrinsic activation maps already present or not changed.\n');
                end
            end
            row_count = row_count + 1;
            rows(row_count) = local_make_row(run_info, row_status, ...
                row_message, post_file, main_png, deconv_png, timescale_png, intrinsic_activation_dir, ...
                n_intrinsic_maps, toc(t_run));
            continue;
        end

        EDMD_outputs = ctx.EDMD_outputs;
        concat_info = ctx.concat_info;
        source_info = ctx.source_info;
        observable_file = ctx.observable_file;
        obs_meta = ctx.obs_meta;
        session = ctx.session;
        dt = ctx.dt;
        time_vec = ctx.time_vec;
        session_border_t = ctx.session_border_t;
        fprintf('Loaded EDMD eigenfunctions: [%d x %d]\n', ...
            size(EDMD_outputs.efuns, 1), size(EDMD_outputs.efuns, 2));
        fprintf('Resolved observable file:\n  %s\n', observable_file);
        fprintf('Resolved dt: %.6g\n', dt);
        fprintf('Resolved sessions: %d\n', numel(session.session_start_idx));

        %% -------------------------
        %  Step 2.2. Run shared pipeline 7 stages
        %  -------------------------
        stage = run_bold_reskoopnet_postprocessing_stages(ctx, params);

        post_opts = stage.post_opts;
        fig_main = stage.fig_main;
        main_png = stage.main_png;
        EDMD_outputs = stage.EDMD_outputs;

        deconv_cfg = stage.deconv_cfg;
        fig_deconv = stage.fig_deconv;
        deconv = stage.deconv;
        deconv_png = stage.deconv_png;

        timescale_cfg = stage.timescale_cfg;
        fig_timescale = stage.fig_timescale;
        timescale_info = stage.timescale_info;
        timescale_png = stage.timescale_png;

        artifacts = stage.artifacts;
        BOLD_POST = stage.BOLD_POST;
        save_result = stage.save_result;
        act_result = stage.act_result;
        post_file = stage.save_result.post_file;
        intrinsic_activation_dir = stage.intrinsic_activation_dir;
        n_intrinsic_maps = stage.n_intrinsic_maps;

        if ~isempty(main_png)
            fprintf('Saved main efun figure:\n  %s\n', main_png);
        end
        if ~isempty(deconv_png)
            fprintf('Saved deconv figure:\n  %s\n', deconv_png);
        end
        if ~isempty(timescale_png)
            fprintf('Saved timescale figure:\n  %s\n', timescale_png);
        end
        fprintf('Saved BOLD_POST MAT:\n  %s\n', save_result.post_file);
        if ~isempty(intrinsic_activation_dir)
            fprintf('Saved intrinsic activation maps:\n  %s\n', intrinsic_activation_dir);
        end

        row_count = row_count + 1;
        rows(row_count) = local_make_row(run_info, 'ok', '', post_file, ...
            main_png, deconv_png, timescale_png, intrinsic_activation_dir, ...
            n_intrinsic_maps, stage.result.runtime_sec);
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

%% -------------------------
%  Step 3. Build run manifest in workspace
%  -------------------------
rows = rows(1:row_count);
manifest = struct();
manifest.cfg = cfg;
manifest.params = params;
manifest.runs = runs;
manifest.rows = rows;
manifest.table = struct2table(rows(:), 'AsArray', true);
manifest.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));

fprintf('\nFinished pipeline 7 for %s.\n', cfg.dataset_id);
fprintf('Processed rows: %d\n', height(manifest.table));
disp(manifest.table);


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
row.observable_mode = run_info.observable_mode;
row.residual_form = run_info.residual_form;
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


function params = local_sync_timescale_struct(params)
if ~isfield(params, 'timescale') || ~isstruct(params.timescale)
    params.timescale = struct();
end
params.timescale.max_modes_all = params.timescale_max_modes_all;
params.timescale.max_modes_sel = params.timescale_max_modes_sel;
params.timescale.maxLag = params.timescale_max_lag;
params.timescale.xlim_time = params.timescale_xlim_time;
params.timescale.match_empirical_scale_to_theoretical = ...
    params.timescale_match_empirical_scale;
end
