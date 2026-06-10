% Canonical interactive BLP-BOLD cross-modal coupling entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_bold_cross_modal_coupling.m');
%
% Optional controls, set before run(...) if needed:
%   observable_modes = {'HP_svd100', 'global_svd100'};
%   residual_forms = {'projected_kv'};
%   run_name_filter = {'mlp_obs_bold_wsl_20260423_e10gb1_projected_kv_HP_svd100'};
%   run_name_contains = {'HP_svd100', 'projected_kv'};
%   current_best_p7_only = true;
%   current_p7_run_names = {'mlp_obs_bold_...'};  % explicit override
%   current_p7_autodl_roots = {'E:\autodl_results_local\bold_wsl'};
%   use_default_density_triplet = true;
%   density_source_kinds = {'event_density', 'raw_eigenfunction_density', 'dimred_eigenfunction_density'};
%   require_all_density_sources = false;
%   output_mode = 'separate';  % 'combined' | 'separate' | 'both'
%   density_sources = struct([]);
%   blp_dimred_method_tags = {'svd_k03', ..., 'umap_k08'};
%   force_recompute = false;
%   max_runs = [];
%   max_lag_sec = 10;
%   border_mask_sec = 10;
%   top_n = 5;
%   make_xcorr_figures = false;
%   make_activation_maps = false;
%   make_roi_summaries = false;
%   force_activation_redraw = false;
%   activation_export_combined = true;
%   activation_export_by_density = true;
%   activation_value_mode = 'abs';
%   activation_slice_list = 1:20;
%   activation_tiles_per_row = 10;
%   activation_feature_reduce = 'mean';
%   activation_highlight_core_rois = true;
%   headless = true;
%   close_figures_after_each_run = true;
%
% Role:
%   - canonical single-dataset entry for pipeline 8
%   - discovers BOLD_POST artifacts from pipeline 7 for one cfg
%   - runs the numeric xcorr stage by default; browse figures are opt-in
%   - leaves cfg / params / candidates / manifest and the most recent
%     run-level outputs in the workspace

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
set(groot, 'defaultFigureVisible', 'off');

%% -------------------------
%  Config and user options
%  -------------------------
cfg_name = char(string(cfg_name));
cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
params = build_bold_cross_modal_coupling_params();
params.dataset_stems = {cfg.file_stem};
params.headless = true;
params.close_figures = true;
params.continue_on_error = true;
params.load_existing_xcorr = true;

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
if exist('current_p7_run_names', 'var') && ~isempty(current_p7_run_names)
    params.current_p7_run_names = cellstr(string(current_p7_run_names(:)).');
end
if exist('current_p7_autodl_roots', 'var') && ~isempty(current_p7_autodl_roots)
    params.current_p7_autodl_roots = cellstr(string(current_p7_autodl_roots(:)).');
end
if exist('current_p7_allow_summary_only_outputs', 'var') && ~isempty(current_p7_allow_summary_only_outputs)
    params.current_p7_allow_summary_only_outputs = logical(current_p7_allow_summary_only_outputs);
end
if exist('use_default_density_triplet', 'var') && ~isempty(use_default_density_triplet)
    params.use_default_density_triplet = logical(use_default_density_triplet);
end
if exist('density_source_kinds', 'var') && ~isempty(density_source_kinds)
    params.density_source_kinds = cellstr(string(density_source_kinds(:)).');
end
if exist('require_all_density_sources', 'var') && ~isempty(require_all_density_sources)
    params.require_all_density_sources = logical(require_all_density_sources);
end
if exist('output_mode', 'var') && ~isempty(output_mode)
    params = apply_bold_cross_modal_output_mode(params, output_mode);
end
if exist('density_sources', 'var') && ~isempty(density_sources)
    params.density_sources = density_sources;
end
if exist('blp_dimred_method_tags', 'var') && ~isempty(blp_dimred_method_tags)
    params.blp_dimred_method_tags = cellstr(string(blp_dimred_method_tags(:)).');
end
if exist('blp_dimred_method_tag', 'var') && ~isempty(blp_dimred_method_tag)
    params.blp_dimred_method_tags = cellstr(string(blp_dimred_method_tag));
    params.blp_dimred_method_tag = params.blp_dimred_method_tags{1};
end
if exist('blp_density_condition_suffix', 'var') && ~isempty(blp_density_condition_suffix)
    params.blp_density_condition_suffix = char(string(blp_density_condition_suffix));
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing_xcorr = ~logical(force_recompute);
    params.activation.skip_existing = ~logical(force_recompute);
    params.roi_summary.skip_existing = ~logical(force_recompute);
    params.load_existing_xcorr = ~logical(force_recompute);
end
if exist('headless', 'var') && ~isempty(headless)
    params.headless = logical(headless);
    if params.headless
        set(groot, 'defaultFigureVisible', 'off');
    else
        set(groot, 'defaultFigureVisible', 'on');
    end
end
if exist('max_runs', 'var') && ~isempty(max_runs)
    params.max_runs = max_runs;
end
if exist('max_lag_sec', 'var') && ~isempty(max_lag_sec)
    params.xcorr.max_lag_sec = max_lag_sec;
end
if exist('border_mask_sec', 'var') && ~isempty(border_mask_sec)
    params.xcorr.border_mask_sec = border_mask_sec;
end
if exist('top_n', 'var') && ~isempty(top_n)
    params.xcorr.top_n = top_n;
    params.activation.top_n = top_n;
    params.roi_summary.top_n = top_n;
end
if exist('xcorr_save_tag', 'var') && ~isempty(xcorr_save_tag)
    params.xcorr.save_tag = char(string(xcorr_save_tag));
end
if exist('make_xcorr_figures', 'var') && ~isempty(make_xcorr_figures)
    params.xcorr.make_figures = logical(make_xcorr_figures);
end
if exist('make_activation_maps', 'var') && ~isempty(make_activation_maps)
    params.activation.enabled = logical(make_activation_maps);
end
if exist('make_roi_summaries', 'var') && ~isempty(make_roi_summaries)
    params.roi_summary.enabled = logical(make_roi_summaries);
end
if exist('force_activation_redraw', 'var') && ~isempty(force_activation_redraw)
    params.activation.skip_existing = ~logical(force_activation_redraw);
end
if exist('activation_export_combined', 'var') && ~isempty(activation_export_combined)
    params.activation.export_combined = logical(activation_export_combined);
end
if exist('activation_export_by_density', 'var') && ~isempty(activation_export_by_density)
    params.activation.export_by_density = logical(activation_export_by_density);
end
if exist('activation_value_mode', 'var') && ~isempty(activation_value_mode)
    params.activation.value_mode = char(string(activation_value_mode));
end
if exist('activation_slice_list', 'var') && ~isempty(activation_slice_list)
    params.activation.slice_list = activation_slice_list;
end
if exist('activation_tiles_per_row', 'var') && ~isempty(activation_tiles_per_row)
    params.activation.tiles_per_row = activation_tiles_per_row;
end
if exist('activation_feature_reduce', 'var') && ~isempty(activation_feature_reduce)
    params.activation.feature_reduce = char(string(activation_feature_reduce));
end
if exist('activation_highlight_core_rois', 'var') && ~isempty(activation_highlight_core_rois)
    params.activation.highlight_core_rois = logical(activation_highlight_core_rois);
end
if exist('xcorr_export_combined', 'var') && ~isempty(xcorr_export_combined)
    params.xcorr.export_combined = logical(xcorr_export_combined);
end
if exist('xcorr_export_by_density', 'var') && ~isempty(xcorr_export_by_density)
    params.xcorr.export_by_density = logical(xcorr_export_by_density);
end
if exist('activation_export_combined', 'var') && ~isempty(activation_export_combined)
    params.activation.export_combined = logical(activation_export_combined);
end
if exist('activation_export_by_density', 'var') && ~isempty(activation_export_by_density)
    params.activation.export_by_density = logical(activation_export_by_density);
end
if exist('close_figures_after_each_run', 'var') && ~isempty(close_figures_after_each_run)
    params.close_figures = logical(close_figures_after_each_run);
end

fprintf('Running pipeline 8 for %s (%s)\n', cfg_name, cfg.dataset_id);
fprintf('Processed root:\n  %s\n', params.processed_root);
fprintf('Output mode: xcorr=%s | activation=%s\n', ...
    local_mode_label(params.xcorr.export_combined, params.xcorr.export_by_density), ...
    local_mode_label(params.activation.export_combined, params.activation.export_by_density));
fprintf('Default outputs: xcorr_figures=%d | activation_maps=%d | roi_summaries=%d\n', ...
    params.xcorr.make_figures, params.activation.enabled, params.roi_summary.enabled);
fprintf('Current-best P7 only: %d\n', params.current_best_p7_only);
fprintf('Pipeline 8 xcorr dir:\n  %s\n', ...
    io_project.get_pipeline_stage_dir(params.processed_root, cfg, 8, ...
    'efun_density_cross_correlation'));
fprintf('Pipeline 8 map dir:\n  %s\n', ...
    io_project.get_pipeline_stage_dir(params.processed_root, cfg, 8, ...
    'figures_bold_top_xcorr_activation_maps'));

%% -------------------------
%  Step 1. Discover pipeline 7 outputs to consume
%  -------------------------
candidates = discover_completed_bold_cross_modal_coupling_runs(params);
if isempty(candidates)
    error('No pipeline 8 candidate runs were found for cfg %s.', cfg_name);
end

if ~isempty(params.max_runs)
    candidates = candidates(1:min(numel(candidates), params.max_runs));
end

fprintf('Discovered %d pipeline 8 candidate run(s).\n', numel(candidates));
for i = 1:numel(candidates)
    fprintf('  %s | %s | %s | %s\n', ...
        candidates(i).dataset_stem, candidates(i).observable_mode, ...
        candidates(i).residual_form, candidates(i).run_name);
end
if isempty(params.density_sources) && params.use_default_density_triplet
    fprintf('Density sources: default triplet = %s\n', ...
        strjoin(params.density_source_kinds, ', '));
elseif ~isempty(params.density_sources)
    fprintf('Density sources: using explicit override (%d source(s)).\n', ...
        numel(params.density_sources));
end

%% -------------------------
%  Step 2. Process each run
%  -------------------------
rows = repmat(local_empty_row(), numel(candidates), 1);
row_count = 0;

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

        xcorr_params = one_result.xcorr_params;
        activation_params = one_result.activation_params;
        xcorr_out = one_result.xcorr_out;
        act_result = one_result.act_result;
        xcorr_mat_file = one_result.xcorr_mat_file;
        activation_dir = one_result.activation_dir;

        fprintf('XCORR status: %s\n', one_result.xcorr_status);
        fprintf('Activation status: %s\n', one_result.activation_status);
        fprintf('XCORR MAT:\n  %s\n', xcorr_mat_file);
        if ~isempty(activation_dir)
            fprintf('Activation maps:\n  %s\n', activation_dir);
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
T = struct2table(rows, 'AsArray', true);
tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest = struct();
manifest.cfg = cfg;
manifest.params = params;
manifest.candidates = candidates;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(params.manifest_dir, ...
    ['pipeline8_bold_cross_modal_coupling_', lower(cfg.file_stem), '_', tag, '.csv']);
manifest.mat_file = fullfile(params.manifest_dir, ...
    ['pipeline8_bold_cross_modal_coupling_', lower(cfg.file_stem), '_', tag, '.mat']);

if exist(params.manifest_dir, 'dir') ~= 7
    mkdir(params.manifest_dir);
end
writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');

fprintf('\nPipeline 8 finished for cfg %s.\n', cfg_name);
fprintf('Manifest:\n  %s\n', manifest.csv_file);


function row = local_empty_row()
row = struct('dataset_stem', '', 'run_name', '', 'observable_mode', '', ...
    'residual_form', '', 'status', '', 'xcorr_status', '', ...
    'activation_status', '', 'message', '', 'bold_post_file', '', ...
    'density_source_names', '', 'n_density_sources', 0, ...
    'xcorr_mat_file', '', 'xcorr_by_density_dir', '', ...
    'activation_dir', '', 'activation_group_dirs', '', ...
    'n_activation_groups', 0, ...
    'n_requested_maps', 0, 'n_saved_maps', 0, 'runtime_sec', NaN);
end


function row = local_make_row(one_result)
row = local_empty_row();
row.dataset_stem = one_result.dataset_stem;
row.run_name = one_result.run_name;
row.observable_mode = one_result.observable_mode;
row.residual_form = one_result.residual_form;
row.status = one_result.status;
row.xcorr_status = one_result.xcorr_status;
row.activation_status = one_result.activation_status;
row.message = one_result.message;
row.bold_post_file = one_result.bold_post_file;
row.density_source_names = one_result.density_source_names;
row.n_density_sources = one_result.n_density_sources;
row.xcorr_mat_file = one_result.xcorr_mat_file;
row.xcorr_by_density_dir = one_result.xcorr_by_density_dir;
row.activation_dir = one_result.activation_dir;
row.activation_group_dirs = one_result.activation_group_dirs;
row.n_activation_groups = one_result.n_activation_groups;
row.n_requested_maps = one_result.n_requested_maps;
row.n_saved_maps = one_result.n_saved_maps;
row.runtime_sec = one_result.runtime_sec;
end


function row = local_error_row(candidate, message)
row = local_empty_row();
row.dataset_stem = candidate.dataset_stem;
row.run_name = candidate.run_name;
row.observable_mode = candidate.observable_mode;
row.residual_form = candidate.residual_form;
row.status = 'error';
row.xcorr_status = 'error';
row.activation_status = 'error';
row.message = message;
row.bold_post_file = candidate.bold_post_file;
row.density_source_names = '';
row.n_density_sources = 0;
row.xcorr_mat_file = candidate.xcorr_file;
row.xcorr_by_density_dir = fullfile(candidate.xcorr_dir, 'by_density');
end


function out = local_mode_label(export_combined, export_by_density)
if export_combined && export_by_density
    out = 'both';
elseif export_combined
    out = 'combined';
elseif export_by_density
    out = 'separate';
else
    out = 'none';
end
end
