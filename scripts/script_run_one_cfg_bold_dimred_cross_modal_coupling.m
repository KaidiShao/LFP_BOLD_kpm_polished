% Canonical interactive dimred-BOLD/LFP-density coupling entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_run_one_cfg_bold_dimred_cross_modal_coupling.m');
%
% Optional controls:
%   observable_modes = {'global_svd100'};
%   residual_forms = {'projected_vlambda'};
%   run_name_filter = {'mlp_obs_bold_...'};
%   run_name_contains = {'global_svd100'};
%   current_best_p7_only = true;
%   current_p7_run_names = {'mlp_obs_bold_...'};  % explicit override
%   current_p7_autodl_roots = {'E:\autodl_results_local\bold_wsl'};
%   feature_names = {'efun_abs','deconv_abs'};
%   method_tags = {'svd_k05','umap_k08'};
%   component_counts = 3:20;
%   use_default_density_triplet = true;
%   density_source_kinds = {'event_density','raw_abs_density','dimred_abs_density'};
%   blp_dimred_method_tags = {'svd_k03', ..., 'umap_k08'};
%   output_mode = 'both'; % 'combined' | 'separate' | 'both'
%   max_lag_sec = 10;
%   border_mask_sec = 10;
%   top_n = 5;
%   make_xcorr_figures = false;
%   make_activation_maps = false;
%   make_roi_summaries = false;
%   force_recompute = false;
%   max_runs = [];
%
% Role:
%   - pipeline 10 entry
%   - consumes pipeline 9 BOLD efun/deconv dimred results
%   - compares dimred BOLD components with LFP/BLP density sources
%   - runs numeric xcorr by default; component-ranked maps/summaries are opt-in

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
params = build_bold_dimred_cross_modal_coupling_params();
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
if exist('current_p7_run_names', 'var') && ~isempty(current_p7_run_names)
    params.current_p7_run_names = cellstr(string(current_p7_run_names(:)).');
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
if exist('method_tags', 'var') && ~isempty(method_tags)
    params.method_tags = cellstr(string(method_tags(:)).');
end
if exist('path_kinds', 'var') && ~isempty(path_kinds)
    params.path_kinds = cellstr(string(path_kinds(:)).');
end
if exist('component_counts', 'var') && ~isempty(component_counts)
    params.component_counts = double(component_counts(:)).';
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
if exist('output_mode', 'var') && ~isempty(output_mode)
    params = local_apply_output_mode(params, output_mode);
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing_xcorr = ~logical(force_recompute);
    params.activation.skip_existing = ~logical(force_recompute);
    params.roi_summary.skip_existing = ~logical(force_recompute);
    params.load_existing_xcorr = ~logical(force_recompute);
end
if exist('max_runs', 'var') && ~isempty(max_runs)
    params.max_runs = double(max_runs);
end
if exist('max_lag_sec', 'var') && ~isempty(max_lag_sec)
    params.xcorr.max_lag_sec = double(max_lag_sec);
end
if exist('border_mask_sec', 'var') && ~isempty(border_mask_sec)
    params.xcorr.border_mask_sec = double(border_mask_sec);
end
if exist('top_n', 'var') && ~isempty(top_n)
    params.xcorr.top_n = double(top_n);
    params.activation.top_n = double(top_n);
    params.roi_summary.top_n = double(top_n);
end
if exist('xcorr_save_tag', 'var') && ~isempty(xcorr_save_tag)
    params.xcorr.save_tag = char(string(xcorr_save_tag));
end
if exist('component_value_modes', 'var') && ~isempty(component_value_modes)
    params.xcorr.component_value_modes = cellstr(string(component_value_modes(:)).');
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
if exist('continue_on_error', 'var') && ~isempty(continue_on_error)
    params.continue_on_error = logical(continue_on_error);
end

fprintf('Running pipeline 10 dimred-BOLD/LFP coupling for %s (%s)\n', ...
    cfg_name, cfg.dataset_id);
fprintf('Processed root:\n  %s\n', params.processed_root);
fprintf('Feature filters: %s\n', strjoin(params.feature_names, ', '));
fprintf('Method filters: %s\n', strjoin(params.method_tags, ', '));
fprintf('Default outputs: xcorr_figures=%d | activation_maps=%d | roi_summaries=%d\n', ...
    params.xcorr.make_figures, params.activation.enabled, params.roi_summary.enabled);
fprintf('Current-best P7 only: %d\n', params.current_best_p7_only);

manifest = process_bold_dimred_cross_modal_coupling_runs(params);


function params = local_apply_output_mode(params, output_mode)
mode = lower(strtrim(char(string(output_mode))));
switch mode
    case 'combined'
        export_combined = true;
        export_by_density = false;
    case 'separate'
        export_combined = false;
        export_by_density = true;
    case 'both'
        export_combined = true;
        export_by_density = true;
    otherwise
        error('Unsupported pipeline 10 output_mode: %s', output_mode);
end
params.output_mode = mode;
params.xcorr.export_combined = export_combined;
params.xcorr.export_by_density = export_by_density;
params.activation.export_combined = export_combined;
params.activation.export_by_density = export_by_density;
params.roi_summary.export_combined = export_combined;
params.roi_summary.export_by_density = export_by_density;
end
