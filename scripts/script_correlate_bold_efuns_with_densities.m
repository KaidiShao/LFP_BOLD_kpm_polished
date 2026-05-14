% Canonical batch pipeline 8 entry.
%
% This script consumes pipeline 7 BOLD_POST artifacts and runs:
%   1) lagged BLP-BOLD cross-correlation
%   2) xcorr-ranked BOLD activation-map export
%
% Optional controls, set before run(...) if needed:
%   dataset_stems = {'e10gb1', 'e10gh1', 'e10fV1'};
%   observable_modes = {'HP_svd100'};
%   residual_forms = {'projected_kv'};
%   run_name_filter = {'mlp_obs_bold_wsl_20260423_e10gb1_projected_kv_HP_svd100'};
%   run_name_contains = {'HP_svd100', 'projected_kv'};
%   use_default_density_triplet = true;
%   density_source_kinds = {'event_density', 'raw_abs_density', 'dimred_abs_density', ...
%       'raw_complex_split_density', 'dimred_complex_split_density'};
%   blp_density_threshold_ratio = 0.70;
%   blp_dimred_method_tag = 'umap_k08';
%   require_all_density_sources = false;
%   output_mode = 'separate';  % 'combined' | 'separate' | 'both'
%   density_sources = struct([]);
%   force_recompute = false;
%   max_runs = [];
%   max_lag_sec = 10;
%   border_mask_sec = 10;
%   top_n = 5;
%   make_xcorr_figures = true;
%   make_activation_maps = true;
%   force_activation_redraw = false;
%   activation_export_combined = true;
%   activation_export_by_density = true;
%   activation_export_by_density_feature = true;
%   activation_value_mode = 'abs';
%   activation_slice_list = 1:20;
%   activation_tiles_per_row = 10;
%   activation_feature_reduce = 'mean';
%   activation_highlight_core_rois = true;
%   make_roi_summary = true;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

params = build_bold_cross_modal_coupling_params();
params.headless = true;
params.close_figures = true;
params.continue_on_error = true;

if ~exist('dataset_stems', 'var') || isempty(dataset_stems)
    dataset_stems = {'e10gb1', 'e10gh1', 'e10fV1'};
end
params.dataset_stems = cellstr(string(dataset_stems(:)).');

if exist('use_default_density_triplet', 'var') && ~isempty(use_default_density_triplet)
    params.use_default_density_triplet = logical(use_default_density_triplet);
end
if exist('density_source_kinds', 'var') && ~isempty(density_source_kinds)
    params.density_source_kinds = cellstr(string(density_source_kinds(:)).');
end
if exist('blp_density_threshold_ratio', 'var') && ~isempty(blp_density_threshold_ratio)
    params.blp_density_threshold_ratio = double(blp_density_threshold_ratio);
end
if exist('blp_dimred_method_tag', 'var') && ~isempty(blp_dimred_method_tag)
    params.blp_dimred_method_tag = char(string(blp_dimred_method_tag));
end
if exist('require_all_density_sources', 'var') && ~isempty(require_all_density_sources)
    params.require_all_density_sources = logical(require_all_density_sources);
end
if exist('output_mode', 'var') && ~isempty(output_mode)
    params = apply_bold_cross_modal_output_mode(params, output_mode);
end
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
if exist('density_sources', 'var') && ~isempty(density_sources)
    params.density_sources = density_sources;
end
if exist('force_recompute', 'var') && ~isempty(force_recompute)
    params.skip_existing_xcorr = ~logical(force_recompute);
    params.activation.skip_existing = ~logical(force_recompute);
    params.roi_summary.skip_existing = ~logical(force_recompute);
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
if exist('make_xcorr_figures', 'var') && ~isempty(make_xcorr_figures)
    params.xcorr.make_figures = logical(make_xcorr_figures);
end
if exist('make_activation_maps', 'var') && ~isempty(make_activation_maps)
    params.activation.enabled = logical(make_activation_maps);
end
if exist('make_roi_summary', 'var') && ~isempty(make_roi_summary)
    params.roi_summary.enabled = logical(make_roi_summary);
end
if exist('force_activation_redraw', 'var') && ~isempty(force_activation_redraw)
    params.activation.skip_existing = ~logical(force_activation_redraw);
end
if exist('force_roi_summary_redraw', 'var') && ~isempty(force_roi_summary_redraw)
    params.roi_summary.skip_existing = ~logical(force_roi_summary_redraw);
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

% Per-stage overrides still win if explicitly set after output_mode.
if exist('xcorr_export_combined', 'var') && ~isempty(xcorr_export_combined)
    params.xcorr.export_combined = logical(xcorr_export_combined);
end
if exist('xcorr_export_by_density', 'var') && ~isempty(xcorr_export_by_density)
    params.xcorr.export_by_density = logical(xcorr_export_by_density);
end
if exist('xcorr_export_by_density_feature', 'var') && ~isempty(xcorr_export_by_density_feature)
    params.xcorr.export_by_density_feature = logical(xcorr_export_by_density_feature);
end
if exist('activation_export_combined', 'var') && ~isempty(activation_export_combined)
    params.activation.export_combined = logical(activation_export_combined);
end
if exist('activation_export_by_density', 'var') && ~isempty(activation_export_by_density)
    params.activation.export_by_density = logical(activation_export_by_density);
end
if exist('activation_export_by_density_feature', 'var') && ~isempty(activation_export_by_density_feature)
    params.activation.export_by_density_feature = logical(activation_export_by_density_feature);
end
if exist('roi_summary_export_combined', 'var') && ~isempty(roi_summary_export_combined)
    params.roi_summary.export_combined = logical(roi_summary_export_combined);
end
if exist('roi_summary_export_by_density', 'var') && ~isempty(roi_summary_export_by_density)
    params.roi_summary.export_by_density = logical(roi_summary_export_by_density);
end
if exist('roi_summary_export_by_density_feature', 'var') && ~isempty(roi_summary_export_by_density_feature)
    params.roi_summary.export_by_density_feature = logical(roi_summary_export_by_density_feature);
end

fprintf('Running pipeline 8 batch for dataset(s): %s\n', strjoin(params.dataset_stems, ', '));
fprintf('Output mode: xcorr=%s | activation=%s\n', ...
    local_mode_label(params.xcorr.export_combined, params.xcorr.export_by_density), ...
    local_mode_label(params.activation.export_combined, params.activation.export_by_density));
if isempty(params.density_sources) && params.use_default_density_triplet
    fprintf('Density sources: default explicit set = %s\n', ...
        strjoin(params.density_source_kinds, ', '));
elseif ~isempty(params.density_sources)
    fprintf('Density sources: using explicit override (%d source(s)).\n', ...
        numel(params.density_sources));
end

manifest = process_bold_cross_modal_coupling_runs(params);


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
