function result = run_one_bold_dimred_cross_modal_coupling_core(candidate, params)
%RUN_ONE_BOLD_DIMRED_CROSS_MODAL_COUPLING_CORE Run pipeline 10 for one P9 result.

if nargin < 1 || isempty(candidate)
    error('candidate is required.');
end
if nargin < 2 || isempty(params)
    params = build_bold_dimred_cross_modal_coupling_params();
else
    params = local_apply_defaults(params);
end

t_run = tic;
density_sources = params.density_sources;
if isempty(density_sources) && params.use_default_density_triplet
    density_sources = resolve_bold_cross_modal_density_sources(candidate, params);
end
[density_source_names, n_density_sources] = local_density_source_summary(density_sources);

xcorr_params = params.xcorr;
xcorr_params.save_dir = candidate.xcorr_dir;
xcorr_params.save_tag = local_get_field(xcorr_params, 'save_tag', 'dimred_xcorr');
candidate.xcorr_file = fullfile(candidate.xcorr_dir, [xcorr_params.save_tag, '.mat']);

xcorr_recomputed = false;
if params.skip_existing_xcorr && exist(candidate.xcorr_file, 'file') == 2 && ...
        local_saved_xcorr_matches_density_sources(candidate.xcorr_file, density_sources, xcorr_params)
    xcorr_status = 'skipped_existing';
    xcorr_input = candidate.xcorr_file;
    if params.load_existing_xcorr
        xcorr_out = load_bold_xcorr_output(candidate.xcorr_file);
    else
        xcorr_out = struct();
    end
else
    xcorr_out = compute_bold_dimred_efun_density_cross_correlation( ...
        candidate.dimred_result_file, density_sources, xcorr_params);
    xcorr_input = xcorr_out;
    xcorr_status = 'ok';
    xcorr_recomputed = true;
end

act_result = struct();
activation_status = 'disabled';
activation_params = params.activation;
if isfield(params, 'activation') && isstruct(params.activation) && ...
        isfield(params.activation, 'enabled') && params.activation.enabled
    activation_params.processed_root = params.processed_root;
    activation_params.datapons_root = params.datapons_root;
    activation_params.output_root = candidate.activation_root;
    if xcorr_recomputed
        activation_params.skip_existing = false;
    end
    act_result = export_bold_dimred_top_xcorr_activation_maps( ...
        candidate.dimred_result_file, candidate.bold_post_file, xcorr_input, activation_params);
    activation_status = act_result.status;
end

roi_summary_result = struct();
roi_summary_status = 'disabled';
roi_summary_params = local_get_field(params, 'roi_summary', struct());
roi_summary_enabled = isstruct(roi_summary_params) && ...
    isfield(roi_summary_params, 'enabled') && roi_summary_params.enabled;
if roi_summary_enabled
    roi_summary_params.processed_root = params.processed_root;
    roi_summary_params.datapons_root = params.datapons_root;
    roi_summary_params.output_root = candidate.activation_root;
    if xcorr_recomputed
        roi_summary_params.skip_existing = false;
    end
    roi_summary_result = export_bold_dimred_top_xcorr_roi_bar_summaries( ...
        candidate.dimred_result_file, candidate.bold_post_file, xcorr_input, roi_summary_params);
    roi_summary_status = roi_summary_result.status;
end

result = struct();
result.candidate = candidate;
result.dataset_stem = candidate.dataset_stem;
result.dataset_id = candidate.dataset_id;
result.run_name = candidate.run_name;
result.run_tag = candidate.run_tag;
result.p10_tag = candidate.p10_tag;
result.observable_mode = candidate.observable_mode;
result.residual_form = candidate.residual_form;
result.feature_name = candidate.feature_name;
result.method_tag = candidate.method_tag;
result.path_kind = candidate.path_kind;
result.n_components = candidate.n_components;
result.bold_post_file = candidate.bold_post_file;
result.dimred_result_file = candidate.dimred_result_file;
result.xcorr_dir = candidate.xcorr_dir;
result.xcorr_mat_file = candidate.xcorr_file;
result.xcorr_peak_csv_file = fullfile(candidate.xcorr_dir, [xcorr_params.save_tag, '_peaks.csv']);
result.xcorr_top_csv_file = fullfile(candidate.xcorr_dir, [xcorr_params.save_tag, '_top.csv']);
result.xcorr_by_density_dir = fullfile(candidate.xcorr_dir, 'density');
result.activation_dir = local_get_field(act_result, 'activation_dir', '');
result.activation_info_file = local_get_field(act_result, 'info_file', '');
result.activation_group_dirs = local_join_group_field(act_result, 'activation_dir');
result.activation_info_files = local_join_group_field(act_result, 'info_file');
result.n_activation_groups = local_count_group_results(act_result);
result.roi_summary_dir = local_get_field(roi_summary_result, 'summary_dir', '');
result.roi_summary_info_file = local_get_field(roi_summary_result, 'info_file', '');
result.roi_summary_group_dirs = local_join_group_field(roi_summary_result, 'summary_dir');
result.roi_summary_info_files = local_join_group_field(roi_summary_result, 'info_file');
result.n_roi_summary_groups = local_count_group_results(roi_summary_result);
result.density_sources = density_sources;
result.density_source_names = density_source_names;
result.n_density_sources = n_density_sources;
result.xcorr_status = xcorr_status;
result.activation_status = activation_status;
result.roi_summary_status = roi_summary_status;
result.status = local_combine_status(xcorr_status, activation_status, params.activation.enabled);
result.status = local_combine_status(result.status, roi_summary_status, roi_summary_enabled);
result.message = local_compose_message(xcorr_status, activation_status, act_result, ...
    roi_summary_status, roi_summary_result);
result.n_requested_maps = local_get_field(act_result, 'n_requested_maps', 0);
result.n_saved_maps = local_get_field(act_result, 'n_saved_maps', 0);
result.n_selected_roi_summary_modes = local_get_field(roi_summary_result, 'n_selected_modes', 0);
result.n_saved_roi_summary_figures = local_get_field(roi_summary_result, 'n_saved_figures', 0);
result.runtime_sec = toc(t_run);
result.xcorr_params = xcorr_params;
result.activation_params = activation_params;
result.roi_summary_params = roi_summary_params;
if local_get_field(params, 'return_heavy_outputs', false)
    result.xcorr_out = xcorr_out;
    result.act_result = act_result;
    result.roi_summary_result = roi_summary_result;
end
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


function [names_text, n] = local_density_source_summary(density_sources)
if isempty(density_sources)
    names_text = '';
    n = 0;
    return;
end
if iscell(density_sources)
    sources = density_sources;
else
    sources = num2cell(density_sources(:));
end
names = cell(1, numel(sources));
for i = 1:numel(sources)
    names{i} = local_get_field(sources{i}, 'name', sprintf('density%d', i));
end
names_text = strjoin(names, ';');
n = numel(sources);
end


function tf = local_saved_xcorr_matches_density_sources(xcorr_file, density_sources, xcorr_params)
tf = false;
if local_saved_xcorr_is_oversized_for_current_policy(xcorr_file, xcorr_params)
    return;
end
try
    out = load_bold_xcorr_output(xcorr_file);
catch
    return;
end
saved_sources = local_get_field(out, 'density_sources', {});
tf = local_density_source_names(saved_sources) == local_density_source_names(density_sources);
end


function tf = local_saved_xcorr_is_oversized_for_current_policy(xcorr_file, xcorr_params)
tf = false;
if ~isstruct(xcorr_params) || ~isfield(xcorr_params, 'save_source_results') || ...
        xcorr_params.save_source_results
    return;
end
max_bytes = local_get_field(xcorr_params, 'max_reusable_mat_bytes', 25 * 1024^2);
info = dir(xcorr_file);
if isempty(info)
    return;
end
tf = info.bytes > max_bytes;
end


function key = local_density_source_names(sources)
if isempty(sources)
    key = "";
    return;
end
if ~iscell(sources)
    sources = num2cell(sources(:));
end
names = strings(1, numel(sources));
for i = 1:numel(sources)
    names(i) = string(local_get_field(sources{i}, 'name', sprintf('density%d', i)));
end
key = strjoin(sort(names), "|");
end


function value = local_join_group_field(act_result, field_name)
value = '';
if ~isstruct(act_result) || ~isfield(act_result, 'group_results') || isempty(act_result.group_results)
    return;
end
raw = {act_result.group_results.(field_name)};
raw = raw(~cellfun(@isempty, raw));
if isempty(raw)
    return;
end
value = strjoin(cellstr(string(raw(:)).'), ';');
end


function n = local_count_group_results(act_result)
n = 0;
if isstruct(act_result) && isfield(act_result, 'group_results') && ~isempty(act_result.group_results)
    n = numel(act_result.group_results);
end
end


function status = local_combine_status(a, b, enabled)
if ~enabled
    status = a;
elseif strcmp(a, 'error') || strcmp(b, 'error')
    status = 'error';
elseif strcmp(a, 'ok') || strcmp(b, 'ok')
    status = 'ok';
elseif strcmp(a, 'skipped_existing') && strcmp(b, 'skipped_existing')
    status = 'skipped_existing';
else
    status = sprintf('%s+%s', a, b);
end
end


function message = local_compose_message(xcorr_status, activation_status, act_result, roi_status, roi_result)
parts = {sprintf('xcorr=%s', xcorr_status), sprintf('activation=%s', activation_status), ...
    sprintf('roi=%s', roi_status)};
if isstruct(act_result) && isfield(act_result, 'message') && ~isempty(act_result.message)
    parts{end + 1} = act_result.message;
end
if isstruct(roi_result) && isfield(roi_result, 'message') && ~isempty(roi_result.message)
    parts{end + 1} = roi_result.message;
end
message = strjoin(parts, ' | ');
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
