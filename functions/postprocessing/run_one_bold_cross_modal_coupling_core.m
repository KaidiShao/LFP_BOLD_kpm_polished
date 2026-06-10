function result = run_one_bold_cross_modal_coupling_core(candidate, params)
%RUN_ONE_BOLD_CROSS_MODAL_COUPLING_CORE Run pipeline 8 for one BOLD_POST artifact.

if nargin < 1 || isempty(candidate)
    error('candidate is required.');
end
if nargin < 2 || isempty(params)
    params = build_bold_cross_modal_coupling_params();
else
    params = local_apply_defaults(params);
end

t_run = tic;
candidate = local_normalize_candidate(candidate, params.processed_root);
density_sources = params.density_sources;
if isempty(density_sources) && params.use_default_density_triplet
    fprintf('[P8 core] Resolving density sources...\n');
    density_sources = resolve_bold_cross_modal_density_sources(candidate, params);
end
[density_source_names, n_density_sources] = local_density_source_summary(density_sources);
fprintf('[P8 core] Density sources resolved (%d): %s\n', ...
    n_density_sources, density_source_names);
xcorr_params = params.xcorr;
xcorr_params.save_dir = candidate.xcorr_dir;
xcorr_params.save_tag = local_get_field(xcorr_params, 'save_tag', 'xcorr');
candidate.xcorr_file = fullfile(candidate.xcorr_dir, [xcorr_params.save_tag, '.mat']);

xcorr_recomputed_reason = '';
xcorr_recomputed = false;
if params.skip_existing_xcorr && exist(candidate.xcorr_file, 'file') == 2 && ...
        local_saved_xcorr_matches_density_sources(candidate.xcorr_file, density_sources, xcorr_params)
    xcorr_status = 'skipped_existing';
    xcorr_input = candidate.xcorr_file;
    if params.load_existing_xcorr
        xcorr_out = local_load_saved_xcorr(candidate.xcorr_file);
    else
        xcorr_out = struct();
    end
else
    fprintf('[P8 core] Computing xcorr...\n');
    xcorr_out = compute_bold_efun_density_cross_correlation( ...
        candidate.bold_post_file, density_sources, xcorr_params);
    fprintf('[P8 core] Xcorr computed.\n');
    xcorr_input = xcorr_out;
    xcorr_status = 'ok';
    xcorr_recomputed = true;
    if params.skip_existing_xcorr && exist(candidate.xcorr_file, 'file') == 2
        xcorr_recomputed_reason = 'recomputed because saved xcorr density sources did not match the current request';
    end
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
    if ~isfield(activation_params, 'top_n') || isempty(activation_params.top_n)
        activation_params.top_n = xcorr_params.top_n;
    end
    fprintf('[P8 core] Exporting activation maps...\n');
    act_result = export_bold_top_xcorr_activation_maps( ...
        candidate.bold_post_file, xcorr_input, activation_params);
    fprintf('[P8 core] Activation map export finished: %s\n', act_result.status);
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
    if ~isfield(roi_summary_params, 'top_n') || isempty(roi_summary_params.top_n)
        roi_summary_params.top_n = xcorr_params.top_n;
    end
    fprintf('[P8 core] Exporting ROI summaries...\n');
    roi_summary_result = export_bold_top_xcorr_roi_bar_summaries( ...
        candidate.bold_post_file, xcorr_input, roi_summary_params);
    fprintf('[P8 core] ROI summary export finished: %s\n', roi_summary_result.status);
    roi_summary_status = roi_summary_result.status;
end

result = struct();
result.candidate = candidate;
result.dataset_stem = candidate.dataset_stem;
result.dataset_id = candidate.dataset_id;
result.run_name = candidate.run_name;
result.run_tag = local_get_field(candidate, 'run_tag', '');
result.observable_mode = candidate.observable_mode;
result.residual_form = candidate.residual_form;
result.bold_post_file = candidate.bold_post_file;
result.xcorr_dir = candidate.xcorr_dir;
result.xcorr_mat_file = candidate.xcorr_file;
result.xcorr_peak_csv_file = fullfile(candidate.xcorr_dir, ...
    [xcorr_params.save_tag, '_peaks.csv']);
result.xcorr_top_csv_file = fullfile(candidate.xcorr_dir, ...
    [xcorr_params.save_tag, '_top.csv']);
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
result.status = local_combine_status(xcorr_status, activation_status, params.activation.enabled);
result.status = local_combine_status(result.status, roi_summary_status, roi_summary_enabled);
result.xcorr_status = xcorr_status;
result.activation_status = activation_status;
result.roi_summary_status = roi_summary_status;
result.message = local_compose_message( ...
    xcorr_status, xcorr_recomputed_reason, activation_status, act_result, ...
    roi_summary_status, roi_summary_result);
result.n_requested_maps = local_get_field(act_result, 'n_requested_maps', 0);
result.n_saved_maps = local_get_field(act_result, 'n_saved_maps', 0);
result.n_selected_roi_summary_modes = local_get_field(roi_summary_result, 'n_selected_modes', 0);
result.n_saved_roi_summary_figures = local_get_field(roi_summary_result, 'n_saved_figures', 0);
result.runtime_sec = toc(t_run);
result.xcorr_params = xcorr_params;
result.activation_params = activation_params;
result.roi_summary_params = roi_summary_params;
result.xcorr_out = xcorr_out;
result.act_result = act_result;
result.roi_summary_result = roi_summary_result;
fprintf('[P8 core] Syncing summary figures...\n');
result.summary_sync = sync_one_pipeline8_cross_modal_summary_figures(candidate, result, params.processed_root);
result.summary_sync_status = local_get_field(result.summary_sync, 'status', '');
result.summary_dir = local_get_field(result.summary_sync, 'summary_dir', '');
result.n_summary_figures = local_get_field(result.summary_sync, 'n_copied_files', 0);
result.activation_summary_sync = sync_one_pipeline8_top_xcorr_activation_summary_figures(candidate, result, params.processed_root);
result.activation_summary_sync_status = local_get_field(result.activation_summary_sync, 'status', '');
result.activation_summary_dir = local_get_field(result.activation_summary_sync, 'summary_dir', '');
result.n_activation_summary_figures = local_get_field(result.activation_summary_sync, 'n_copied_files', 0);
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


function candidate = local_normalize_candidate(candidate, processed_root)
if ischar(candidate) || isstring(candidate)
    candidate = local_candidate_from_bold_post_file(char(string(candidate)), processed_root);
elseif isstruct(candidate)
    if ~isfield(candidate, 'bold_post_file') || isempty(candidate.bold_post_file)
        error('candidate.bold_post_file is required.');
    end
    template = local_candidate_from_bold_post_file(candidate.bold_post_file, processed_root);
    names = fieldnames(template);
    for i = 1:numel(names)
        name = names{i};
        if ~isfield(candidate, name) || isempty(candidate.(name))
            candidate.(name) = template.(name);
        end
    end
else
    error('candidate must be a BOLD_POST path or candidate struct.');
end
end


function candidate = local_candidate_from_bold_post_file(bold_post_file, processed_root)
[mat_dir, base_name, ~] = fileparts(bold_post_file);
run_root = fileparts(mat_dir);
post_root = fileparts(run_root);
dataset_dir = fileparts(post_root);
[~, dataset_stem] = fileparts(dataset_dir);
run_name = char(erase(string(base_name), "_bold_post"));
save_tag = 'xcorr';

candidate = struct();
candidate.dataset_stem = dataset_stem;
candidate.dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
candidate.run_name = run_name;
candidate.observable_mode = local_parse_observable_mode(run_name);
candidate.residual_form = local_parse_residual_form(run_name);
candidate.run_tag = make_bold_cross_modal_run_tag(candidate);
candidate.bold_post_file = bold_post_file;
candidate.xcorr_dir = fullfile( ...
    io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 8, ...
    'efun_density_cross_correlation'), ...
    candidate.run_tag);
candidate.xcorr_file = fullfile(candidate.xcorr_dir, [save_tag, '.mat']);
candidate.activation_root = fullfile( ...
    io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 8, ...
    'figures_bold_top_xcorr_activation_maps'), ...
    candidate.run_tag);
end


function out = local_load_saved_xcorr(xcorr_file)
S = load(xcorr_file, 'out');
if ~isfield(S, 'out')
    error('Saved xcorr file is missing variable ''out'': %s', xcorr_file);
end
out = S.out;
out.source_file = xcorr_file;
end


function status = local_combine_status(xcorr_status, activation_status, activation_enabled)
if activation_enabled && any(strcmp(activation_status, {'no_top_table', 'no_raw_indices'}))
    status = activation_status;
    return;
end
if strcmp(xcorr_status, 'ok') || strcmp(activation_status, 'ok')
    status = 'ok';
elseif strcmp(xcorr_status, 'skipped_existing') && ...
        (~activation_enabled || any(strcmp(activation_status, {'skipped_existing', 'disabled'})))
    status = 'skipped_existing';
else
    status = xcorr_status;
end
end


function message = local_compose_message( ...
        xcorr_status, xcorr_recomputed_reason, activation_status, act_result, ...
        roi_summary_status, roi_summary_result)
parts = {['xcorr=' xcorr_status]};
if ~isempty(xcorr_recomputed_reason)
    parts{end + 1} = xcorr_recomputed_reason; %#ok<AGROW>
end
if ~strcmp(activation_status, 'disabled')
    parts{end + 1} = ['activation=' activation_status]; %#ok<AGROW>
    act_message = local_get_field(act_result, 'message', '');
    if ~isempty(act_message)
        parts{end + 1} = act_message; %#ok<AGROW>
    end
end
if ~strcmp(roi_summary_status, 'disabled')
    parts{end + 1} = ['roi_summary=' roi_summary_status]; %#ok<AGROW>
    roi_message = local_get_field(roi_summary_result, 'message', '');
    if ~isempty(roi_message)
        parts{end + 1} = roi_message; %#ok<AGROW>
    end
end
message = strjoin(parts, ' | ');
end


function mode = local_parse_observable_mode(run_name)
known = {'global_slow_band_power_svd100', 'roi_mean_slow_band_power', ...
    'slow_band_power_svd', 'gsvd100_ds', 'global_svd100', 'HP_svd100', 'slow_band_power', ...
    'roi_mean', 'eleHP', 'HP', 'svd'};
mode = '';
for i = 1:numel(known)
    if contains(run_name, ['_' known{i}], 'IgnoreCase', true) || ...
            endsWith(run_name, known{i}, 'IgnoreCase', true)
        mode = known{i};
        return;
    end
end
end


function form = local_parse_residual_form(run_name)
forms = {'projected_vlambda', 'projected_kv'};
form = '';
for i = 1:numel(forms)
    if contains(run_name, forms{i}, 'IgnoreCase', true)
        form = forms{i};
        return;
    end
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function tf = local_saved_xcorr_matches_density_sources(xcorr_file, density_sources, xcorr_params)
tf = false;
if local_saved_xcorr_is_oversized_for_current_policy(xcorr_file, xcorr_params)
    return;
end
try
    S = load(xcorr_file, 'out');
    if ~isfield(S, 'out') || ~isfield(S.out, 'density_sources')
        return;
    end
    saved_signature = local_density_source_signature(S.out.density_sources);
    current_signature = local_density_source_signature(density_sources);
    tf = isequal(saved_signature, current_signature);
catch
    tf = false;
end
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


function [names_text, n_sources] = local_density_source_summary(density_sources)
names_text = '';
n_sources = 0;
if isempty(density_sources)
    return;
end
if iscell(density_sources)
    sources = density_sources;
else
    sources = num2cell(density_sources(:));
end
names = cell(numel(sources), 1);
for i = 1:numel(sources)
    source = sources{i};
    if isstruct(source) && isfield(source, 'name') && ~isempty(source.name)
        names{i} = char(string(source.name));
    else
        names{i} = sprintf('density%d', i);
    end
end
names_text = strjoin(names, '|');
n_sources = numel(names);
end


function signature = local_density_source_signature(density_sources)
if isempty(density_sources)
    signature = strings(0, 1);
    return;
end
if iscell(density_sources)
    sources = density_sources;
else
    sources = num2cell(density_sources(:));
end

signature = strings(numel(sources), 1);
for i = 1:numel(sources)
    source = sources{i};
    name = lower(string(local_get_field(source, 'name', sprintf('density%d', i))));
    type = lower(string(local_get_field(source, 'type', 'generic_density')));
    file = lower(string(local_get_field(source, 'file', '')));
    signature(i) = name + "|" + type + "|" + file;
end
signature = sort(signature);
end
