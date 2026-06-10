function result = export_bold_top_xcorr_activation_maps(bold_post_input, xcorr_input, params)
%EXPORT_BOLD_TOP_XCORR_ACTIVATION_MAPS Export xcorr-ranked BOLD mode maps.
%
% This pipeline 8 exporter intentionally reuses the pipeline 7 spatial
% plotting helpers so cross-modal activation maps follow the same
% observable-to-voxel mapping contract as intrinsic maps.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(xcorr_input)
    error('xcorr_input is required.');
end
if nargin < 3 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

t_run = tic;
[plot_ctx, B] = build_bold_activation_plot_context(bold_post_input, struct( ...
    'datapons_root', params.datapons_root, ...
    'roi_ts_file', params.roi_ts_file, ...
    'background_limits', params.background_limits, ...
    'feature_reduce', params.feature_reduce));
highlight_spec = build_bold_core_roi_highlight_spec(plot_ctx, params);
xcorr_out = load_bold_xcorr_output(xcorr_input);
run_info = plot_ctx.run_info;

group_specs = build_bold_top_xcorr_group_specs(xcorr_out, params);
if isempty(group_specs)
    result = local_make_result(run_info, 'no_top_table', ...
        'XCORR top tables are empty.', '', '', 0, 0, toc(t_run));
    return;
end

if isempty(params.output_root)
    out_root = fullfile( ...
        io_project.get_pipeline_stage_dir(params.processed_root, ...
        run_info.dataset_stem, 8, 'figures_bold_top_xcorr_activation_maps'), ...
        run_info.run_name);
else
    out_root = char(string(params.output_root));
end
mat_dir = fullfile(out_root, 'mat');
fig_dir = fullfile(out_root, 'fig');
if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end
group_results = repmat(local_empty_group_result(), numel(group_specs), 1);
primary_info_file = '';
primary_activation_dir = fig_dir;
total_requested = 0;
total_saved = 0;

for i_group = 1:numel(group_specs)
    group_spec = group_specs(i_group);
    act_dir = local_group_output_dir(fig_dir, 'activation_maps', params.top_n, group_spec);
    if exist(act_dir, 'dir') ~= 7, mkdir(act_dir); end

    [raw_indices, top_row_indices] = resolve_bold_top_xcorr_raw_indices( ...
        group_spec.top_table, B.EDMD_outputs, params.top_n);
    if isempty(raw_indices)
        group_results(i_group) = local_make_group_result(group_spec, ...
            'no_raw_indices', ...
            sprintf('%s top table did not map to any raw Koopman basis index.', ...
            group_spec.display_name), ...
            act_dir, '', 0, 0, 0, [], table(), repmat(local_empty_map_row(), 0, 1));
        continue;
    end

    selected_top_table = group_spec.top_table(top_row_indices, :);
    map_rows = repmat(local_empty_map_row(), numel(raw_indices), 1);
    n_saved = 0;
    n_existing = 0;
    for i_idx = 1:numel(raw_indices)
        raw_idx = raw_indices(i_idx);
        top_row = selected_top_table(i_idx, :);
        png_file = fullfile(act_dir, sprintf( ...
            '%s_top%02d_raw%04d_%s.png', ...
            group_spec.slug, i_idx, raw_idx, params.value_mode));
        fig_file = fullfile(act_dir, sprintf( ...
            '%s_top%02d_raw%04d_%s.fig', ...
            group_spec.slug, i_idx, raw_idx, params.value_mode));

        if params.skip_existing && exist(png_file, 'file') == 2 && ...
                (~params.save_fig || exist(fig_file, 'file') == 2)
            n_existing = n_existing + 1;
            map_rows(i_idx) = local_make_map_row(top_row, i_idx, raw_idx, ...
                png_file, fig_file, 'skipped_existing', struct(), group_spec);
            continue;
        end

        fig_params = struct();
        fig_params.value_mode = params.value_mode;
        fig_params.slice_list = params.slice_list;
        fig_params.tiles_per_row = params.tiles_per_row;
        fig_params.overlay_alpha = params.overlay_alpha;
        fig_params.feature_reduce = params.feature_reduce;
        fig_params.highlight_spec = highlight_spec;
        fig_params.annotate_regions = params.annotate_regions;
        fig_params.annotation_observable_modes = params.annotation_observable_modes;
        fig_params.annotation_exact_names = params.annotation_exact_names;
        fig_params.annotation_contains_tokens = params.annotation_contains_tokens;
        fig_params.annotation_repeat_each_slice = params.annotation_repeat_each_slice;
        fig_params.annotation_min_voxels_per_slice = params.annotation_min_voxels_per_slice;
        fig_params.annotation_text_color = params.annotation_text_color;
        fig_params.annotation_font_size = params.annotation_font_size;
        fig_params.annotation_font_weight = params.annotation_font_weight;
        fig_params.annotation_background_color = params.annotation_background_color;
        fig_params.annotation_margin = params.annotation_margin;
        fig_params.annotation_text_offset = params.annotation_text_offset;
        fig_params.annotation_min_text_spacing = params.annotation_min_text_spacing;
        fig_params.annotation_show_marker = params.annotation_show_marker;
        fig_params.outline_regions = params.outline_regions;
        fig_params.outline_observable_modes = params.outline_observable_modes;
        fig_params.outline_exact_names = params.outline_exact_names;
        fig_params.outline_contains_tokens = params.outline_contains_tokens;
        fig_params.outline_use_annotation_names = params.outline_use_annotation_names;
        fig_params.outline_line_width = params.outline_line_width;
        fig_params.outline_colors = params.outline_colors;
        fig_params.outline_show_legend = params.outline_show_legend;
        fig_params.title_text = local_compose_title(plot_ctx, top_row, raw_idx, ...
            i_idx, params, group_spec);
        [fig, plot_info] = plot_bold_activation_map_reference_style(plot_ctx, raw_idx, fig_params);

        if params.save_png
            local_export_png(fig, png_file, params.resolution);
        end
        if params.save_fig
            savefig(fig, fig_file);
        end
        close(fig);
        drawnow limitrate;

        n_saved = n_saved + 1;
        map_rows(i_idx) = local_make_map_row(top_row, i_idx, raw_idx, ...
            png_file, fig_file, 'ok', plot_info, group_spec);
    end

    info_dir = local_group_info_dir(mat_dir, group_spec);
    if exist(info_dir, 'dir') ~= 7, mkdir(info_dir); end
    info_file = fullfile(info_dir, sprintf( ...
        'act_info__%s.mat', group_spec.slug));
    activation_info = struct();
    activation_info.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    activation_info.run_info = run_info;
    activation_info.bold_post_file = local_get_field(B, 'source_file', '');
    activation_info.xcorr_file = local_get_field(xcorr_out, 'source_file', '');
    activation_info.observable_file = local_get_field(plot_ctx, 'observable_file', '');
    activation_info.roi_ts_file = plot_ctx.roi_ts_file;
    activation_info.target_region_names = plot_ctx.target_region_names;
    activation_info.region_labels = plot_ctx.region_labels;
    activation_info.n_var_region = plot_ctx.n_var_region;
    activation_info.volume_size = size(plot_ctx.ana);
    activation_info.source_label_origin = plot_ctx.source_label_origin;
    activation_info.source_space_kind = plot_ctx.source_space_kind;
    activation_info.feature_reduce = params.feature_reduce;
    activation_info.used_backprojection = ~plot_ctx.direct_voxel_mode;
    activation_info.mapping_summary = plot_ctx.mapping_summary;
    activation_info.value_mode = params.value_mode;
    activation_info.slice_list = params.slice_list;
    activation_info.tiles_per_row = params.tiles_per_row;
    activation_info.overlay_alpha = params.overlay_alpha;
    activation_info.top_n = params.top_n;
    activation_info.core_roi_highlight = highlight_spec.info;
    activation_info.annotate_regions = params.annotate_regions;
    activation_info.annotation_exact_names = cellstr(string(params.annotation_exact_names(:)).');
    activation_info.annotation_contains_tokens = cellstr(string(params.annotation_contains_tokens(:)).');
    activation_info.outline_regions = params.outline_regions;
    activation_info.outline_exact_names = cellstr(string(params.outline_exact_names(:)).');
    activation_info.outline_contains_tokens = cellstr(string(params.outline_contains_tokens(:)).');
    activation_info.outline_show_legend = params.outline_show_legend;
    activation_info.selection_scope = group_spec.scope;
    activation_info.selection_name = group_spec.display_name;
    activation_info.selection_slug = group_spec.slug;
    activation_info.selection_density_name = group_spec.density_name;
    activation_info.selection_feature_name = group_spec.feature_name;
    activation_info.raw_indices = raw_indices;
    activation_info.selected_top_table = selected_top_table;
    activation_info.map_rows = map_rows;
    activation_info.activation_dir = act_dir;
    save_mat_variable_atomic(info_file, 'activation_info', activation_info);

    if isempty(primary_info_file) || strcmp(group_spec.scope, 'combined')
        primary_info_file = info_file;
    end
    if strcmp(group_spec.scope, 'combined')
        primary_activation_dir = act_dir;
    end

    total_requested = total_requested + numel(raw_indices);
    total_saved = total_saved + n_saved;
    if n_saved > 0
        group_status = 'ok';
        group_message = sprintf('%s: saved %d map(s); %d already existed.', ...
            group_spec.display_name, n_saved, n_existing);
    else
        group_status = 'skipped_existing';
        group_message = sprintf('%s: all %d requested map(s) already existed.', ...
            group_spec.display_name, numel(raw_indices));
    end
    group_results(i_group) = local_make_group_result(group_spec, ...
        group_status, group_message, act_dir, info_file, ...
        numel(raw_indices), n_saved, n_existing, raw_indices, ...
        selected_top_table, map_rows);
end

status = local_combine_group_status(group_results);
message = local_compose_group_message(group_results);
result = local_make_result(run_info, status, message, fig_dir, primary_info_file, ...
    total_requested, total_saved, toc(t_run));
result.output_root = out_root;
result.primary_activation_dir = primary_activation_dir;
result.group_results = group_results;
result.activation_group_dirs = {group_results.activation_dir};
result.info_files = {group_results.info_file};
result.observable_file = local_get_field(plot_ctx, 'observable_file', '');
result.roi_ts_file = plot_ctx.roi_ts_file;
result.xcorr_file = local_get_field(xcorr_out, 'source_file', '');
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'processed_root', io_project.get_project_processed_root());
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'output_root', '');
params = local_set_default(params, 'roi_ts_file', '');
params = local_set_default(params, 'background_limits', []);
params = local_set_default(params, 'top_n', 5);
params = local_set_default(params, 'value_mode', 'abs');
params = local_set_default(params, 'slice_list', 1:20);
params = local_set_default(params, 'tiles_per_row', 10);
params = local_set_default(params, 'overlay_alpha', 0.86);
params = local_set_default(params, 'feature_reduce', 'mean');
params = local_set_default(params, 'skip_existing', true);
params = local_set_default(params, 'export_combined', true);
params = local_set_default(params, 'export_by_density', true);
params = local_set_default(params, 'export_by_density_feature', false);
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'resolution', 220);
params = local_set_default(params, 'highlight_core_rois', true);
params = local_set_default(params, 'core_roi_observable_modes', ...
    {'svd', 'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100', ...
    'roi_mean', 'roi_mean_slow_band_power'});
params = local_set_default(params, 'core_roi_exact_names', ...
    {'ParabrachialN', 'Brainstem', 'HP', 'Tha', 'LGN'});
params = local_set_default(params, 'core_roi_contains_tokens', {'V1'});
params = local_set_default(params, 'core_roi_outline_color', [1 1 1]);
params = local_set_default(params, 'core_roi_line_width', 1.6);
params = local_set_default(params, 'annotate_regions', false);
params = local_set_default(params, 'annotation_observable_modes', {});
params = local_set_default(params, 'annotation_exact_names', {});
params = local_set_default(params, 'annotation_contains_tokens', {});
params = local_set_default(params, 'annotation_repeat_each_slice', false);
params = local_set_default(params, 'annotation_min_voxels_per_slice', 5);
params = local_set_default(params, 'annotation_text_color', [1 1 1]);
params = local_set_default(params, 'annotation_font_size', 8);
params = local_set_default(params, 'annotation_font_weight', 'bold');
params = local_set_default(params, 'annotation_background_color', [0 0 0]);
params = local_set_default(params, 'annotation_margin', 1);
params = local_set_default(params, 'annotation_text_offset', [0 -2]);
params = local_set_default(params, 'annotation_min_text_spacing', 5);
params = local_set_default(params, 'annotation_show_marker', false);
params = local_set_default(params, 'outline_regions', false);
params = local_set_default(params, 'outline_observable_modes', {});
params = local_set_default(params, 'outline_exact_names', {});
params = local_set_default(params, 'outline_contains_tokens', {});
params = local_set_default(params, 'outline_use_annotation_names', true);
params = local_set_default(params, 'outline_line_width', 1.6);
params = local_set_default(params, 'outline_colors', []);
params = local_set_default(params, 'outline_show_legend', true);
end


function out = local_load_xcorr_output(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load(file, 'out');
    if ~isfield(S, 'out')
        error('XCORR variable ''out'' missing in %s.', file);
    end
    out = S.out;
    out.source_file = file;
elseif isstruct(input)
    if isfield(input, 'top_table')
        out = input;
    elseif isfield(input, 'out')
        out = input.out;
    else
        error('xcorr_input struct must be an xcorr result or contain field ''out''.');
    end
    if ~isfield(out, 'source_file')
        out.source_file = '';
    end
else
    error('xcorr_input must be a path or struct.');
end
end


function group_specs = local_build_selection_groups(xcorr_out, params)
group_specs = repmat(local_empty_group_spec(), 0, 1);
if params.export_combined && isfield(xcorr_out, 'top_table') && ~isempty(xcorr_out.top_table)
    group_spec = local_empty_group_spec();
    group_spec.scope = 'combined';
    group_spec.display_name = 'combined';
    group_spec.slug = 'combined';
    group_spec.density_name = '';
    group_spec.top_table = xcorr_out.top_table;
    group_specs(end + 1, 1) = group_spec; %#ok<AGROW>
end

if ~params.export_by_density
    return;
end

[density_names, field_names, top_tables] = local_density_group_top_tables(xcorr_out, params.top_n);
for i_group = 1:numel(density_names)
    top_table_i = top_tables{i_group};
    if isempty(top_table_i)
        continue;
    end
    group_spec = local_empty_group_spec();
    group_spec.scope = 'density';
    group_spec.display_name = density_names{i_group};
    group_spec.slug = field_names{i_group};
    group_spec.density_name = density_names{i_group};
    group_spec.top_table = top_table_i;
    group_specs(end + 1, 1) = group_spec; %#ok<AGROW>
end
end


function [density_names, field_names, top_tables] = local_density_group_top_tables(xcorr_out, top_n)
density_names = {};
field_names = {};
top_tables = {};
if isfield(xcorr_out, 'top_table_by_density') && isstruct(xcorr_out.top_table_by_density) && ...
        ~isempty(fieldnames(xcorr_out.top_table_by_density))
    top_struct = xcorr_out.top_table_by_density;
    if isfield(xcorr_out, 'density_group_names') && isfield(xcorr_out, 'density_group_fields') && ...
            numel(xcorr_out.density_group_names) == numel(xcorr_out.density_group_fields)
        density_names = cellstr(string(xcorr_out.density_group_names(:)).');
        field_names = cellstr(string(xcorr_out.density_group_fields(:)).');
    else
        field_names = fieldnames(top_struct).';
        density_names = field_names;
    end
    top_tables = cell(size(field_names));
    for i_group = 1:numel(field_names)
        top_tables{i_group} = top_struct.(field_names{i_group});
    end
    return;
end

if ~isfield(xcorr_out, 'peak_table') || isempty(xcorr_out.peak_table)
    return;
end

[top_struct, density_names, field_names] = local_build_density_group_top_struct( ...
    xcorr_out.peak_table, top_n);
top_tables = cell(size(field_names));
for i_group = 1:numel(field_names)
    top_tables{i_group} = top_struct.(field_names{i_group});
end
end


function [top_struct, density_names, field_names] = ...
        local_build_density_group_top_struct(peak_table, top_n)
top_struct = struct();
density_names = {};
field_names = {};
if isempty(peak_table) || ~ismember('density_name', peak_table.Properties.VariableNames)
    return;
end
all_names = cellstr(string(peak_table.density_name));
[density_names, ~] = unique(all_names, 'stable');
density_names = density_names(:).';
field_names = cell(size(density_names));
for i_group = 1:numel(density_names)
    density_name = density_names{i_group};
    field_name = local_density_group_field_name(density_name);
    mask = strcmp(all_names, density_name);
    group_table = local_sort_peak_table(peak_table(mask, :));
    if isempty(group_table)
        top_struct.(field_name) = group_table;
    else
        top_struct.(field_name) = group_table(1:min(height(group_table), top_n), :);
    end
    field_names{i_group} = field_name;
end
end


function [raw_indices, top_row_indices] = local_resolve_top_raw_indices(top_table, EDMD_outputs, top_n)
raw_indices = nan(1, 0);
top_row_indices = nan(1, 0);
if isempty(top_table)
    raw_indices = raw_indices(:).';
    top_row_indices = top_row_indices(:).';
    return;
end
n = height(top_table);
for ii = 1:n
    if numel(raw_indices) >= top_n
        break;
    end
    mode_idx = top_table.bold_mode_index(ii);
    raw_idx = local_mode_to_raw_index(mode_idx, EDMD_outputs);
    if ~isfinite(raw_idx) || raw_idx < 1
        continue;
    end
    if raw_idx > size(EDMD_outputs.kpm_modes, 1)
        continue;
    end
    if any(raw_indices == raw_idx)
        continue;
    end
    raw_indices(end + 1) = round(raw_idx); %#ok<AGROW>
    top_row_indices(end + 1) = ii; %#ok<AGROW>
end
raw_indices = raw_indices(:).';
top_row_indices = top_row_indices(:).';
end


function raw_idx = local_mode_to_raw_index(mode_idx, EDMD_outputs)
raw_idx = double(mode_idx);
if isfield(EDMD_outputs, 'idx_final_in_original') && ...
        numel(EDMD_outputs.idx_final_in_original) >= mode_idx && ...
        ~isempty(EDMD_outputs.idx_final_in_original(mode_idx))
    raw_idx = double(EDMD_outputs.idx_final_in_original(mode_idx));
end
end


function title_text = local_compose_title(plot_ctx, top_row, raw_idx, top_rank, params, group_spec)
run_name = local_get_field(plot_ctx.run_info, 'run_name', 'BOLD');
density_name = local_table_char(top_row, 'density_name', 'density');
density_label = local_table_char(top_row, 'density_label', '');
bold_feature = local_table_char(top_row, 'bold_feature', 'feature');
mode_idx = local_table_scalar(top_row, 'bold_mode_index', raw_idx);
peak_corr = local_table_scalar(top_row, 'peak_corr', NaN);
peak_lag_sec = local_table_scalar(top_row, 'peak_lag_sec', NaN);

if isempty(density_label) || strcmpi(density_label, density_name)
    density_text = density_name;
else
    density_text = sprintf('%s:%s', density_name, density_label);
end

title_text = sprintf('%s | %s | top%02d | %s | %s | mode %d raw %d | corr %.3f lag %.1fs | xcorr %s', ...
    run_name, group_spec.display_name, top_rank, density_text, bold_feature, ...
    mode_idx, raw_idx, peak_corr, peak_lag_sec, params.value_mode);
end


function out_dir = local_group_output_dir(fig_dir, prefix, top_n, group_spec)
if contains(lower(char(string(prefix))), 'activation')
    dir_name = 'act';
else
    dir_name = sprintf('%s_top%d', prefix, top_n);
end
family = local_group_feature_family(group_spec);
if isempty(family)
    out_dir = fullfile(fig_dir, dir_name);
else
    out_dir = fullfile(fig_dir, family, dir_name);
end
end


function info_dir = local_group_info_dir(mat_dir, group_spec)
family = local_group_feature_family(group_spec);
if isempty(family)
    info_dir = mat_dir;
else
    info_dir = fullfile(mat_dir, family);
end
end


function family = local_group_feature_family(group_spec)
family = '';
feature_name = local_get_field(group_spec, 'feature_name', '');
if isempty(feature_name)
    return;
end
if contains(lower(char(string(feature_name))), 'deconv')
    family = 'deconv_efun';
else
    family = 'efun';
end
end


function group = local_empty_group_spec()
group = struct( ...
    'scope', '', ...
    'display_name', '', ...
    'slug', '', ...
    'density_name', '', ...
    'feature_name', '', ...
    'top_table', table());
end


function group_result = local_empty_group_result()
group_result = struct( ...
    'scope', '', ...
    'display_name', '', ...
    'slug', '', ...
    'density_name', '', ...
    'feature_name', '', ...
    'status', '', ...
    'message', '', ...
    'activation_dir', '', ...
    'info_file', '', ...
    'n_requested_maps', 0, ...
    'n_saved_maps', 0, ...
    'n_existing_maps', 0, ...
    'raw_indices', [], ...
    'selected_top_table', table(), ...
    'map_rows', repmat(local_empty_map_row(), 0, 1));
end


function group_result = local_make_group_result(group_spec, status, message, ...
        activation_dir, info_file, n_requested_maps, n_saved_maps, ...
        n_existing_maps, raw_indices, selected_top_table, map_rows)
group_result = local_empty_group_result();
group_result.scope = group_spec.scope;
group_result.display_name = group_spec.display_name;
group_result.slug = group_spec.slug;
group_result.density_name = group_spec.density_name;
group_result.feature_name = group_spec.feature_name;
group_result.status = status;
group_result.message = message;
group_result.activation_dir = activation_dir;
group_result.info_file = info_file;
group_result.n_requested_maps = n_requested_maps;
group_result.n_saved_maps = n_saved_maps;
group_result.n_existing_maps = n_existing_maps;
group_result.raw_indices = raw_indices;
group_result.selected_top_table = selected_top_table;
group_result.map_rows = map_rows;
end


function status = local_combine_group_status(group_results)
if isempty(group_results)
    status = 'no_top_table';
    return;
end
statuses = {group_results.status};
if any(strcmp(statuses, 'ok'))
    status = 'ok';
elseif all(strcmp(statuses, 'skipped_existing'))
    status = 'skipped_existing';
elseif all(strcmp(statuses, 'no_raw_indices'))
    status = 'no_raw_indices';
elseif any(strcmp(statuses, 'skipped_existing'))
    status = 'skipped_existing';
else
    status = statuses{1};
end
end


function message = local_compose_group_message(group_results)
if isempty(group_results)
    message = 'No activation-map groups were available.';
    return;
end
parts = cell(1, numel(group_results));
for i_group = 1:numel(group_results)
    parts{i_group} = group_results(i_group).message;
end
message = strjoin(parts, ' | ');
end


function row = local_empty_map_row()
row = struct( ...
    'top_rank', NaN, ...
    'raw_index', NaN, ...
    'mode_index', NaN, ...
    'selection_scope', '', ...
    'selection_name', '', ...
    'selection_slug', '', ...
    'density_name', '', ...
    'selection_feature_name', '', ...
    'density_label', '', ...
    'bold_feature', '', ...
    'peak_corr', NaN, ...
    'peak_lag_sec', NaN, ...
    'png_file', '', ...
    'fig_file', '', ...
    'status', '', ...
    'map_space', '', ...
    'slice_list', []);
end


function row = local_make_map_row(top_row, top_rank, raw_idx, png_file, fig_file, ...
        status, plot_info, group_spec)
row = local_empty_map_row();
row.top_rank = top_rank;
row.raw_index = raw_idx;
row.mode_index = local_table_scalar(top_row, 'bold_mode_index', raw_idx);
row.selection_scope = group_spec.scope;
row.selection_name = group_spec.display_name;
row.selection_slug = group_spec.slug;
row.selection_feature_name = group_spec.feature_name;
row.density_name = local_table_char(top_row, 'density_name', '');
row.density_label = local_table_char(top_row, 'density_label', '');
row.bold_feature = local_table_char(top_row, 'bold_feature', '');
row.peak_corr = local_table_scalar(top_row, 'peak_corr', NaN);
row.peak_lag_sec = local_table_scalar(top_row, 'peak_lag_sec', NaN);
row.png_file = png_file;
row.fig_file = fig_file;
row.status = status;
if isstruct(plot_info)
    row.map_space = local_get_field(plot_info, 'map_space', '');
    row.slice_list = local_get_field(plot_info, 'slice_list', []);
end
end


function T = local_sort_peak_table(T)
if isempty(T)
    return;
end
finite_peak = isfinite(T.peak_abs_corr);
T = [ ...
    sortrows(T(finite_peak, :), 'peak_abs_corr', 'descend'); ...
    T(~finite_peak, :)];
end


function field_name = local_density_group_field_name(name)
field_name = char(matlab.lang.makeValidName(char(string(name))));
if isempty(field_name)
    field_name = 'density_group';
end
end


function result = local_make_result(run_info, status, message, activation_dir, ...
        info_file, n_requested_maps, n_saved_maps, runtime_sec)
result = struct();
result.dataset_stem = run_info.dataset_stem;
result.run_name = run_info.run_name;
result.observable_mode = local_get_field(run_info, 'observable_mode', '');
result.residual_form = local_get_field(run_info, 'residual_form', '');
result.status = status;
result.message = message;
result.activation_dir = activation_dir;
result.info_file = info_file;
result.n_requested_maps = n_requested_maps;
result.n_saved_maps = n_saved_maps;
result.runtime_sec = runtime_sec;
end


function value = local_table_char(T, name, default_value)
value = default_value;
if ~istable(T) || ~ismember(name, T.Properties.VariableNames) || isempty(T.(name))
    return;
end
raw = T.(name);
if iscell(raw)
    raw = raw{1};
elseif isstring(raw)
    raw = raw(1);
elseif numel(raw) > 1
    raw = raw(1);
end
value = char(string(raw));
end


function value = local_table_scalar(T, name, default_value)
value = default_value;
if ~istable(T) || ~ismember(name, T.Properties.VariableNames) || isempty(T.(name))
    return;
end
raw = T.(name);
if isnumeric(raw) || islogical(raw)
    value = double(raw(1));
elseif iscell(raw) && ~isempty(raw{1}) && isnumeric(raw{1})
    value = double(raw{1});
else
    value = default_value;
end
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function local_export_png(fig, png_file, resolution)
[png_dir, ~, ~] = fileparts(png_file);
if exist(png_dir, 'dir') ~= 7
    mkdir(png_dir);
end

tmp_dir = fullfile(tempdir, 'koopman_png_shortpath');
if exist(tmp_dir, 'dir') ~= 7
    mkdir(tmp_dir);
end
tmp_png = [tempname(tmp_dir), '.png'];
cleanup_tmp = onCleanup(@() local_delete_if_exists(tmp_png));

exportgraphics(fig, tmp_png, 'Resolution', resolution);
[ok, msg] = copyfile(tmp_png, png_file, 'f');
if ~ok
    error('export_bold_top_xcorr_activation_maps:CopyPngFailed', ...
        'Unable to copy exported PNG to target path:\n%s\n%s', png_file, msg);
end

clear cleanup_tmp
end


function local_delete_if_exists(path_in)
if exist(path_in, 'file') == 2
    try
        delete(path_in);
    catch
    end
end
end
