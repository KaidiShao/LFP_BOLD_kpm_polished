function [result, B] = export_bold_intrinsic_activation_maps(bold_post_input, params)
%EXPORT_BOLD_INTRINSIC_ACTIVATION_MAPS Export intrinsic BOLD mode activation maps for pipeline 7.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

t_run = tic;
[plot_ctx, B] = build_bold_activation_plot_context(bold_post_input, params);
run_info = plot_ctx.run_info;

if isempty(params.output_root)
    out_root = fullfile( ...
        io_project.get_pipeline_stage_dir(params.processed_root, ...
        run_info.dataset_stem, 7, 'bold_postprocessing'), ...
        run_info.run_name);
else
    out_root = char(string(params.output_root));
end
mat_dir = fullfile(out_root, 'mat');
fig_dir = fullfile(out_root, 'fig');
act_dir = fullfile(fig_dir, 'intrinsic_activation_maps');
if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end
if exist(act_dir, 'dir') ~= 7, mkdir(act_dir); end

[raw_indices, selection_labels] = local_resolve_raw_indices( ...
    plot_ctx.evalues, params.selection_mode, params.basis_indices);
if isempty(raw_indices)
    result = local_make_result(run_info, 'no_indices', ...
        'No intrinsic basis indices were resolved.', act_dir, '', 0, 0, toc(t_run));
    return;
end

n_saved = 0;
n_existing = 0;
for i_idx = 1:numel(raw_indices)
    raw_idx = raw_indices(i_idx);
    selection_label = selection_labels{i_idx};
    stem = sprintf('intrinsic_%s_raw%04d_%s', ...
        selection_label, raw_idx, params.value_mode);
    png_file = fullfile(act_dir, [stem, '_activation_reference.png']);
    fig_file = fullfile(act_dir, [stem, '_activation_reference.fig']);

    if params.skip_existing && exist(png_file, 'file') == 2 && ...
            (~params.save_fig || exist(fig_file, 'file') == 2)
        n_existing = n_existing + 1;
        continue;
    end

    title_text = sprintf('%s | %s -> raw %d | |lambda|=%.4g | intrinsic %s', ...
        run_info.run_name, selection_label, raw_idx, abs(plot_ctx.evalues(raw_idx)), ...
        params.value_mode);
    fig_params = struct();
    fig_params.value_mode = params.value_mode;
    fig_params.slice_list = params.slice_list;
    fig_params.tiles_per_row = params.tiles_per_row;
    fig_params.overlay_alpha = params.overlay_alpha;
    fig_params.feature_reduce = params.feature_reduce;
    fig_params.title_text = title_text;
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
    fig = plot_bold_activation_map_reference_style(plot_ctx, raw_idx, fig_params);

    if params.save_png
        exportgraphics(fig, png_file, 'Resolution', params.resolution);
    end
    if params.save_fig
        savefig(fig, fig_file);
    end
    close(fig);
    drawnow limitrate;
    n_saved = n_saved + 1;
end

info_file = fullfile(mat_dir, 'intrinsic_activation_map_info.mat');
activation_info = struct();
activation_info.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
activation_info.run_info = run_info;
activation_info.bold_post_file = local_get_field(B, 'source_file', '');
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
activation_info.selection_mode = params.selection_mode;
activation_info.basis_indices = params.basis_indices(:).';
activation_info.selection_labels = selection_labels;
activation_info.raw_indices = raw_indices;
activation_info.value_mode = params.value_mode;
activation_info.slice_list = params.slice_list;
activation_info.tiles_per_row = params.tiles_per_row;
activation_info.annotate_regions = params.annotate_regions;
activation_info.annotation_exact_names = cellstr(string(params.annotation_exact_names(:)).');
activation_info.annotation_contains_tokens = cellstr(string(params.annotation_contains_tokens(:)).');
activation_info.outline_regions = params.outline_regions;
activation_info.outline_exact_names = cellstr(string(params.outline_exact_names(:)).');
activation_info.outline_contains_tokens = cellstr(string(params.outline_contains_tokens(:)).');
activation_info.outline_show_legend = params.outline_show_legend;
activation_info.activation_dir = act_dir;
save(info_file, 'activation_info', '-v7.3');

B.activation_context = struct();
B.activation_context.roi_ts_file = plot_ctx.roi_ts_file;
B.activation_context.target_region_names = plot_ctx.target_region_names;
B.activation_context.region_labels = plot_ctx.region_labels;
B.activation_context.n_var_region = plot_ctx.n_var_region;
B.activation_context.volume_size = size(plot_ctx.ana);
B.activation_context.source_label_origin = plot_ctx.source_label_origin;
B.activation_context.source_space_kind = plot_ctx.source_space_kind;
B.activation_context.feature_reduce = params.feature_reduce;
B.activation_context.used_backprojection = ~plot_ctx.direct_voxel_mode;
B.activation_context.mapping_summary = plot_ctx.mapping_summary;
B.artifacts.intrinsic_activation_dir = act_dir;
B.artifacts.intrinsic_activation_info_file = info_file;

if params.update_bold_post && isfield(B, 'source_file') && ~isempty(B.source_file)
    BOLD_POST = B;
    save_bold_post_mat(B.source_file, BOLD_POST);
end

if n_saved > 0
    status = 'ok';
    message = sprintf('Saved %d intrinsic map(s); %d already existed.', n_saved, n_existing);
else
    status = 'skipped_existing';
    message = sprintf('All %d intrinsic map(s) already existed.', numel(raw_indices));
end

result = local_make_result(run_info, status, message, act_dir, info_file, ...
    numel(raw_indices), n_saved, toc(t_run));
result.raw_indices = raw_indices;
result.selection_labels = selection_labels;
result.roi_ts_file = plot_ctx.roi_ts_file;
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'processed_root', io_project.get_project_processed_root());
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'output_root', '');
params = local_set_default(params, 'selection_mode', 'sorted');
params = local_set_default(params, 'basis_indices', 1:5);
params = local_set_default(params, 'value_mode', 'abs');
params = local_set_default(params, 'slice_list', 1:20);
params = local_set_default(params, 'tiles_per_row', 10);
params = local_set_default(params, 'overlay_alpha', 0.86);
params = local_set_default(params, 'feature_reduce', 'mean');
params = local_set_default(params, 'skip_existing', true);
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'resolution', 220);
params = local_set_default(params, 'update_bold_post', true);
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


function [raw_indices, selection_labels] = local_resolve_raw_indices(evalues, selection_mode, basis_indices)
basis_indices = unique(round(double(basis_indices(:).')), 'stable');
basis_indices = basis_indices(isfinite(basis_indices) & basis_indices >= 1);
selection_labels = {};
if isempty(basis_indices)
    raw_indices = [];
    return;
end

switch lower(char(string(selection_mode)))
    case 'sorted'
        [~, order] = sort(abs(evalues), 'descend');
        valid = basis_indices(basis_indices <= numel(order));
        raw_indices = order(valid);
        selection_labels = arrayfun(@(x) sprintf('sorted%03d', x), valid, 'UniformOutput', false);
    case 'raw'
        valid = basis_indices(basis_indices <= numel(evalues));
        raw_indices = valid;
        selection_labels = arrayfun(@(x) sprintf('raw%03d', x), valid, 'UniformOutput', false);
    otherwise
        error('selection_mode must be ''sorted'' or ''raw''.');
end

raw_indices = raw_indices(:).';
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
