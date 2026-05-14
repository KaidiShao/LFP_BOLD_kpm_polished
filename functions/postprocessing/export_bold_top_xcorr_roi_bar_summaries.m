function result = export_bold_top_xcorr_roi_bar_summaries(bold_post_input, xcorr_input, params)
%EXPORT_BOLD_TOP_XCORR_ROI_BAR_SUMMARIES Export pipeline 8 top-xcorr ROI bar summaries.

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
    'feature_reduce', params.feature_reduce));
xcorr_out = load_bold_xcorr_output(xcorr_input);
run_info = plot_ctx.run_info;
group_specs = build_bold_top_xcorr_group_specs(xcorr_out, params);
if isempty(group_specs)
    result = local_make_result(run_info, 'no_top_table', ...
        'XCORR top tables are empty.', '', '', 0, 0, toc(t_run));
    result.group_results = repmat(local_empty_group_result(), 0, 1);
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
primary_summary_dir = '';
total_selected_modes = 0;
total_saved_figures = 0;

for i_group = 1:numel(group_specs)
    group_spec = group_specs(i_group);
    summary_dir = local_group_output_dir(fig_dir, 'roi_bar_summaries', params.top_n, group_spec);
    if exist(summary_dir, 'dir') ~= 7, mkdir(summary_dir); end

    [raw_indices, top_row_indices] = resolve_bold_top_xcorr_raw_indices( ...
        group_spec.top_table, B.EDMD_outputs, params.top_n);
    if isempty(raw_indices)
        group_results(i_group) = local_make_group_result(group_spec, ...
            'no_raw_indices', ...
            sprintf('%s top table did not map to any raw Koopman basis index.', ...
            group_spec.display_name), ...
            summary_dir, '', '', 0, 0, [], table());
        continue;
    end

    selected_top_table = group_spec.top_table(top_row_indices, :);
    stem = local_build_summary_stem(group_spec, params);
    png_file = fullfile(summary_dir, [stem, '.png']);
    fig_file = fullfile(summary_dir, [stem, '.fig']);
    info_dir = local_group_info_dir(mat_dir, group_spec);
    if exist(info_dir, 'dir') ~= 7, mkdir(info_dir); end
    info_file = fullfile(info_dir, [stem, '_info.mat']);

    plot_out = struct();
    if params.skip_existing && exist(png_file, 'file') == 2 && ...
            (~params.save_fig || exist(fig_file, 'file') == 2) && ...
            exist(info_file, 'file') == 2
        plot_out = local_load_existing_plot_out(info_file);
        n_saved = 0;
        group_status = 'skipped_existing';
        group_message = sprintf('%s ROI bar summary already existed.', group_spec.display_name);
    else
        plot_params = params;
        plot_params.selection_mode = 'raw';
        plot_params.basis_indices = raw_indices;
        plot_params.mode_title_texts = local_build_mode_title_texts(selected_top_table, raw_indices, group_spec);
        plot_params.figure_title_text = local_build_figure_title_text(run_info, group_spec, plot_ctx, params);
        [fig, plot_out] = plot_bold_mode_roi_bar_summary(plot_ctx, plot_params);
        if params.save_png
            local_export_png(fig, png_file, params.resolution);
        end
        if params.save_fig
            savefig(fig, fig_file);
        end
        close(fig);
        drawnow limitrate;
        n_saved = 1;
        group_status = 'ok';
        group_message = sprintf('%s ROI bar summary saved with %d mode(s).', ...
            group_spec.display_name, numel(raw_indices));
    end

    summary_info = struct();
    summary_info.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    summary_info.run_info = run_info;
    summary_info.bold_post_file = local_get_field(B, 'source_file', '');
    summary_info.xcorr_file = local_get_field(xcorr_out, 'source_file', '');
    summary_info.observable_file = local_get_field(plot_ctx, 'observable_file', '');
    summary_info.roi_ts_file = plot_ctx.roi_ts_file;
    summary_info.selection_scope = group_spec.scope;
    summary_info.selection_name = group_spec.display_name;
    summary_info.selection_slug = group_spec.slug;
    summary_info.selection_density_name = group_spec.density_name;
    summary_info.selection_feature_name = group_spec.feature_name;
    summary_info.selected_top_table = selected_top_table;
    summary_info.raw_indices = raw_indices;
    summary_info.png_file = png_file;
    summary_info.fig_file = fig_file;
    summary_info.summary_dir = summary_dir;
    summary_info.params = params;
    summary_info.plot_out = plot_out;
    save_mat_variable_atomic(info_file, 'summary_info', summary_info);

    if isempty(primary_info_file) || strcmp(group_spec.scope, 'combined')
        primary_info_file = info_file;
        primary_summary_dir = summary_dir;
    end
    total_selected_modes = total_selected_modes + numel(raw_indices);
    total_saved_figures = total_saved_figures + n_saved;
    group_results(i_group) = local_make_group_result(group_spec, group_status, ...
        group_message, summary_dir, png_file, info_file, n_saved, ...
        numel(raw_indices), raw_indices, selected_top_table);
end

status = local_combine_group_status(group_results);
message = local_compose_group_message(group_results);
result = local_make_result(run_info, status, message, primary_summary_dir, ...
    primary_info_file, total_selected_modes, total_saved_figures, toc(t_run));
result.output_root = out_root;
result.group_results = group_results;
result.summary_group_dirs = {group_results.summary_dir};
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
params = local_set_default(params, 'top_n', 5);
params = local_set_default(params, 'feature_reduce', 'mean');
params = local_set_default(params, 'roi_reduce', 'mean');
params = local_set_default(params, 'roi_value_mode', 'mean_abs');
params = local_set_default(params, 'mode_normalization', 'range');
params = local_set_default(params, 'flip_roi_order', true);
params = local_set_default(params, 'layout_mode', 'row');
params = local_set_default(params, 'skip_existing', true);
params = local_set_default(params, 'export_combined', true);
params = local_set_default(params, 'export_by_density', true);
params = local_set_default(params, 'export_by_density_feature', false);
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'resolution', 220);
end


function texts = local_build_mode_title_texts(selected_top_table, raw_indices, group_spec)
n_mode = numel(raw_indices);
texts = cell(1, n_mode);
for i_mode = 1:n_mode
    top_row = selected_top_table(i_mode, :);
    density_name = local_table_char(top_row, 'density_name', group_spec.density_name);
    density_label = local_table_char(top_row, 'density_label', '');
    bold_feature = local_table_char(top_row, 'bold_feature', 'feature');
    mode_idx = local_table_scalar(top_row, 'bold_mode_index', raw_indices(i_mode));
    peak_corr = local_table_scalar(top_row, 'peak_corr', NaN);
    peak_lag_sec = local_table_scalar(top_row, 'peak_lag_sec', NaN);
    density_text = local_short_density_text(density_name, density_label);
    feature_text = local_short_feature_text(bold_feature);
    texts{i_mode} = sprintf('top%02d | %s | %s | raw %d\\newlinecorr %.3f | lag %.1fs', ...
        i_mode, density_text, feature_text, raw_indices(i_mode), peak_corr, peak_lag_sec);
end
end


function title_text = local_build_figure_title_text(run_info, group_spec, plot_ctx, params)
title_text = sprintf('%s | %s | ROI-average top-xcorr summary | %s | ROI order %s', ...
    local_compact_run_label(run_info), group_spec.display_name, ...
    local_describe_plot_space(plot_ctx), ...
    ternary(local_get_field(params, 'flip_roi_order', true), 'flipped', 'native'));
end


function out_dir = local_group_output_dir(fig_dir, prefix, top_n, group_spec)
if contains(lower(char(string(prefix))), 'roi')
    dir_name = 'roi';
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


function desc = local_describe_plot_space(plot_ctx)
kind = local_get_field(plot_ctx, 'source_space_kind', 'unknown');
reduce = local_get_field(plot_ctx, 'feature_reduce', 'mean');
used_backprojection = ~plot_ctx.direct_voxel_mode;
switch kind
    case 'voxel'
        if used_backprojection
            desc = 'SVD back-projected to voxel observables';
        else
            desc = 'direct voxel observables';
        end
    case 'region_mean'
        if used_backprojection
            desc = 'SVD back-projected to region means';
        else
            desc = 'region means';
        end
    case 'multi_feature_voxel'
        if used_backprojection
            desc = sprintf('SVD back-projected to multi-feature voxels (%s by voxel)', reduce);
        else
            desc = sprintf('multi-feature voxel observables (%s by voxel)', reduce);
        end
    otherwise
        desc = sprintf('%s source space', kind);
end
end


function stem = local_build_summary_stem(group_spec, params)
stem = sprintf('roi_%s_top%d_%s_%s_%s', ...
    group_spec.slug, params.top_n, char(string(params.roi_value_mode)), ...
    char(string(params.mode_normalization)), char(string(params.layout_mode)));
end


function plot_out = local_load_existing_plot_out(info_file)
plot_out = struct();
try
    S = load(info_file, 'summary_info');
    if isfield(S, 'summary_info') && isfield(S.summary_info, 'plot_out')
        plot_out = S.summary_info.plot_out;
    end
catch
    plot_out = struct();
end
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
    'summary_dir', '', ...
    'png_file', '', ...
    'info_file', '', ...
    'n_saved_figures', 0, ...
    'n_selected_modes', 0, ...
    'raw_indices', [], ...
    'selected_top_table', table());
end


function group_result = local_make_group_result(group_spec, status, message, ...
        summary_dir, png_file, info_file, n_saved_figures, ...
        n_selected_modes, raw_indices, selected_top_table)
group_result = local_empty_group_result();
group_result.scope = group_spec.scope;
group_result.display_name = group_spec.display_name;
group_result.slug = group_spec.slug;
group_result.density_name = group_spec.density_name;
group_result.feature_name = group_spec.feature_name;
group_result.status = status;
group_result.message = message;
group_result.summary_dir = summary_dir;
group_result.png_file = png_file;
group_result.info_file = info_file;
group_result.n_saved_figures = n_saved_figures;
group_result.n_selected_modes = n_selected_modes;
group_result.raw_indices = raw_indices;
group_result.selected_top_table = selected_top_table;
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
    message = 'No ROI bar-summary groups were available.';
    return;
end
parts = cell(1, numel(group_results));
for i_group = 1:numel(group_results)
    parts{i_group} = group_results(i_group).message;
end
message = strjoin(parts, ' | ');
end


function result = local_make_result(run_info, status, message, summary_dir, ...
        info_file, n_selected_modes, n_saved_figures, runtime_sec)
result = struct();
result.dataset_stem = run_info.dataset_stem;
result.run_name = run_info.run_name;
result.observable_mode = local_get_field(run_info, 'observable_mode', '');
result.residual_form = local_get_field(run_info, 'residual_form', '');
result.status = status;
result.message = message;
result.summary_dir = summary_dir;
result.info_file = info_file;
result.n_selected_modes = n_selected_modes;
result.n_saved_figures = n_saved_figures;
result.runtime_sec = runtime_sec;
end


function label = local_compact_run_label(run_info)
dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', '')));
observable_mode = char(string(local_get_field(run_info, 'observable_mode', '')));
residual_form = char(string(local_get_field(run_info, 'residual_form', '')));
run_name = char(string(local_get_field(run_info, 'run_name', '')));
parts = {};
if ~isempty(dataset_stem), parts{end + 1} = dataset_stem; end %#ok<AGROW>
if ~isempty(observable_mode), parts{end + 1} = observable_mode; end %#ok<AGROW>
if ~isempty(residual_form), parts{end + 1} = residual_form; end %#ok<AGROW>
if isempty(parts)
    label = run_name;
else
    label = strjoin(parts, ' | ');
end
end


function out = ternary(cond, a, b)
if cond
    out = a;
else
    out = b;
end
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


function text = local_short_density_text(density_name, density_label)
if ~isempty(density_label)
    text = density_label;
    return;
end
text = density_name;
text = strrep(text, 'blp_event_density', 'event');
text = strrep(text, 'blp_raw_eigenfunction_density', 'raw');
text = strrep(text, 'blp_dimred_eigenfunction_density', 'dimred');
text = strrep(text, '_', ' ');
end


function text = local_short_feature_text(bold_feature)
text = char(string(bold_feature));
switch lower(text)
    case 'efun_abs'
        text = 'efun|abs|';
    case 'efun_real'
        text = 'efun real';
    case 'deconv_abs'
        text = 'deconv|abs|';
    case 'deconv_real'
        text = 'deconv real';
end
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
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
    error('export_bold_top_xcorr_roi_bar_summaries:CopyPngFailed', ...
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


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
