function result = export_bold_dimred_top_xcorr_roi_bar_summaries(dimred_result_input, bold_post_input, xcorr_input, params)
%EXPORT_BOLD_DIMRED_TOP_XCORR_ROI_BAR_SUMMARIES Export P10 component ROI summaries.

if nargin < 1 || isempty(dimred_result_input)
    error('dimred_result_input is required.');
end
if nargin < 2
    bold_post_input = [];
end
if nargin < 3 || isempty(xcorr_input)
    error('xcorr_input is required.');
end
if nargin < 4 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);
t_run = tic;

[component_ctx, component_info, ~, dimred_result] = ...
    build_bold_dimred_component_plot_context(dimred_result_input, bold_post_input, params);
xcorr_out = load_bold_xcorr_output(xcorr_input);
group_specs = build_bold_top_xcorr_group_specs(xcorr_out, params);
xcorr_source_file = local_get_field(xcorr_out, 'source_file', '');
clear xcorr_out
run_info = component_ctx.run_info;

if isempty(group_specs)
    result = local_make_result(run_info, 'no_top_table', ...
        'XCORR top tables are empty.', '', '', 0, 0, toc(t_run));
    result.group_results = repmat(local_empty_group_result(), 0, 1);
    return;
end

out_root = local_output_root(params, run_info);
mat_dir = fullfile(out_root, 'mat');
fig_dir = fullfile(out_root, 'fig');
if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end

group_results = repmat(local_empty_group_result(), numel(group_specs), 1);
primary_info_file = '';
primary_summary_dir = '';
total_selected = 0;
total_saved = 0;

for i_group = 1:numel(group_specs)
    group_spec = group_specs(i_group);
    summary_dir = local_group_output_dir(fig_dir, group_spec);
    if exist(summary_dir, 'dir') ~= 7, mkdir(summary_dir); end

    [component_indices, top_row_indices] = resolve_bold_dimred_top_xcorr_component_indices( ...
        group_spec.top_table, component_info.n_components, params.top_n);
    if isempty(component_indices)
        group_results(i_group) = local_make_group_result(group_spec, ...
            'no_component_indices', ...
            sprintf('%s top table did not map to any component index.', group_spec.display_name), ...
            summary_dir, '', '', 0, 0, [], table());
        continue;
    end

    selected_top_table = group_spec.top_table(top_row_indices, :);
    stem = sprintf('dimred_roi__%s_top%d', group_spec.slug, params.top_n);
    png_file = fullfile(summary_dir, [stem, '.png']);
    fig_file = fullfile(summary_dir, [stem, '.fig']);
    info_dir = local_group_info_dir(mat_dir, group_spec);
    if exist(info_dir, 'dir') ~= 7, mkdir(info_dir); end
    info_file = fullfile(info_dir, [stem, '_info.mat']);

    plot_out = struct();
    if params.skip_existing && exist(png_file, 'file') == 2 && ...
            (~params.save_fig || exist(fig_file, 'file') == 2) && ...
            exist(info_file, 'file') == 2
        n_saved = 0;
        group_status = 'skipped_existing';
        group_message = sprintf('%s ROI summary already existed.', group_spec.display_name);
    else
        plot_params = params;
        plot_params.selection_mode = 'raw';
        plot_params.basis_indices = component_indices;
        plot_params.mode_title_texts = local_build_mode_title_texts( ...
            selected_top_table, component_indices, group_spec, dimred_result);
        plot_params.figure_title_text = local_build_figure_title_text( ...
            run_info, group_spec, dimred_result);
        [fig, plot_out] = plot_bold_mode_roi_bar_summary(component_ctx, plot_params);
        if params.save_png
            exportgraphics(fig, png_file, 'Resolution', params.resolution);
        end
        if params.save_fig
            savefig(fig, fig_file);
        end
        close(fig);
        drawnow limitrate;
        n_saved = 1;
        group_status = 'ok';
        group_message = sprintf('%s ROI summary saved with %d component(s).', ...
            group_spec.display_name, numel(component_indices));
    end

    summary_info = struct();
    summary_info.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
    summary_info.run_info = run_info;
    summary_info.dimred_result_file = component_info.dimred_result_file;
    summary_info.xcorr_file = xcorr_source_file;
    summary_info.selection_scope = group_spec.scope;
    summary_info.selection_name = group_spec.display_name;
    summary_info.selection_slug = group_spec.slug;
    summary_info.component_indices = component_indices;
    summary_info.selected_top_table = selected_top_table;
    summary_info.png_file = png_file;
    summary_info.fig_file = fig_file;
    summary_info.summary_dir = summary_dir;
    summary_info.component_info = component_info;
    summary_info.plot_out = plot_out;
    save_mat_variable_atomic(info_file, 'summary_info', summary_info);

    if isempty(primary_info_file) || strcmp(group_spec.scope, 'combined')
        primary_info_file = info_file;
        primary_summary_dir = summary_dir;
    end
    total_selected = total_selected + numel(component_indices);
    total_saved = total_saved + n_saved;
    group_results(i_group) = local_make_group_result(group_spec, group_status, ...
        group_message, summary_dir, png_file, info_file, n_saved, ...
        numel(component_indices), component_indices, selected_top_table);
end

status = local_combine_group_status(group_results);
message = strjoin(cellstr(string({group_results.message})), ' | ');
result = local_make_result(run_info, status, message, primary_summary_dir, ...
    primary_info_file, total_selected, total_saved, toc(t_run));
result.output_root = out_root;
result.group_results = group_results;
result.summary_group_dirs = {group_results.summary_dir};
result.info_files = {group_results.info_file};
result.dimred_result_file = component_info.dimred_result_file;
result.xcorr_file = xcorr_source_file;
end


function params = local_apply_defaults(params)
base = build_bold_dimred_cross_modal_coupling_params();
params = local_merge_defaults(base.roi_summary, params);
params = local_set_default(params, 'processed_root', io_project.get_project_processed_root());
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'output_root', '');
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
    end
end
end


function out_root = local_output_root(params, run_info)
if isfield(params, 'output_root') && ~isempty(params.output_root)
    out_root = char(string(params.output_root));
else
    out_root = fullfile(io_project.get_pipeline_stage_dir(params.processed_root, ...
        run_info.dataset_stem, 10, 'figures_bold_dimred_top_xcorr_activation_maps'), ...
        local_get_field(run_info, 'run_tag', run_info.run_name));
end
end


function texts = local_build_mode_title_texts(selected_top_table, component_indices, group_spec, dimred_result)
texts = cell(1, numel(component_indices));
feature_name = local_get_field(local_get_field(dimred_result, 'feature', struct()), 'name', '');
method_name = local_get_field(local_get_field(dimred_result, 'meta', struct()), 'method', '');
for i = 1:numel(component_indices)
    row = selected_top_table(i, :);
    peak_corr = local_table_scalar(row, 'peak_corr', NaN);
    peak_lag = local_table_scalar(row, 'peak_lag_sec', NaN);
    texts{i} = sprintf('top%02d | comp %d\\newline%s/%s\\newlinecorr %.3f | lag %.1fs', ...
        i, component_indices(i), feature_name, method_name, peak_corr, peak_lag);
end
if nargin >= 3 && ~isempty(group_spec.display_name)
    texts = cellfun(@(x) sprintf('%s\\newline%s', group_spec.display_name, x), ...
        texts, 'UniformOutput', false);
end
end


function title_text = local_build_figure_title_text(run_info, group_spec, dimred_result)
feature_name = local_get_field(local_get_field(dimred_result, 'feature', struct()), 'name', '');
method_name = local_get_field(local_get_field(dimred_result, 'meta', struct()), 'method', '');
title_text = sprintf('%s | %s/%s | P10 ROI component summary | %s', ...
    local_get_field(run_info, 'run_name', 'BOLD'), feature_name, method_name, ...
    group_spec.display_name);
end


function value = local_table_scalar(T, name, default_value)
if istable(T) && ismember(name, T.Properties.VariableNames)
    value = double(T.(name)(1));
else
    value = default_value;
end
end


function out_dir = local_group_output_dir(fig_dir, group_spec)
family = local_group_feature_family(group_spec);
if isempty(family)
    out_dir = fullfile(fig_dir, 'roi');
else
    out_dir = fullfile(fig_dir, family, 'roi');
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
if contains(lower(char(string(feature_name))), 'deconv')
    family = 'deconv_dimred';
elseif ~isempty(feature_name)
    family = 'efun_dimred';
end
end


function result = local_make_result(run_info, status, message, summary_dir, info_file, n_selected, n_saved, runtime_sec)
result = struct();
result.run_info = run_info;
result.dataset_stem = local_get_field(run_info, 'dataset_stem', '');
result.run_name = local_get_field(run_info, 'run_name', '');
result.status = status;
result.message = message;
result.summary_dir = summary_dir;
result.info_file = info_file;
result.n_selected_modes = n_selected;
result.n_saved_figures = n_saved;
result.runtime_sec = runtime_sec;
end


function group_result = local_empty_group_result()
group_result = struct('scope', '', 'display_name', '', 'slug', '', ...
    'density_name', '', 'feature_name', '', 'status', '', 'message', '', ...
    'summary_dir', '', 'png_file', '', 'info_file', '', ...
    'n_saved_figures', 0, 'n_selected_modes', 0, ...
    'component_indices', [], 'selected_top_table', table());
end


function group_result = local_make_group_result(group_spec, status, message, ...
        summary_dir, png_file, info_file, n_saved, n_selected, component_indices, selected_top_table)
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
group_result.n_saved_figures = n_saved;
group_result.n_selected_modes = n_selected;
group_result.component_indices = component_indices;
group_result.selected_top_table = selected_top_table;
end


function status = local_combine_group_status(group_results)
statuses = string({group_results.status});
if any(statuses == "ok")
    status = 'ok';
elseif any(statuses == "skipped_existing")
    status = 'skipped_existing';
elseif any(statuses == "no_component_indices")
    status = 'no_component_indices';
else
    status = 'empty';
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
