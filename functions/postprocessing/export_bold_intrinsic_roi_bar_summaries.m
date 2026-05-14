function [result, B] = export_bold_intrinsic_roi_bar_summaries(bold_post_input, params)
%EXPORT_BOLD_INTRINSIC_ROI_BAR_SUMMARIES Export pipeline 7 ROI-average Koopman summaries.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

t_run = tic;
[plot_ctx, B] = build_bold_activation_plot_context(bold_post_input, struct( ...
    'datapons_root', params.datapons_root, ...
    'roi_ts_file', params.roi_ts_file, ...
    'feature_reduce', params.feature_reduce));
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
summary_dir = fullfile(fig_dir, 'intrinsic_roi_bar_summaries');
if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end
if exist(summary_dir, 'dir') ~= 7, mkdir(summary_dir); end

stem = local_build_summary_stem(params);
png_file = fullfile(summary_dir, [stem, '.png']);
fig_file = fullfile(summary_dir, [stem, '.fig']);
info_file = fullfile(mat_dir, [stem, '_info.mat']);

if params.skip_existing && exist(png_file, 'file') == 2 && ...
        (~params.save_fig || exist(fig_file, 'file') == 2) && ...
        exist(info_file, 'file') == 2
    plot_out = local_load_existing_plot_out(info_file);
    n_saved = 0;
    status = 'skipped_existing';
    message = 'Intrinsic ROI bar summary already existed.';
else
    plot_params = params;
    [fig, plot_out] = plot_bold_mode_roi_bar_summary(plot_ctx, plot_params);
    if params.save_png
        exportgraphics(fig, png_file, 'Resolution', params.resolution);
    end
    if params.save_fig
        savefig(fig, fig_file);
    end
    close(fig);
    drawnow limitrate;
    n_saved = 1;
    status = 'ok';
    message = sprintf('Saved intrinsic ROI bar summary with %d mode(s).', ...
        numel(local_get_field(plot_out, 'raw_indices', [])));
end

summary_info = struct();
summary_info.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
summary_info.run_info = run_info;
summary_info.bold_post_file = local_get_field(B, 'source_file', '');
summary_info.observable_file = local_get_field(plot_ctx, 'observable_file', '');
summary_info.roi_ts_file = plot_ctx.roi_ts_file;
summary_info.png_file = png_file;
summary_info.fig_file = fig_file;
summary_info.summary_dir = summary_dir;
summary_info.params = params;
summary_info.plot_out = plot_out;
save(info_file, 'summary_info', '-v7.3');

B.artifacts.intrinsic_roi_summary_dir = summary_dir;
B.artifacts.intrinsic_roi_summary_png = png_file;
B.artifacts.intrinsic_roi_summary_info_file = info_file;
if params.update_bold_post && isfield(B, 'source_file') && ~isempty(B.source_file)
    BOLD_POST = B;
    save_bold_post_mat(B.source_file, BOLD_POST);
end

result = struct();
result.dataset_stem = run_info.dataset_stem;
result.run_name = run_info.run_name;
result.observable_mode = local_get_field(run_info, 'observable_mode', '');
result.residual_form = local_get_field(run_info, 'residual_form', '');
result.status = status;
result.message = message;
result.summary_dir = summary_dir;
result.png_file = png_file;
result.fig_file = fig_file;
result.info_file = info_file;
result.selection_mode = params.selection_mode;
result.basis_indices = params.basis_indices(:).';
result.raw_indices = local_get_field(plot_out, 'raw_indices', []);
result.n_selected_modes = numel(result.raw_indices);
result.n_saved_figures = n_saved;
result.runtime_sec = toc(t_run);
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'processed_root', io_project.get_project_processed_root());
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'output_root', '');
params = local_set_default(params, 'roi_ts_file', '');
params = local_set_default(params, 'selection_mode', 'sorted');
params = local_set_default(params, 'basis_indices', 1:5);
params = local_set_default(params, 'feature_reduce', 'mean');
params = local_set_default(params, 'roi_reduce', 'mean');
params = local_set_default(params, 'roi_value_mode', 'mean_abs');
params = local_set_default(params, 'mode_normalization', 'range');
params = local_set_default(params, 'flip_roi_order', true);
params = local_set_default(params, 'layout_mode', 'row');
params = local_set_default(params, 'skip_existing', true);
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'resolution', 220);
params = local_set_default(params, 'update_bold_post', true);
end


function stem = local_build_summary_stem(params)
stem = sprintf('intrinsic_roi_bar_summary__%s__%s__%s__%s', ...
    local_build_basis_tag(params.selection_mode, params.basis_indices), ...
    char(string(params.roi_value_mode)), ...
    char(string(params.mode_normalization)), ...
    char(string(params.layout_mode)));
end


function tag = local_build_basis_tag(selection_mode, basis_indices)
basis_indices = unique(round(double(basis_indices(:).')), 'stable');
basis_indices = basis_indices(isfinite(basis_indices) & basis_indices >= 1);
if isempty(basis_indices)
    tag = [char(string(selection_mode)) '_idxnone'];
    return;
end
body = strjoin(arrayfun(@(x) sprintf('%03d', x), basis_indices, 'UniformOutput', false), '-');
tag = sprintf('%s_idx%s', char(string(selection_mode)), body);
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
