function P = export_top_consensus_event_diversity_window_plots(cfg, output_root, W, params)
%EXPORT_TOP_CONSENSUS_EVENT_DIVERSITY_WINDOW_PLOTS Save top-window figures.

if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 3 || isempty(W) || ~isstruct(W)
    error('W must be provided as an event-diversity result struct.');
end

if nargin < 4
    params = struct();
end

params = apply_default_params(cfg, output_root, params);
validate_window_result(W);

top_windows = W.top_windows_table;
if isempty(top_windows)
    error('Top-window table is empty in %s.', W.save_file);
end

prep_cfg = struct();
prep_cfg.show_events = true;
prep_cfg.event_input = 'auto';
prep_cfg.band_colors = params.event_colors;
prep_cfg.include_spectrogram = true;
prep_cfg.show_consensus = params.show_consensus;
prep_cfg.consensus_input = 'auto';

plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg);
t_raw = plot_data.t_raw;

if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

plot_manifest = strings(height(top_windows), 1);
for i = 1:height(top_windows)
    plot_manifest(i) = string(local_export_one_top_window_plot( ...
        plot_data, top_windows, params, i, t_raw));
end

manifest_table = top_windows;
manifest_table.plot_file = plot_manifest;
manifest_file = fullfile(params.save_dir, 'plot_manifest.csv');
writetable(manifest_table, manifest_file);

P = struct();
P.save_dir = params.save_dir;
P.plot_manifest_file = manifest_file;
P.plot_manifest_table = manifest_table;
P.top_windows_table = top_windows;
end


function params = apply_default_params(cfg, output_root, params)
params_block = build_blp_consensus_state_pipeline_params(struct('plot_params', params));
params = params_block.plot_params;

if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', 'top_window_plots');
end

if ~isfield(params, 'show_consensus')
    params.show_consensus = true;
end
end


function validate_window_result(W)
required_vars = {'diversity_rank', 'global_start_idx', 'global_end_idx', 'total_event_count'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table is missing required column "%s".', required_vars{i});
    end
end
end


function png_file = local_export_one_top_window_plot(plot_data, top_windows, params, row_idx, t_raw)
idx1 = double(top_windows.global_start_idx(row_idx));
idx2 = double(top_windows.global_end_idx(row_idx));
png_file = fullfile(params.save_dir, [local_build_file_stub(top_windows, row_idx), '.png']);

if params.skip_existing && exist(png_file, 'file') == 2
    return;
end

if idx1 < 1 || idx2 > numel(t_raw) || idx2 < idx1
    error('Invalid global sample range [%d, %d] for row %d.', idx1, idx2, row_idx);
end

time_range_sec = [t_raw(idx1), t_raw(idx2)];
base_plot_cache = build_blp_plot_window_cache(plot_data, time_range_sec, params.freq_range_to_plot);
hfig = plot_blp_segment_with_spectrogram(base_plot_cache, params.spec_colormap, params.color_limits);

[title_line_1, title_line_2, fig_title] = local_build_plot_titles(top_windows, row_idx, idx1, idx2);
tl = findobj(hfig, 'Type', 'tiledlayout');
if ~isempty(tl)
    title(tl(1), {title_line_1, title_line_2}, 'Interpreter', 'none', ...
        'FontSize', 24, 'FontWeight', 'bold');
else
    set(hfig, 'Name', fig_title, 'NumberTitle', 'off');
end

set(hfig, 'Name', fig_title, 'NumberTitle', 'off');
target_size_px = local_get_export_size(params);
set(hfig, 'Units', 'pixels');
fig_pos = get(hfig, 'Position');
fig_pos(3:4) = target_size_px;
set(hfig, 'Position', fig_pos);
drawnow;

if params.save_png
    local_export_fixed_canvas_png(hfig, png_file, target_size_px, params.resolution);
end

if params.close_after_save
    close(hfig);
end
end


function file_stub = local_build_file_stub(top_windows, row_idx)
if ismember('global_window_idx', top_windows.Properties.VariableNames)
    file_stub = sprintf('rank_%02d_globalwin_%03d', ...
        double(top_windows.diversity_rank(row_idx)), ...
        double(top_windows.global_window_idx(row_idx)));
else
    file_stub = sprintf('rank_%02d_session_%02d_window_%03d', ...
        double(top_windows.diversity_rank(row_idx)), ...
        double(top_windows.session_id(row_idx)), ...
        double(top_windows.window_idx_in_session(row_idx)));
end
end


function [title_line_1, title_line_2, fig_title] = local_build_plot_titles(top_windows, row_idx, idx1, idx2)
if ismember('global_window_idx', top_windows.Properties.VariableNames)
    title_line_1 = sprintf('Rank %d | globalwin %d | sessions %d->%d | samples [%d,%d]', ...
        double(top_windows.diversity_rank(row_idx)), ...
        double(top_windows.global_window_idx(row_idx)), ...
        double(top_windows.start_session_id(row_idx)), ...
        double(top_windows.end_session_id(row_idx)), ...
        idx1, idx2);
else
    title_line_1 = sprintf('Rank %d | session %d | window %d | samples [%d,%d]', ...
        double(top_windows.diversity_rank(row_idx)), ...
        double(top_windows.session_id(row_idx)), ...
        double(top_windows.window_idx_in_session(row_idx)), ...
        idx1, idx2);
end

title_line_2 = sprintf('theta=%d | gamma=%d | ripple=%d | total=%d', ...
    double(top_windows.theta_count(row_idx)), ...
    double(top_windows.gamma_count(row_idx)), ...
    double(top_windows.ripple_count(row_idx)), ...
    double(top_windows.total_event_count(row_idx)));
fig_title = sprintf('%s | %s', title_line_1, title_line_2);
end


function target_size_px = local_get_export_size(params)
target_size_px = params.fallback_size_px;

reference_dir = char(string(params.reference_dir));
if isempty(reference_dir) || exist(reference_dir, 'dir') ~= 7
    return;
end

L = dir(fullfile(reference_dir, '*_ref.png'));
if isempty(L)
    return;
end

L = sortrows(struct2table(L), {'datenum', 'name'}, {'descend', 'ascend'});
reference_png = fullfile(L.folder{1}, L.name{1});
info = imfinfo(reference_png);
target_size_px = [double(info.Width), double(info.Height)];
end


function local_export_fixed_canvas_png(hfig, png_file, target_size_px, resolution)
paper_size_in = target_size_px ./ resolution;
set(hfig, 'PaperUnits', 'inches');
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition', [0, 0, paper_size_in(1), paper_size_in(2)]);
set(hfig, 'PaperSize', paper_size_in);
print(hfig, png_file, '-dpng', sprintf('-r%d', resolution));
end
