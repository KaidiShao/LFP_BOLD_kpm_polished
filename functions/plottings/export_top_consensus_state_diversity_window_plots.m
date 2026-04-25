function P = export_top_consensus_state_diversity_window_plots(cfg, output_root, W, params)
%EXPORT_TOP_CONSENSUS_STATE_DIVERSITY_WINDOW_PLOTS Save top-window figures.

if nargin < 2 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 3 || isempty(W) || ~isstruct(W)
    error('W must be provided as a state-diversity result struct.');
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
prep_cfg.show_consensus = true;
prep_cfg.consensus_input = 'auto';

plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg);
t_raw = plot_data.t_raw;

if exist(params.save_dir, 'dir') ~= 7
    mkdir(params.save_dir);
end

plot_manifest = strings(height(top_windows), 1);

for i = 1:height(top_windows)
    idx1 = double(top_windows.global_start_idx(i));
    idx2 = double(top_windows.global_end_idx(i));

    if ismember('global_window_idx', top_windows.Properties.VariableNames)
        file_stub = sprintf('rank_%02d_globalwin_%03d', ...
            double(top_windows.state_diversity_rank(i)), ...
            double(top_windows.global_window_idx(i)));
    else
        file_stub = sprintf('rank_%02d_session_%02d_window_%03d', ...
            double(top_windows.state_diversity_rank(i)), ...
            double(top_windows.session_id(i)), ...
            double(top_windows.window_idx_in_session(i)));
    end

    png_file = fullfile(params.save_dir, [file_stub, '.png']);
    if params.skip_existing && exist(png_file, 'file') == 2
        plot_manifest(i) = string(png_file);
        continue;
    end

    if idx1 < 1 || idx2 > numel(t_raw) || idx2 < idx1
        error('Invalid global sample range [%d, %d] for row %d.', idx1, idx2, i);
    end

    time_range_sec = [t_raw(idx1), t_raw(idx2)];
    base_plot_cache = build_blp_plot_window_cache(plot_data, time_range_sec, params.freq_range_to_plot);
    hfig = plot_blp_segment_with_spectrogram(base_plot_cache, params.spec_colormap, params.color_limits);

    if ismember('global_window_idx', top_windows.Properties.VariableNames)
        title_line_1 = sprintf('Rank %d | globalwin %d | sessions %d->%d | samples [%d,%d]', ...
            double(top_windows.state_diversity_rank(i)), ...
            double(top_windows.global_window_idx(i)), ...
            double(top_windows.start_session_id(i)), ...
            double(top_windows.end_session_id(i)), ...
            idx1, idx2);
    else
        title_line_1 = sprintf('Rank %d | session %d | window %d | samples [%d,%d]', ...
            double(top_windows.state_diversity_rank(i)), ...
            double(top_windows.session_id(i)), ...
            double(top_windows.window_idx_in_session(i)), ...
            idx1, idx2);
    end

    title_line_2 = sprintf('richness=%d | total states=%d | dominant=%s | Hnorm=%.3f', ...
        double(top_windows.active_state_richness(i)), ...
        double(top_windows.total_state_window_count(i)), ...
        char(string(top_windows.dominant_state(i))), ...
        double(top_windows.normalized_state_entropy(i)));
    fig_title = sprintf('%s | %s', title_line_1, title_line_2);

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
    plot_manifest(i) = string(png_file);

    if params.close_after_save
        close(hfig);
    end
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
if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = fullfile(output_root, cfg.file_stem, 'consensus_state_diversity_windows', 'top_window_plots');
end

if ~isfield(params, 'save_png')
    params.save_png = true;
end

if ~isfield(params, 'close_after_save')
    params.close_after_save = true;
end

if ~isfield(params, 'skip_existing')
    params.skip_existing = false;
end

if ~isfield(params, 'fallback_size_px') || isempty(params.fallback_size_px)
    params.fallback_size_px = [4979, 2888];
end

if ~isfield(params, 'reference_dir')
    params.reference_dir = '';
end

if ~isfield(params, 'resolution') || isempty(params.resolution)
    params.resolution = 220;
end

if ~isfield(params, 'event_colors') || isempty(params.event_colors)
    params.event_colors = [ ...
        0.0000, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.4660, 0.6740, 0.1880];
end

if ~isfield(params, 'freq_range_to_plot') || isempty(params.freq_range_to_plot)
    params.freq_range_to_plot = [0, 250];
end

if ~isfield(params, 'color_limits')
    params.color_limits = [];
end

if ~isfield(params, 'spec_colormap') || isempty(params.spec_colormap)
    if exist('othercolor', 'file') == 2
        params.spec_colormap = flipud(othercolor('Spectral10'));
    else
        params.spec_colormap = flipud(turbo(256));
    end
end
end


function validate_window_result(W)
if ~isfield(W, 'window_mode') || ~strcmpi(char(string(W.window_mode)), 'global')
    error('State-diversity plotting currently expects a global-window result.');
end

required_vars = {'state_diversity_rank', 'global_start_idx', 'global_end_idx'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table is missing required column "%s".', required_vars{i});
    end
end
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
