% Canonical e10gb1 diversity-plotting entry point.
% This script always uses the 6000-sample global-window diversity result:
%   e10gb1_event_diversity_windows_6000samp_globalwin.mat
% and writes plots into:
%   E:\DataPons_processed\e10gb1\event_diversity_windows\top_window_plots

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

cfg = cfg_E10gb1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = get_project_processed_root();
window_result_file = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', ...
    [cfg.file_stem, '_event_diversity_windows_6000samp_globalwin.mat']);

save_cfg = struct();
save_cfg.save_dir = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', 'top_window_plots');
save_cfg.save_png = true;
save_cfg.close_after_save = true;
save_cfg.skip_existing = false;
save_cfg.fallback_size_px = [4979, 2888];
save_cfg.reference_dir = fullfile(save_cfg.save_dir, 'ref_figs');

event_colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];

freq_range_to_plot = [0, 250];
color_limits = [];

if exist('othercolor', 'file') == 2
    spec_colormap = flipud(othercolor('Spectral10'));
else
    spec_colormap = flipud(turbo(256));
end

S = load(window_result_file, 'W');
W = S.W;
local_validate_e10gb1_diversity_result(W, window_result_file);
top_windows = W.top_windows_table;

if isempty(top_windows)
    error('Top-window table is empty in %s.', window_result_file);
end

prep_cfg = struct();
prep_cfg.show_events = true;
prep_cfg.event_input = 'auto';
prep_cfg.band_colors = event_colors;
prep_cfg.include_spectrogram = true;
prep_cfg.show_consensus = true;
prep_cfg.consensus_input = 'auto';

plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg);
t_raw = plot_data.t_raw;

if exist(save_cfg.save_dir, 'dir') ~= 7
    mkdir(save_cfg.save_dir);
end

plot_manifest = strings(height(top_windows), 1);

for i = 1:height(top_windows)
    idx1 = double(top_windows.global_start_idx(i));
    idx2 = double(top_windows.global_end_idx(i));

    file_stub = sprintf('rank_%02d_globalwin_%03d', ...
        double(top_windows.diversity_rank(i)), ...
        double(top_windows.global_window_idx(i)));
    png_file = fullfile(save_cfg.save_dir, [file_stub, '.png']);

    if save_cfg.skip_existing && exist(png_file, 'file') == 2
        plot_manifest(i) = string(png_file);
        continue;
    end

    if idx1 < 1 || idx2 > numel(t_raw) || idx2 < idx1
        error('Invalid global sample range [%d, %d] for row %d.', idx1, idx2, i);
    end

    time_range_sec = [t_raw(idx1), t_raw(idx2)];
    base_plot_cache = build_blp_plot_window_cache(plot_data, time_range_sec, freq_range_to_plot);

    hfig = plot_blp_segment_with_spectrogram(base_plot_cache, spec_colormap, color_limits);

    title_line_1 = sprintf('Rank %d | globalwin %d | sessions %d->%d | samples [%d,%d]', ...
        double(top_windows.diversity_rank(i)), ...
        double(top_windows.global_window_idx(i)), ...
        double(top_windows.start_session_id(i)), ...
        double(top_windows.end_session_id(i)), ...
        idx1, ...
        idx2);
    title_line_2 = sprintf('theta=%d | gamma=%d | ripple=%d | total=%d', ...
        double(top_windows.theta_count(i)), ...
        double(top_windows.gamma_count(i)), ...
        double(top_windows.ripple_count(i)), ...
        double(top_windows.total_event_count(i)));
    fig_title = sprintf('%s | %s', title_line_1, title_line_2);

    tl = findobj(hfig, 'Type', 'tiledlayout');
    if ~isempty(tl)
        title(tl(1), {title_line_1, title_line_2}, 'Interpreter', 'none', ...
            'FontSize', 24, 'FontWeight', 'bold');
    else
        set(hfig, 'Name', fig_title, 'NumberTitle', 'off');
    end

    set(hfig, 'Name', fig_title, 'NumberTitle', 'off');
    target_size_px = local_get_reference_export_size(save_cfg);
    set(hfig, 'Units', 'pixels');
    fig_pos = get(hfig, 'Position');
    fig_pos(3:4) = target_size_px;
    set(hfig, 'Position', fig_pos);
    drawnow;

    local_export_fixed_canvas_png(hfig, png_file, target_size_px, 220);
    plot_manifest(i) = string(png_file);

    if save_cfg.close_after_save
        close(hfig);
    end
end

manifest_table = top_windows;
manifest_table.plot_file = plot_manifest;
writetable(manifest_table, fullfile(save_cfg.save_dir, 'plot_manifest.csv'));

fprintf('Saved %d top-window plots to:\n  %s\n', height(top_windows), save_cfg.save_dir);


function target_size_px = local_get_reference_export_size(save_cfg)
target_size_px = save_cfg.fallback_size_px;

reference_png = local_find_reference_png(save_cfg.reference_dir);
if ~isempty(reference_png)
    info = imfinfo(reference_png);
    target_size_px = [double(info.Width), double(info.Height)];
end
end


function reference_png = local_find_reference_png(reference_dir)
reference_png = '';
if exist(reference_dir, 'dir') ~= 7
    return;
end

L = dir(fullfile(reference_dir, '*_ref.png'));
if isempty(L)
    return;
end

L = sortrows(struct2table(L), {'datenum', 'name'}, {'descend', 'ascend'});
reference_png = fullfile(L.folder{1}, L.name{1});
end


function local_validate_e10gb1_diversity_result(W, window_result_file)
if ~isfield(W, 'window_mode') || ~strcmpi(char(string(W.window_mode)), 'global')
    error(['%s is not a global-window diversity result. ', ...
        'Use the canonical 6000-sample global-window file for e10gb1.'], window_result_file);
end

if ~isfield(W, 'window_length_samples') || double(W.window_length_samples) ~= 6000
    error(['%s does not use 6000-sample windows. ', ...
        'Use the canonical 6000-sample global-window file for e10gb1.'], window_result_file);
end

required_vars = {'global_window_idx', 'start_session_id', 'end_session_id'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table in %s is missing required column "%s".', ...
            window_result_file, required_vars{i});
    end
end
end


function local_export_fixed_canvas_png(hfig, png_file, target_size_px, resolution)
paper_size_in = target_size_px ./ resolution;
set(hfig, 'PaperUnits', 'inches');
set(hfig, 'PaperPositionMode', 'manual');
set(hfig, 'PaperPosition', [0, 0, paper_size_in(1), paper_size_in(2)]);
set(hfig, 'PaperSize', paper_size_in);
print(hfig, png_file, '-dpng', sprintf('-r%d', resolution));
end
