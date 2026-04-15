this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
results_root = get_project_results_root(repo_root);
addpath(genpath(repo_root));

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

cfg = cfg_F12m01();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = 'D:\DataPons_processed\';
window_result_file = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', ...
    [cfg.file_stem, '_event_diversity_windows_5000samp.mat']);

save_cfg = struct();
save_cfg.save_dir = fullfile(results_root, 'top_consensus_event_diversity_windows');
save_cfg.save_png = true;
save_cfg.close_after_save = true;
save_cfg.skip_existing = true;

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
top_windows = W.top_windows_table;

if isempty(top_windows)
    error('Top-window table is empty in %s.', window_result_file);
end

prep_cfg = struct();
prep_cfg.show_events = true;
prep_cfg.event_input = 'auto';
prep_cfg.band_colors = event_colors;
prep_cfg.include_spectrogram = true;

plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg);
t_raw = plot_data.t_raw;

if exist(save_cfg.save_dir, 'dir') ~= 7
    mkdir(save_cfg.save_dir);
end

plot_manifest = strings(height(top_windows), 1);

for i = 1:height(top_windows)
    idx1 = double(top_windows.global_start_idx(i));
    idx2 = double(top_windows.global_end_idx(i));

    file_stub = sprintf('rank_%02d_session_%02d_window_%03d', ...
        double(top_windows.diversity_rank(i)), ...
        double(top_windows.session_id(i)), ...
        double(top_windows.window_idx_in_session(i)));
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

    fig_title = sprintf(['Rank %d | session %d window %d | samples [%d,%d] | ' ...
        'theta=%d gamma=%d ripple=%d | total=%d'], ...
        double(top_windows.diversity_rank(i)), ...
        double(top_windows.session_id(i)), ...
        double(top_windows.window_idx_in_session(i)), ...
        idx1, ...
        idx2, ...
        double(top_windows.theta_count(i)), ...
        double(top_windows.gamma_count(i)), ...
        double(top_windows.ripple_count(i)), ...
        double(top_windows.total_event_count(i)));

    tl = findobj(hfig, 'Type', 'tiledlayout');
    if ~isempty(tl)
        title(tl(1), fig_title, 'Interpreter', 'none');
    else
        set(hfig, 'Name', fig_title, 'NumberTitle', 'off');
    end

    set(hfig, 'Name', fig_title, 'NumberTitle', 'off');
    drawnow;

    exportgraphics(hfig, png_file, 'Resolution', 220);
    plot_manifest(i) = string(png_file);

    if save_cfg.close_after_save
        close(hfig);
    end
end

manifest_table = top_windows;
manifest_table.plot_file = plot_manifest;
writetable(manifest_table, fullfile(save_cfg.save_dir, 'plot_manifest.csv'));

fprintf('Saved %d top-window plots to:\n  %s\n', height(top_windows), save_cfg.save_dir);
