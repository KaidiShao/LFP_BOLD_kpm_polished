close all force;

cfg = cfg_F12m01();

cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;
cfg.plot.event_linewidth = 1.1;

output_root = get_project_processed_root();
time_range_sec = [200, 210];
band_colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];

prep_cfg = struct();
prep_cfg.show_events = true;
prep_cfg.event_input = 'auto';
prep_cfg.band_colors = band_colors;
prep_cfg.include_spectrogram = false;

plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg);
base_plot_cache = build_blp_plot_window_cache(plot_data, time_range_sec);

plot_blp_segment_with_events(base_plot_cache);
