addpath(genpath('D:\Onedrive\util_functions\othercolor\'));

close all force;
clear functions;

cfg = cfg_F12m01();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

% Plot settings
output_root = 'D:\DataPons_processed\';
time_range_sec = [200, 210];
spec_colormap = flipud(othercolor('Spectral10'));
freq_range_to_plot = [0, 250];
color_limits = [];

% Plot
plot_blp_segment_with_spectrogram( ...
    cfg, ...
    output_root, ...
    time_range_sec, ...
    spec_colormap, ...
    freq_range_to_plot, ...
    color_limits);