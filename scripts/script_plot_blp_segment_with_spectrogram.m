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
time_range_sec = [300, 320];
spec_colormap = flipud(othercolor('Spectral10'));
freq_range_to_plot = [0, 250];
color_limits = [];
event_input = [];
event_colors = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880];

% Plot
plot_blp_segment_with_spectrogram( ...
    cfg, ...
    output_root, ...
    time_range_sec, ...
    spec_colormap, ...
    freq_range_to_plot, ...
    color_limits, ...
    event_input, ...
    event_colors);
