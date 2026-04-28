function base_plot_cache = build_blp_plot_window_cache(plot_data, time_range_sec, freq_lim)
%BUILD_BLP_PLOT_WINDOW_CACHE Slice preloaded plot data for one time window.

if nargin < 2 || isempty(time_range_sec)
    error('time_range_sec must be provided as [t_start, t_end].');
end

if nargin < 3 || isempty(freq_lim)
    freq_lim = [0 250];
end

if numel(time_range_sec) ~= 2 || time_range_sec(2) <= time_range_sec(1)
    error('time_range_sec must be [t_start, t_end] with t_end > t_start.');
end

if ~isstruct(plot_data) || ~isfield(plot_data, 'D') || ~isfield(plot_data, 't_raw')
    error('plot_data must come from prepare_blp_plot_data.');
end

t1 = double(time_range_sec(1));
t2 = double(time_range_sec(2));

[D, t_raw, raw_idx, x_seg, t_seg, idx_plot_1, idx_plot_2, border_t_raw] = ...
    local_slice_raw_window(plot_data, t1, t2);
[spec_file, regions, spec_seg, t_spec, f_spec, spec_global_clim] = ...
    local_slice_spectrogram_window(plot_data, raw_idx, t_seg, freq_lim);
[n_channels, x_disp, offsets, ytick_pos, ytick_labels] = ...
    local_build_trace_layout(plot_data, D, x_seg);
[show_consensus, consensus_windows, consensus_state_catalog, consensus_colors, consensus_face_alpha] = ...
    local_collect_consensus_windows(plot_data, idx_plot_1, idx_plot_2);
plot_settings = plot_data.plot_settings;

base_plot_cache = struct();
base_plot_cache.cfg = plot_data.cfg;
base_plot_cache.output_root = plot_data.output_root;
base_plot_cache.plot_settings = plot_settings;
base_plot_cache.time_range_sec = [t1, t2];

base_plot_cache.D = D;
base_plot_cache.t_raw = t_raw;
base_plot_cache.raw_idx = raw_idx;
base_plot_cache.x_seg = x_seg;
base_plot_cache.t_seg = t_seg;
base_plot_cache.idx_plot_1 = idx_plot_1;
base_plot_cache.idx_plot_2 = idx_plot_2;

base_plot_cache.R_events = plot_data.R_events;
base_plot_cache.n_event_bands = plot_data.n_event_bands;
base_plot_cache.band_colors = plot_data.band_colors;
base_plot_cache.show_events = plot_data.show_events;

base_plot_cache.show_consensus = show_consensus;
base_plot_cache.consensus_windows = consensus_windows;
base_plot_cache.consensus_state_catalog = consensus_state_catalog;
base_plot_cache.consensus_colors = consensus_colors;
base_plot_cache.consensus_face_alpha = consensus_face_alpha;

base_plot_cache.spec_file = spec_file;
base_plot_cache.regions = regions;
base_plot_cache.spec_seg = spec_seg;
base_plot_cache.t_spec = t_spec;
base_plot_cache.f_spec = f_spec;
base_plot_cache.freq_lim = freq_lim;
base_plot_cache.spec_global_clim = spec_global_clim;

base_plot_cache.n_channels = n_channels;
base_plot_cache.x_disp = x_disp;
base_plot_cache.offsets = offsets;
base_plot_cache.ytick_pos = ytick_pos;
base_plot_cache.ytick_labels = ytick_labels;
base_plot_cache.border_t_raw = border_t_raw;
end


function [D, t_raw, raw_idx, x_seg, t_seg, idx_plot_1, idx_plot_2, border_t_raw] = local_slice_raw_window(plot_data, t1, t2)
D = plot_data.D;
t_raw = double(plot_data.t_raw(:));
raw_idx = find(t_raw >= t1 & t_raw <= t2);
if isempty(raw_idx)
    error('The requested time range does not overlap the raw data.');
end

x_seg = double(D.data(raw_idx, :));
t_seg = double(t_raw(raw_idx));
idx_plot_1 = raw_idx(1);
idx_plot_2 = raw_idx(end);

border_t_raw = double(t_raw(D.border_idx));
border_t_raw = border_t_raw(border_t_raw >= t1 & border_t_raw <= t2);
end


function [spec_file, regions, spec_seg, t_spec, f_spec, spec_global_clim] = local_slice_spectrogram_window(plot_data, raw_idx, t_seg, freq_lim)
spec_file = plot_data.spec_file;
regions = {};
spec_seg = [];
t_spec = [];
f_spec = [];
spec_global_clim = [];

if ~plot_data.include_spectrogram
    return;
end

freqs = double(plot_data.spec_freqs(:));
regions = plot_data.spec_regions;
freq_idx = find(freqs >= freq_lim(1) & freqs <= freq_lim(2));
if isempty(freq_idx)
    error('No frequencies fall inside the requested freq_lim.');
end

M = matfile(spec_file);
spec_seg = double(M.tmpall_mean_abs(freq_idx, raw_idx, :));
t_spec = double(t_seg(:).');
f_spec = double(freqs(freq_idx));

if isfield(plot_data, 'spec_global_clim')
    spec_global_clim = plot_data.spec_global_clim;
end
end


function [n_channels, x_disp, offsets, ytick_pos, ytick_labels] = local_build_trace_layout(plot_data, D, x_seg)
plot_settings = plot_data.plot_settings;
channel_sites = plot_data.cfg.channels.sites(D.selected_channels);
n_channels = numel(D.selected_channels);

x_disp = zeros(size(x_seg));
for c = 1:n_channels
    xc = x_seg(:, c);
    xc = xc - median(xc, 'omitnan');

    s = mad(xc, 1);
    if ~isfinite(s) || s == 0
        s = std(xc, 0, 'omitnan');
    end
    if ~isfinite(s) || s == 0
        s = 1;
    end

    xc = xc / s;
    xc = max(min(xc, plot_settings.trace_clip), -plot_settings.trace_clip);
    x_disp(:, c) = plot_settings.trace_scale * xc;
end

offsets = zeros(1, n_channels);
offset_val = 0;
offsets(n_channels) = offset_val;

for c = n_channels-1:-1:1
    if strcmp(channel_sites{c}, channel_sites{c+1})
        offset_val = offset_val + plot_settings.within_gap;
    else
        offset_val = offset_val + plot_settings.between_gap;
    end
    offsets(c) = offset_val;
end

offsets = double(offsets);
channel_labels = cell(1, n_channels);
for c = 1:n_channels
    channel_labels{c} = sprintf('Ch. %d', D.selected_channels(c));
end

[ytick_pos, tick_order] = sort(offsets, 'ascend');
ytick_labels = channel_labels(tick_order);
end


function [show_consensus, consensus_windows, consensus_state_catalog, consensus_colors, consensus_face_alpha] = ...
    local_collect_consensus_windows(plot_data, idx_plot_1, idx_plot_2)
show_consensus = isfield(plot_data, 'show_consensus') && plot_data.show_consensus;
consensus_windows = struct([]);
consensus_state_catalog = struct([]);
consensus_colors = [];
consensus_face_alpha = [];

if ~show_consensus
    return;
end

C_consensus = plot_data.C_consensus;
consensus_state_catalog = C_consensus.state_catalog;
consensus_colors = plot_data.consensus_colors;
consensus_face_alpha = plot_data.consensus_face_alpha;

if isempty(C_consensus.state_windows)
    return;
end

all_windows = C_consensus.state_windows(:);
keep_mask = false(numel(all_windows), 1);
for i = 1:numel(all_windows)
    win = double(all_windows(i).win_global);
    keep_mask(i) = win(2) >= idx_plot_1 && win(1) <= idx_plot_2;
end

consensus_windows = all_windows(keep_mask);
end
