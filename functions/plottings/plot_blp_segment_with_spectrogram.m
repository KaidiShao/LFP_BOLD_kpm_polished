function hfig = plot_blp_segment_with_spectrogram(cfg, output_root, time_range_sec, cmap_in, freq_lim, clim)
% Plot raw BLP traces and saved region-mean spectrograms for a selected time range.

%% =========================
%  Input defaults
%  =========================
if nargin < 2 || isempty(output_root)
    output_root = 'D:\DataPons_processed\';
end

if nargin < 3 || isempty(time_range_sec)
    error('time_range_sec must be provided as [t_start, t_end].');
end

if nargin < 4 || isempty(cmap_in)
    cmap_in = parula(256);
end

if nargin < 5 || isempty(freq_lim)
    freq_lim = [0 250];
end

if nargin < 6
    clim = [];
end

if numel(time_range_sec) ~= 2 || time_range_sec(2) <= time_range_sec(1)
    error('time_range_sec must be [t_start, t_end] with t_end > t_start.');
end

t1 = double(time_range_sec(1));
t2 = double(time_range_sec(2));

%% =========================
%  Plot settings
%  =========================
trace_scale = 0.18;
trace_clip = 4;
within_gap = 1.4;
between_gap = 2.2;
trace_linewidth = 0.4;

if isfield(cfg, 'plot')
    if isfield(cfg.plot, 'trace_scale')
        trace_scale = cfg.plot.trace_scale;
    end
    if isfield(cfg.plot, 'trace_clip')
        trace_clip = cfg.plot.trace_clip;
    end
    if isfield(cfg.plot, 'within_gap')
        within_gap = cfg.plot.within_gap;
    end
    if isfield(cfg.plot, 'between_gap')
        between_gap = cfg.plot.between_gap;
    end
    if isfield(cfg.plot, 'trace_linewidth')
        trace_linewidth = cfg.plot.trace_linewidth;
    end
end

%% =========================
%  Load raw data
%  =========================
D = load_blp_dataset(cfg);

t_raw = build_global_time_axis(D.session_lengths, D.session_dx);
t_raw = double(t_raw(:));

raw_idx = find(t_raw >= t1 & t_raw <= t2);
if isempty(raw_idx)
    error('The requested time range does not overlap the raw data.');
end

x_seg = double(D.data(raw_idx, :));
t_seg = double(t_raw(raw_idx));

%% =========================
%  Locate saved spectrogram file
%  =========================
pad_sec = 20;
pad_mode = 'mirror';

if isfield(cfg, 'spectrogram')
    if isfield(cfg.spectrogram, 'pad_sec')
        pad_sec = cfg.spectrogram.pad_sec;
    end
    if isfield(cfg.spectrogram, 'pad_mode')
        pad_mode = cfg.spectrogram.pad_mode;
    end
end

spec_file = find_regionmean_spectrogram_file(cfg, output_root, pad_mode, pad_sec);

%% =========================
%  Load only the needed spectrogram segment
%  =========================
meta = load(spec_file, 'freqs', 'regions');
freqs = double(meta.freqs(:));
regions = meta.regions;

freq_idx = find(freqs >= freq_lim(1) & freqs <= freq_lim(2));
if isempty(freq_idx)
    error('No frequencies fall inside the requested freq_lim.');
end

spec_idx = raw_idx;

M = matfile(spec_file);
spec_seg = double(M.tmpall_mean_abs(freq_idx, spec_idx, :));

t_spec = double(t_seg(:).');
f_spec = double(freqs(freq_idx));

%% =========================
%  Prepare raw-trace display
%  =========================
channel_sites = cfg.channels.sites(D.selected_channels);
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
    xc = max(min(xc, trace_clip), -trace_clip);

    x_disp(:, c) = trace_scale * xc;
end

offsets = zeros(1, n_channels);
offset_val = 0;
offsets(n_channels) = offset_val;

for c = n_channels-1:-1:1
    if strcmp(channel_sites{c}, channel_sites{c+1})
        offset_val = offset_val + within_gap;
    else
        offset_val = offset_val + between_gap;
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

%% =========================
%  Build figure
%  =========================
n_regions = size(spec_seg, 3);

hfig = figure('Color', 'w', 'Position', [100, 100, 1000, 260 + 220 * n_regions]);
tl = tiledlayout(n_regions + 1, 1, 'TileSpacing', 'compact', 'Padding', 'compact');

ax_all = gobjects(n_regions + 1, 1);

%% =========================
%  Plot raw traces
%  =========================
ax1 = nexttile(tl, 1);
ax_all(1) = ax1;
disableDefaultInteractivity(ax1);
hold(ax1, 'on');

for c = 1:n_channels
    plot(ax1, t_seg, x_disp(:, c) + offsets(c), 'k', 'LineWidth', trace_linewidth);
end

border_t_raw = double(t_raw(D.border_idx));
border_t_raw = border_t_raw(border_t_raw >= t1 & border_t_raw <= t2);

for i = 1:numel(border_t_raw)
    xline(ax1, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end

title(ax1, sprintf('Original Local Field Potential Signals: %.3f-%.3f s', t1, t2));
xlim(ax1, [t1, t2]);
ylim(ax1, [ytick_pos(1) - 1, ytick_pos(end) + 1]);
set(ax1, 'YTick', ytick_pos, 'YTickLabel', ytick_labels);
set(ax1, 'Box', 'off');

if n_regions >= 1
    set(ax1, 'XTickLabel', []);
end

%% =========================
%  Plot region spectrograms
%  =========================
for r = 1:n_regions
    ax = nexttile(tl, r + 1);
    ax_all(r + 1) = ax;
    disableDefaultInteractivity(ax);

    imagesc(ax, t_spec, f_spec, spec_seg(:, :, r));
    set(ax, 'YDir', 'normal');

    apply_colormap(ax, cmap_in);

    if ~isempty(clim)
        caxis(ax, clim);
    end

    hold(ax, 'on');
    for i = 1:numel(border_t_raw)
        xline(ax, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end

    title(ax, sprintf('Power Spectrogram (%s)', regions{r}));
    ylabel(ax, 'Frequency (Hz)');
    xlim(ax, [t1, t2]);
    ylim(ax, freq_lim);
    set(ax, 'Box', 'off');

    if r < n_regions
        set(ax, 'XTickLabel', []);
    else
        xlabel(ax, 'Time (s)');
    end

    colorbar(ax);
end

linkaxes(ax_all, 'x');
drawnow;

end


function t = build_global_time_axis(session_lengths, session_dx)
% Build a global time axis by concatenating sessions in time.

n_sessions = numel(session_lengths);
t_cells = cell(n_sessions, 1);

for k = 1:n_sessions
    n = double(session_lengths(k));
    dx = double(session_dx(k));

    t_local = (0:n-1) * dx;

    if k == 1
        t_global = t_local;
    else
        t_prev = t_cells{k-1};
        t_global = t_local + t_prev(end) + dx;
    end

    t_cells{k} = t_global;
end

t = cat(2, t_cells{:});
t = double(t(:));
end


function spec_file = find_regionmean_spectrogram_file(cfg, output_root, pad_mode, pad_sec)
% Find the saved region-mean spectrogram file.

if strcmpi(pad_mode, 'mirror')
    pad_tag = sprintf('_mirrorpad_%gs', pad_sec);
else
    pad_tag = '_nopad';
end
pad_tag = strrep(pad_tag, '.', 'p');

target_name = [cfg.file_stem, pad_tag, '_regionmean_spectrograms_abs.mat'];

search_dirs = { ...
    fullfile(output_root, cfg.file_stem, 'spectrograms'), ...
    fullfile(output_root, 'spectrograms', cfg.file_stem)};

for i = 1:numel(search_dirs)
    f = fullfile(search_dirs{i}, target_name);
    if exist(f, 'file') == 2
        spec_file = f;
        return;
    end
end

candidate_files = {};

for i = 1:numel(search_dirs)
    if exist(search_dirs{i}, 'dir') ~= 7
        continue;
    end

    L = dir(fullfile(search_dirs{i}, [cfg.file_stem, '*_regionmean_spectrograms_abs.mat']));
    for j = 1:numel(L)
        candidate_files{end+1, 1} = fullfile(L(j).folder, L(j).name);
    end
end

candidate_files = unique(candidate_files);

if isempty(candidate_files)
    msg = sprintf(['No saved region-mean spectrogram file was found.\n\n' ...
        'Checked paths:\n%s\n%s\n\n' ...
        'Expected file name:\n%s'], ...
        search_dirs{1}, search_dirs{2}, target_name);
    error(msg);
elseif numel(candidate_files) == 1
    spec_file = candidate_files{1};
else
    msg = sprintf('Multiple candidate files were found:\n');
    for i = 1:numel(candidate_files)
        msg = sprintf('%s  %s\n', msg, candidate_files{i});
    end
    error(msg);
end
end


function apply_colormap(ax, cmap_in)
% Apply colormap from string, function handle, or numeric array.

if ischar(cmap_in) || isstring(cmap_in)
    cmap_name = char(cmap_in);
    colormap(ax, feval(cmap_name, 256));
elseif isa(cmap_in, 'function_handle')
    colormap(ax, cmap_in(256));
elseif isnumeric(cmap_in)
    colormap(ax, cmap_in);
else
    error('Unsupported colormap input.');
end
end