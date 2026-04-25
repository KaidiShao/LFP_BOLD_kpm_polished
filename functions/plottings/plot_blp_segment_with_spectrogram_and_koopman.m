function [hfig, plot_info] = plot_blp_segment_with_spectrogram_and_koopman( ...
    base_plot_cache, cmap_in, clim, koopman_input)
% Plot raw BLP traces, region spectrograms, Koopman eigenfunctions, and Koopman residuals.
%
% This extends plot_blp_segment_with_spectrogram.m by adding two extra
% heatmaps for a preloaded time window:
%   1) Koopman eigenfunctions
%   2) Koopman residuals u_t = phi_t - lambda * phi_{t-1}
%
% koopman_input fields (required unless already in workspace):
%   .EDMD_outputs       struct with fields evalues/efuns (or original_sorted.*)
%   .efuns              [T x K] complex eigenfunctions, alternative to .EDMD_outputs
%   .evalues            [K x 1] eigenvalues, alternative to .EDMD_outputs
%   .residual           [T x K] optional precomputed residuals
%   .time_vec           [T x 1] optional explicit Koopman time axis
%   .mode_idx           indices of modes to plot (default 1:max_modes)
%   .max_modes          default 20
%   .use_original_sorted logical, default false
%   .feature            'abs' (default) or 'real'
%   .normalize_scope    'window' (default) or 'global'
%   .first_u_mode       'phi1' (default) or 'zero'
%   .koopman_cmap       colormap for Koopman heatmaps (default parula(256))
%   .koopman_clim       optional clim for eigenfunction heatmap
%   .residual_clim      optional clim for residual heatmap
%   .base_plot_cache    no longer used; pass the prebuilt cache as arg 1

if nargin < 1 || isempty(base_plot_cache)
    error('base_plot_cache must be provided. Build it first with build_blp_plot_window_cache.');
end

if nargin < 2 || isempty(cmap_in)
    cmap_in = parula(256);
end

if nargin < 3
    clim = [];
end

if nargin < 4 || isempty(koopman_input)
    error('koopman_input must be provided.');
end

t1 = double(base_plot_cache.time_range_sec(1));
t2 = double(base_plot_cache.time_range_sec(2));
freq_lim = base_plot_cache.freq_lim;
used_cached_base = true;

t_raw = base_plot_cache.t_raw;
raw_idx = base_plot_cache.raw_idx;
t_seg = base_plot_cache.t_seg;
idx_plot_1 = base_plot_cache.idx_plot_1;
idx_plot_2 = base_plot_cache.idx_plot_2;
R_events = base_plot_cache.R_events;
n_event_bands = base_plot_cache.n_event_bands;
band_colors = base_plot_cache.band_colors;
spec_file = base_plot_cache.spec_file;
regions = base_plot_cache.regions;
spec_seg = base_plot_cache.spec_seg;
t_spec = base_plot_cache.t_spec;
f_spec = base_plot_cache.f_spec;
n_channels = base_plot_cache.n_channels;
x_disp = base_plot_cache.x_disp;
offsets = base_plot_cache.offsets;
ytick_pos = base_plot_cache.ytick_pos;
ytick_labels = base_plot_cache.ytick_labels;
border_t_raw = base_plot_cache.border_t_raw;
show_events = base_plot_cache.show_events;
trace_linewidth = base_plot_cache.plot_settings.trace_linewidth;

% -------------------- koopman data --------------------
if ~isfield(koopman_input, 'koopman_cmap') || isempty(koopman_input.koopman_cmap)
    koopman_input.koopman_cmap = cmap_in;
end
K = local_prepare_koopman_segment(koopman_input, t1, t2, t_raw);

% -------------------- figure --------------------
n_regions = size(spec_seg, 3);
n_tiles = n_regions + 3;
hfig = figure('Color', 'w', 'Position', [100, 100, 1080, 400 + 220 * n_tiles]);
tl = tiledlayout(n_tiles, 1, 'TileSpacing', 'compact', 'Padding', 'compact');
ax_all = gobjects(n_tiles, 1);

% Raw traces
ax1 = nexttile(tl, 1);
ax_all(1) = ax1;
disableDefaultInteractivity(ax1);
hold(ax1, 'on');

for c = 1:n_channels
    plot(ax1, t_seg, x_disp(:, c) + offsets(c), 'k', 'LineWidth', trace_linewidth);
end

if show_events
    for c = 1:n_channels
        for b = 1:n_event_bands
            win = R_events.DetectResults{b, c}.event_win;
            if isempty(win), continue; end
            overlap_mask = win(:, 2) >= idx_plot_1 & win(:, 1) <= idx_plot_2;
            win = win(overlap_mask, :);

            for i = 1:size(win, 1)
                g1 = max(win(i, 1), idx_plot_1);
                g2 = min(win(i, 2), idx_plot_2);
                if g2 < g1, continue; end

                local1 = g1 - idx_plot_1 + 1;
                local2 = g2 - idx_plot_1 + 1;
                plot(ax1, t_seg(local1:local2), x_disp(local1:local2, c) + offsets(c), ...
                    'Color', band_colors(b, :), 'LineWidth', max(trace_linewidth, 1.0));
            end
        end
    end
end

for i = 1:numel(border_t_raw)
    xline(ax1, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end

if show_events
    title(ax1, sprintf('Original Local Field Potential Signals with Detected Events: %.3f-%.3f s', t1, t2));
else
    title(ax1, sprintf('Original Local Field Potential Signals: %.3f-%.3f s', t1, t2));
end
xlim(ax1, [t1, t2]);
ylim(ax1, [ytick_pos(1) - 1, ytick_pos(end) + 1]);
set(ax1, 'YTick', ytick_pos, 'YTickLabel', ytick_labels, 'Box', 'off');
set(ax1, 'XTickLabel', []);

if show_events
    legend_handles = gobjects(n_event_bands, 1);
    legend_labels = cell(n_event_bands, 1);
    for b = 1:n_event_bands
        legend_handles(b) = plot(ax1, nan, nan, 'Color', band_colors(b, :), 'LineWidth', 1.4);
        if isfield(R_events, 'params') && isfield(R_events.params, 'band_labels') && numel(R_events.params.band_labels) >= b
            legend_labels{b} = R_events.params.band_labels{b};
        else
            legend_labels{b} = sprintf('Band %d', b);
        end
    end
    legend(ax1, legend_handles, legend_labels, 'Location', 'northeast', 'Box', 'off');
end

% Spectrograms
for r = 1:n_regions
    ax = nexttile(tl, r + 1);
    ax_all(r + 1) = ax;
    disableDefaultInteractivity(ax);

    imagesc(ax, t_spec, f_spec, spec_seg(:, :, r));
    set(ax, 'YDir', 'normal');
    local_apply_colormap(ax, cmap_in);
    if ~isempty(clim), set(ax, 'CLim', clim); end

    hold(ax, 'on');
    for i = 1:numel(border_t_raw)
        xline(ax, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
    end
    title(ax, sprintf('Power Spectrogram (%s)', regions{r}));
    ylabel(ax, 'Frequency (Hz)');
    xlim(ax, [t1, t2]);
    ylim(ax, freq_lim);
    set(ax, 'Box', 'off', 'XTickLabel', []);
    colorbar(ax);
end

% Koopman eigenfunction heatmap
ax_phi = nexttile(tl, n_regions + 2);
ax_all(n_regions + 2) = ax_phi;
disableDefaultInteractivity(ax_phi);
imagesc(ax_phi, K.t_seg, 1:numel(K.mode_idx), K.efun_plot);
set(ax_phi, 'YDir', 'reverse');
local_apply_colormap(ax_phi, K.koopman_cmap);
if ~isempty(K.koopman_clim), set(ax_phi, 'CLim', K.koopman_clim); end
hold(ax_phi, 'on');
for i = 1:numel(border_t_raw)
    xline(ax_phi, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
title(ax_phi, K.efun_title, 'Interpreter', 'none');
ylabel(ax_phi, 'Mode #');
set(ax_phi, 'YTick', 1:numel(K.mode_idx), 'YTickLabel', string(K.mode_idx), 'Box', 'off', 'XTickLabel', []);
colorbar(ax_phi);

% Koopman residual heatmap
ax_res = nexttile(tl, n_regions + 3);
ax_all(n_regions + 3) = ax_res;
disableDefaultInteractivity(ax_res);
imagesc(ax_res, K.t_seg, 1:numel(K.mode_idx), K.residual_plot);
set(ax_res, 'YDir', 'reverse');
local_apply_colormap(ax_res, K.koopman_cmap);
if ~isempty(K.residual_clim), set(ax_res, 'CLim', K.residual_clim); end
hold(ax_res, 'on');
for i = 1:numel(border_t_raw)
    xline(ax_res, border_t_raw(i), '--', 'Color', [0.6 0.6 0.6], 'LineWidth', 0.8);
end
title(ax_res, K.residual_title, 'Interpreter', 'none');
ylabel(ax_res, 'Mode #');
xlabel(ax_res, 'Time (s)');
set(ax_res, 'YTick', 1:numel(K.mode_idx), 'YTickLabel', string(K.mode_idx), 'Box', 'off');
colorbar(ax_res);

linkaxes(ax_all, 'x');
drawnow;

plot_info = struct();
plot_info.raw_idx = raw_idx;
plot_info.koopman = K;
plot_info.spec_file = spec_file;
plot_info.time_range_sec = [t1, t2];
plot_info.base_plot_cache = base_plot_cache;
plot_info.used_cached_base = used_cached_base;
end


function base_plot_cache = local_prepare_base_plot_cache( ...
    cfg, output_root, t1, t2, freq_lim, event_input, band_colors, ...
    trace_scale, trace_clip, within_gap, between_gap)
% Load raw/event/spectrogram data once for reuse across multiple Koopman views.

show_events = ~isempty(event_input);

D = load_blp_dataset(cfg);
t_raw = local_build_global_time_axis(D.session_lengths, D.session_dx);
t_raw = double(t_raw(:));

raw_idx = find(t_raw >= t1 & t_raw <= t2);
if isempty(raw_idx)
    error('The requested time range does not overlap the raw data.');
end

x_seg = double(D.data(raw_idx, :));
t_seg = double(t_raw(raw_idx));
idx_plot_1 = raw_idx(1);
idx_plot_2 = raw_idx(end);

R_events = [];
n_event_bands = 0;
if show_events
    R_events = local_load_event_results(cfg, output_root, event_input);
    n_event_bands = size(R_events.DetectResults, 1);

    if size(R_events.DetectResults, 2) ~= size(D.data, 2)
        error('Event result channel count (%d) does not match loaded data (%d).', ...
            size(R_events.DetectResults, 2), size(D.data, 2));
    end

    if isempty(band_colors)
        band_colors = local_default_band_colors(n_event_bands);
    end

    if size(band_colors, 1) < n_event_bands || size(band_colors, 2) ~= 3
        error('band_colors must be an N-by-3 array with at least one row per event band.');
    end
end

pad_sec = 20;
pad_mode = 'mirror';
if isfield(cfg, 'spectrogram')
    if isfield(cfg.spectrogram, 'pad_sec'), pad_sec = cfg.spectrogram.pad_sec; end
    if isfield(cfg.spectrogram, 'pad_mode'), pad_mode = cfg.spectrogram.pad_mode; end
end

spec_file = local_find_regionmean_spectrogram_file(cfg, output_root, pad_mode, pad_sec);
meta = load(spec_file, 'freqs', 'regions');
freqs = double(meta.freqs(:));
regions = meta.regions;

freq_idx = find(freqs >= freq_lim(1) & freqs <= freq_lim(2));
if isempty(freq_idx)
    error('No frequencies fall inside the requested freq_lim.');
end

M = matfile(spec_file);
spec_seg = double(M.tmpall_mean_abs(freq_idx, raw_idx, :));
t_spec = double(t_seg(:).');
f_spec = double(freqs(freq_idx));

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

border_t_raw = double(t_raw(D.border_idx));
border_t_raw = border_t_raw(border_t_raw >= t1 & border_t_raw <= t2);

base_plot_cache = struct();
base_plot_cache.D = D;
base_plot_cache.t_raw = t_raw;
base_plot_cache.raw_idx = raw_idx;
base_plot_cache.x_seg = x_seg;
base_plot_cache.t_seg = t_seg;
base_plot_cache.idx_plot_1 = idx_plot_1;
base_plot_cache.idx_plot_2 = idx_plot_2;
base_plot_cache.R_events = R_events;
base_plot_cache.n_event_bands = n_event_bands;
base_plot_cache.band_colors = band_colors;
base_plot_cache.spec_file = spec_file;
base_plot_cache.regions = regions;
base_plot_cache.spec_seg = spec_seg;
base_plot_cache.t_spec = t_spec;
base_plot_cache.f_spec = f_spec;
base_plot_cache.n_channels = n_channels;
base_plot_cache.x_disp = x_disp;
base_plot_cache.offsets = offsets;
base_plot_cache.ytick_pos = ytick_pos;
base_plot_cache.ytick_labels = ytick_labels;
base_plot_cache.border_t_raw = border_t_raw;
base_plot_cache.show_events = show_events;
end


function K = local_prepare_koopman_segment(koopman_input, t1, t2, t_raw)
% Resolve Koopman data, select a window, and prepare heatmaps.
koopman_input = local_set_default(koopman_input, 'mode_idx', []);
koopman_input = local_set_default(koopman_input, 'max_modes', 20);
koopman_input = local_set_default(koopman_input, 'use_original_sorted', false);
koopman_input = local_set_default(koopman_input, 'feature', 'abs');
koopman_input = local_set_default(koopman_input, 'normalize_scope', 'window');
koopman_input = local_set_default(koopman_input, 'first_u_mode', 'phi1');
koopman_input = local_set_default(koopman_input, 'koopman_cmap', parula(256));
koopman_input = local_set_default(koopman_input, 'koopman_clim', []);
koopman_input = local_set_default(koopman_input, 'residual_clim', []);

[efuns_full, evalues_full, source_label] = local_resolve_koopman_data(koopman_input);
K0 = size(efuns_full, 2);

if isempty(koopman_input.mode_idx)
    mode_idx = 1:min(K0, koopman_input.max_modes);
else
    mode_idx = koopman_input.mode_idx(:).';
    mode_idx = mode_idx(mode_idx >= 1 & mode_idx <= K0);
end
if isempty(mode_idx)
    error('No valid Koopman mode indices were selected.');
end

efuns_use = efuns_full(:, mode_idx);
evalues_use = evalues_full(mode_idx);

if isfield(koopman_input, 'time_vec') && ~isempty(koopman_input.time_vec)
    t_kpm = double(koopman_input.time_vec(:));
elseif size(efuns_full, 1) == numel(t_raw)
    t_kpm = double(t_raw(:));
else
    error(['Koopman time axis could not be inferred. Provide koopman_input.time_vec ', ...
        'or ensure size(efuns,1) matches the loaded BLP time axis.']);
end

seg_idx = find(t_kpm >= t1 & t_kpm <= t2);
if isempty(seg_idx)
    error('The requested time range does not overlap the Koopman time axis.');
end

if isfield(koopman_input, 'residual') && ~isempty(koopman_input.residual)
    residual_full = koopman_input.residual(:, mode_idx);
else
    residual_full = local_compute_residual(efuns_use, evalues_use, koopman_input.first_u_mode);
end

efun_plot = local_prepare_mode_heatmap(efuns_use, seg_idx, koopman_input.feature, koopman_input.normalize_scope);
residual_plot = local_prepare_mode_heatmap(residual_full, seg_idx, koopman_input.feature, koopman_input.normalize_scope);

K = struct();
K.source_label = source_label;
K.mode_idx = mode_idx;
K.t_seg = double(t_kpm(seg_idx));
K.feature = koopman_input.feature;
K.normalize_scope = koopman_input.normalize_scope;
K.koopman_cmap = koopman_input.koopman_cmap;
K.koopman_clim = koopman_input.koopman_clim;
K.residual_clim = koopman_input.residual_clim;
K.evalues = evalues_use;
K.efun_plot = efun_plot;
K.residual_plot = residual_plot;
K.efun_title = sprintf('Koopman Eigenfunction Heatmap (%s, %s, %s)', ...
    local_feature_label(koopman_input.feature), koopman_input.normalize_scope, source_label);
K.residual_title = sprintf('Koopman Residual Heatmap (%s, %s, %s)', ...
    local_feature_label(koopman_input.feature), koopman_input.normalize_scope, source_label);
end


function residual = local_compute_residual(efuns, evalues, first_u_mode)
% Compute residual u_t = phi_t - lambda * phi_{t-1}.
[T, K] = size(efuns);
residual = zeros(T, K, 'like', efuns);

switch lower(first_u_mode)
    case 'phi1'
        residual(1, :) = efuns(1, :);
    case 'zero'
        residual(1, :) = 0;
    otherwise
        error('Unknown first_u_mode = %s. Use ''phi1'' or ''zero''.', first_u_mode);
end

if T >= 2
    residual(2:end, :) = efuns(2:end, :) - efuns(1:end-1, :) .* reshape(evalues(:).', 1, []);
end
end


function Xplot = local_prepare_mode_heatmap(Xfull, seg_idx, feature, normalize_scope)
% Prepare a [K x Tseg] heatmap from complex mode series.
if strcmpi(normalize_scope, 'global')
    Xn = local_normalize_mode_matrix(Xfull, feature);
    Xplot = Xn(seg_idx, :).';
else
    Xn = local_normalize_mode_matrix(Xfull(seg_idx, :), feature);
    Xplot = Xn.';
end
end


function Xn = local_normalize_mode_matrix(X, feature)
% Convert to abs/real and normalize each column by its max abs value.
if exist('normalize_efun', 'file') == 2
    Xn = normalize_efun(X, feature);
    return;
end

if strcmpi(feature, 'real')
    Y = real(X);
else
    Y = abs(X);
end

Xn = Y;
for j = 1:size(Y, 2)
    y = Y(:, j);
    mx = max(abs(y));
    if mx > 0
        y = y ./ mx;
    end
    Xn(:, j) = y;
end
end


function [efuns_full, evalues_full, source_label] = local_resolve_koopman_data(koopman_input)
% Resolve eigenfunctions/eigenvalues from input struct.
if isfield(koopman_input, 'EDMD_outputs') && ~isempty(koopman_input.EDMD_outputs)
    E = koopman_input.EDMD_outputs;
    use_original_sorted = isfield(koopman_input, 'use_original_sorted') && koopman_input.use_original_sorted;

    if use_original_sorted && isfield(E, 'original_sorted') && ...
            isfield(E.original_sorted, 'efuns') && isfield(E.original_sorted, 'evalues')
        efuns_full = E.original_sorted.efuns;
        evalues_full = E.original_sorted.evalues(:);
        source_label = 'original_sorted';
    elseif isfield(E, 'efuns') && isfield(E, 'evalues')
        efuns_full = E.efuns;
        evalues_full = E.evalues(:);
        source_label = 'selected';
    else
        error('koopman_input.EDMD_outputs must contain efuns/evalues or original_sorted.*.');
    end
elseif isfield(koopman_input, 'efuns') && isfield(koopman_input, 'evalues')
    efuns_full = koopman_input.efuns;
    evalues_full = koopman_input.evalues(:);
    source_label = 'direct_input';
else
    error('Provide either koopman_input.EDMD_outputs or both koopman_input.efuns and koopman_input.evalues.');
end

if size(efuns_full, 2) ~= numel(evalues_full)
    error('Koopman input mismatch: efuns columns must equal numel(evalues).');
end
end


function label = local_feature_label(feature)
switch lower(feature)
    case 'real'
        label = 'Re';
    otherwise
        label = '|.|';
end
end


function cfg = local_set_default(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end


function t = local_build_global_time_axis(session_lengths, session_dx)
t = build_global_time_axis_from_sessions(session_lengths, session_dx);
end


function spec_file = local_find_regionmean_spectrogram_file(cfg, output_root, pad_mode, pad_sec)
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
    if exist(search_dirs{i}, 'dir') ~= 7, continue; end
    L = dir(fullfile(search_dirs{i}, [cfg.file_stem, '*_regionmean_spectrograms_abs.mat']));
    for j = 1:numel(L)
        candidate_files{end+1, 1} = fullfile(L(j).folder, L(j).name); %#ok<AGROW>
    end
end

candidate_files = unique(candidate_files);
if isempty(candidate_files)
    error('No saved region-mean spectrogram file was found for %s.', cfg.file_stem);
elseif isscalar(candidate_files)
    spec_file = candidate_files{1};
else
    error('Multiple candidate region-mean spectrogram files found for %s.', cfg.file_stem);
end
end


function local_apply_colormap(ax, cmap_in)
if ischar(cmap_in) || isstring(cmap_in)
    colormap(ax, feval(char(cmap_in), 256));
elseif isa(cmap_in, 'function_handle')
    colormap(ax, cmap_in(256));
elseif isnumeric(cmap_in)
    colormap(ax, cmap_in);
else
    error('Unsupported colormap input.');
end
end


function R = local_load_event_results(cfg, output_root, event_input)
if isempty(event_input) || (ischar(event_input) || isstring(event_input)) && strcmpi(string(event_input), "auto")
    event_input = [];
end

if isstruct(event_input)
    if isfield(event_input, 'DetectResults')
        R = event_input;
        return;
    end
    error('event_input struct must contain DetectResults.');
end

if ischar(event_input) || isstring(event_input)
    S = load(char(event_input));
    if ~isfield(S, 'R')
        error('The event-result file does not contain variable R.');
    end
    R = S.R;
    return;
end

search_dir = fullfile(output_root, cfg.file_stem, 'event_detection');
pattern = fullfile(search_dir, [cfg.file_stem, '_bandpass_events_*.mat']);
L = dir(pattern);
if isempty(L)
    error('No event-result file matching %s was found.', pattern);
end

if numel(L) > 1
    is_short_name = false(numel(L), 1);
    for i = 1:numel(L)
        is_short_name(i) = ~isempty(regexp(L(i).name, ...
            ['^', regexptranslate('escape', cfg.file_stem), '_bandpass_events_[0-9]+bands\.mat$'], 'once'));
    end
    if sum(is_short_name) == 1
        L = L(is_short_name);
    else
        error('Multiple event-result files found. Please pass event_input explicitly.');
    end
end

S = load(fullfile(L(1).folder, L(1).name));
if ~isfield(S, 'R')
    error('The event-result file does not contain variable R.');
end
R = S.R;
end


function colors = local_default_band_colors(n_bands)
base = [ ...
    0.0000, 0.4470, 0.7410; ...
    0.8500, 0.3250, 0.0980; ...
    0.4660, 0.6740, 0.1880; ...
    0.4940, 0.1840, 0.5560; ...
    0.9290, 0.6940, 0.1250];

if n_bands <= size(base, 1)
    colors = base(1:n_bands, :);
else
    colors = lines(n_bands);
end
end
