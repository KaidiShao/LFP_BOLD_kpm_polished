function plot_data = prepare_blp_plot_data(cfg, output_root, prep_cfg)
%PREPARE_BLP_PLOT_DATA Load reusable data once for BLP plotting scripts.

if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 3
    prep_cfg = struct();
end

prep_cfg = local_apply_default_prep_cfg(prep_cfg);
plot_settings = local_resolve_plot_settings(cfg);

D = io_raw.load_blp_dataset(cfg);
t_raw = io_utils.build_global_time_axis_from_sessions(D.session_lengths, D.session_dx);
t_raw = double(t_raw(:));

[R_events, n_event_bands, band_colors] = local_load_event_support(cfg, output_root, prep_cfg, D);
C_consensus = local_load_consensus_support(cfg, output_root, prep_cfg);
[spec_file, spec_freqs, spec_regions, spec_global_clim] = local_load_spectrogram_support(cfg, output_root, prep_cfg);

plot_data = struct();
plot_data.cfg = cfg;
plot_data.output_root = output_root;
plot_data.prep_cfg = prep_cfg;
plot_data.plot_settings = plot_settings;

plot_data.D = D;
plot_data.t_raw = t_raw;

plot_data.show_events = logical(prep_cfg.show_events);
plot_data.R_events = R_events;
plot_data.n_event_bands = n_event_bands;
plot_data.band_colors = band_colors;

plot_data.show_consensus = logical(prep_cfg.show_consensus);
plot_data.C_consensus = C_consensus;
plot_data.consensus_colors = prep_cfg.consensus_colors;
plot_data.consensus_face_alpha = prep_cfg.consensus_face_alpha;

plot_data.include_spectrogram = logical(prep_cfg.include_spectrogram);
plot_data.spec_file = spec_file;
plot_data.spec_freqs = spec_freqs;
plot_data.spec_regions = spec_regions;
plot_data.spec_global_clim = spec_global_clim;
end


function prep_cfg = local_apply_default_prep_cfg(prep_cfg)
if ~isfield(prep_cfg, 'show_events')
    prep_cfg.show_events = false;
end

if ~isfield(prep_cfg, 'event_input')
    prep_cfg.event_input = [];
end

if ~isfield(prep_cfg, 'band_colors')
    prep_cfg.band_colors = [];
end

if ~isfield(prep_cfg, 'include_spectrogram')
    prep_cfg.include_spectrogram = true;
end

if ~isfield(prep_cfg, 'show_consensus')
    prep_cfg.show_consensus = false;
end

if ~isfield(prep_cfg, 'consensus_input')
    prep_cfg.consensus_input = [];
end

if ~isfield(prep_cfg, 'consensus_colors') || isempty(prep_cfg.consensus_colors)
    prep_cfg.consensus_colors = [ ...
        0.86, 0.74, 0.42; ... % theta
        0.45, 0.64, 0.90; ... % gamma
        0.86, 0.55, 0.75; ... % ripple
        0.63, 0.53, 0.88; ... % theta-gamma
        0.78, 0.45, 0.60];    % sharp-wave-ripple
end

if ~isfield(prep_cfg, 'consensus_face_alpha') || isempty(prep_cfg.consensus_face_alpha)
    prep_cfg.consensus_face_alpha = 0.32;
end
end


function plot_settings = local_resolve_plot_settings(cfg)
plot_settings = struct();
plot_settings.trace_scale = 0.18;
plot_settings.trace_clip = 4;
plot_settings.within_gap = 1.4;
plot_settings.between_gap = 2.2;
plot_settings.trace_linewidth = 0.4;
plot_settings.event_linewidth = 1.1;
plot_settings.border_color = [0.6 0.6 0.6];

if isfield(cfg, 'plot')
    if isfield(cfg.plot, 'trace_scale')
        plot_settings.trace_scale = cfg.plot.trace_scale;
    end
    if isfield(cfg.plot, 'trace_clip')
        plot_settings.trace_clip = cfg.plot.trace_clip;
    end
    if isfield(cfg.plot, 'within_gap')
        plot_settings.within_gap = cfg.plot.within_gap;
    end
    if isfield(cfg.plot, 'between_gap')
        plot_settings.between_gap = cfg.plot.between_gap;
    end
    if isfield(cfg.plot, 'trace_linewidth')
        plot_settings.trace_linewidth = cfg.plot.trace_linewidth;
    end
    if isfield(cfg.plot, 'event_linewidth')
        plot_settings.event_linewidth = cfg.plot.event_linewidth;
    end
end
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


function [R_events, n_event_bands, band_colors] = local_load_event_support(cfg, output_root, prep_cfg, D)
R_events = [];
n_event_bands = 0;
band_colors = prep_cfg.band_colors;

if ~prep_cfg.show_events
    return;
end

event_input = prep_cfg.event_input;
if isempty(event_input) || (ischar(event_input) || isstring(event_input)) && strcmpi(string(event_input), "auto")
    event_input = [];
end

[R_events, ~] = io_results.load_event_results(cfg, output_root, event_input);
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


function C_consensus = local_load_consensus_support(cfg, output_root, prep_cfg)
C_consensus = [];

if ~prep_cfg.show_consensus
    return;
end

consensus_input = prep_cfg.consensus_input;
if isempty(consensus_input) || (ischar(consensus_input) || isstring(consensus_input)) && strcmpi(string(consensus_input), "auto")
    consensus_input = [];
end

[C_consensus, ~] = io_results.load_consensus_state_results(cfg, output_root, consensus_input);
end


function [spec_file, spec_freqs, spec_regions, spec_global_clim] = local_load_spectrogram_support(cfg, output_root, prep_cfg)
spec_file = '';
spec_freqs = [];
spec_regions = {};
spec_global_clim = [];

if ~prep_cfg.include_spectrogram
    return;
end

spec_support = require_saved_blp_spectrogram_support(cfg, output_root);
spec_file = spec_support.abs_file;

meta = load(spec_file, 'freqs', 'regions');
spec_freqs = double(meta.freqs(:));
spec_regions = meta.regions;
spec_global_clim = local_get_or_compute_global_spectrogram_clim(spec_file, numel(spec_regions));
end


function spec_global_clim = local_get_or_compute_global_spectrogram_clim(spec_file, n_regions)
cache_file = [spec_file(1:end-4), '_global_clim.mat'];

if exist(cache_file, 'file') == 2
    S = load(cache_file, 'spec_global_clim');
    if isfield(S, 'spec_global_clim') && isnumeric(S.spec_global_clim) ...
            && size(S.spec_global_clim, 2) == 2 && size(S.spec_global_clim, 1) == n_regions
        spec_global_clim = double(S.spec_global_clim);
        return;
    end
end

M = matfile(spec_file);
var_info = whos('-file', spec_file, 'tmpall_mean_abs');
if isempty(var_info)
    error('tmpall_mean_abs was not found in spectrogram file: %s', spec_file);
end

sz = var_info.size;
if numel(sz) ~= 3 || sz(3) ~= n_regions
    error('Unexpected tmpall_mean_abs size in %s.', spec_file);
end

chunk_size = 100000;
spec_global_clim = nan(n_regions, 2);

for r = 1:n_regions
    region_min = inf;
    region_max = -inf;

    for idx1 = 1:chunk_size:sz(2)
        idx2 = min(idx1 + chunk_size - 1, sz(2));
        block = double(M.tmpall_mean_abs(:, idx1:idx2, r));
        finite_mask = isfinite(block);
        if ~any(finite_mask(:))
            continue;
        end

        block_vals = block(finite_mask);
        region_min = min(region_min, min(block_vals));
        region_max = max(region_max, max(block_vals));
    end

    if ~isfinite(region_min) || ~isfinite(region_max)
        region_min = 0;
        region_max = 1;
    elseif region_min == region_max
        region_max = region_min + eps(region_min + 1);
    end

    spec_global_clim(r, :) = [region_min, region_max];
end

save(cache_file, 'spec_global_clim', 'spec_file', 'chunk_size', '-v7');
end
