function E = compute_blp_event_density(cfg, output_root, R, params, source_event_file)
% Compute event-density time series from saved event-detection results.
%
% Inputs
%   cfg         Dataset config struct
%   output_root Root folder for processed outputs
%   R           Loaded event-detection result struct
%   params      Struct with density parameters
%   source_event_file Optional source file path for metadata
%
% Output
%   E           Struct containing density results and metadata

%% =========================
%  Input defaults
%  =========================
if nargin < 2 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 3
    error('R must be provided. Load it first with load_event_results.');
end

if nargin < 4
    params = struct();
end

if nargin < 5
    source_event_file = '';
end

params = apply_default_params(params);

%% =========================
%  Validate event-detection results
%  =========================
if ~isfield(R, 'DetectResults') || isempty(R.DetectResults)
    error('The loaded event result does not contain DetectResults.');
end

if isempty(source_event_file) && isfield(R, 'save_file') && ~isempty(R.save_file)
    source_event_file = R.save_file;
end
source_event_file_signature = build_file_signature(source_event_file);

%% =========================
%  Prepare save path / cache
%  =========================
save_dir = fullfile(output_root, cfg.file_stem, 'event_density');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_tag = build_save_tag(cfg.file_stem, params.bin_sec);
save_file = fullfile(save_dir, [save_tag, '.mat']);

if exist(save_file, 'file') == 2 && ~params.force_recompute
    S = load(save_file);
    if isfield(S, 'E') && can_reuse_saved_result(S.E, params, source_event_file_signature)
        E = S.E;
        return;
    end
end

n_bands = size(R.DetectResults, 1);
n_channels = size(R.DetectResults, 2);

if isfield(R, 'session_lengths') && isfield(R, 'session_dx') && ...
        ~isempty(R.session_lengths) && ~isempty(R.session_dx)
    [t_raw, ~, ~, t_end] = build_global_time_axis_from_sessions(R.session_lengths, R.session_dx);
    dx_ref = median(double(R.session_dx(:)));
else
    [dx_ref, L_total] = infer_dt_and_length(R.DetectResults);
    t_raw = (0:L_total-1)' * dx_ref;
    t_end = t_raw(end) + dx_ref;
end

edges = 0:params.bin_sec:t_end;
if edges(end) < t_end
    edges = [edges, t_end];
end
t_centers = edges(1:end-1) + diff(edges) / 2;
n_bins = numel(t_centers);

density_by_chan = zeros(n_bins, n_channels, n_bands);
smoothed_density_by_chan = zeros(n_bins, n_channels, n_bands);
density_mean = zeros(n_bins, n_bands);
smoothed_density_mean = zeros(n_bins, n_bands);
counts_by_chan = zeros(n_bins, n_channels, n_bands);
total_event_count = zeros(n_channels, n_bands);

if isfield(R, 'params') && isfield(R.params, 'band_labels') && ...
        numel(R.params.band_labels) >= n_bands
    band_labels = R.params.band_labels(:);
else
    band_labels = cell(n_bands, 1);
    for b = 1:n_bands
        band_labels{b} = sprintf('band%d', b);
    end
end

g = build_gaussian_kernel(params.smooth_sigma_sec, params.bin_sec);

%% =========================
%  Density computation
%  =========================
for b = 1:n_bands
    for c = 1:n_channels
        S = R.DetectResults{b, c};

        if isempty(S) || ~isfield(S, 'loc_peak') || isempty(S.loc_peak)
            density_by_chan(:, c, b) = 0;
            smoothed_density_by_chan(:, c, b) = 0;
            counts_by_chan(:, c, b) = 0;
            total_event_count(c, b) = 0;
            continue;
        end

        loc_peak = double(S.loc_peak(:));
        loc_peak = loc_peak(loc_peak >= 1 & loc_peak <= numel(t_raw));

        if isempty(loc_peak)
            density_by_chan(:, c, b) = 0;
            smoothed_density_by_chan(:, c, b) = 0;
            counts_by_chan(:, c, b) = 0;
            total_event_count(c, b) = 0;
            continue;
        end

        t_peaks = double(t_raw(loc_peak));
        counts = histcounts(t_peaks, edges);
        dens = counts(:) ./ diff(edges(:));

        counts_by_chan(:, c, b) = counts(:);
        density_by_chan(:, c, b) = dens;
        smoothed_density_by_chan(:, c, b) = conv(dens, g, 'same');
        total_event_count(c, b) = numel(loc_peak);
    end

    density_mean(:, b) = mean(density_by_chan(:, :, b), 2);
    smoothed_density_mean(:, b) = mean(smoothed_density_by_chan(:, :, b), 2);
end

%% =========================
%  Pack output
%  =========================
E = struct();
E.save_file = save_file;
E.source_event_file = source_event_file;
E.source_event_file_signature = source_event_file_signature;
E.dataset_id = cfg.dataset_id;
E.file_stem = cfg.file_stem;

E.bin_sec = params.bin_sec;
E.bin_hz = 1 / params.bin_sec;
E.smooth_sigma_sec = params.smooth_sigma_sec;
if params.smooth_sigma_sec > 0
    E.smooth_sigma_hz = 1 / params.smooth_sigma_sec;
else
    E.smooth_sigma_hz = [];
end

E.t_edges = edges(:);
E.t_centers = t_centers(:);
E.band_labels = band_labels;
E.selected_channels = get_optional_field(R, 'selected_channels', []);
E.channel_sites = get_optional_field(R, 'channel_sites', {});
E.regions = get_optional_field(R, 'regions', {});
E.session_ids = get_optional_field(R, 'session_ids', []);
E.session_lengths = get_optional_field(R, 'session_lengths', []);
E.session_dx = get_optional_field(R, 'session_dx', []);
E.dx_ref = dx_ref;

E.counts_by_chan = counts_by_chan;
E.density_by_chan = density_by_chan;
E.smoothed_density_by_chan = smoothed_density_by_chan;
E.density_mean = density_mean;
E.smoothed_density_mean = smoothed_density_mean;
E.total_event_count = total_event_count;
E.params = params;

E.bands = struct([]);
for b = 1:n_bands
    E.bands(b).name = band_labels{b};
    E.bands(b).band_index = b;
    E.bands(b).t_edges = edges(:);
    E.bands(b).t_centers = t_centers(:);
    E.bands(b).counts_by_chan = counts_by_chan(:, :, b);
    E.bands(b).density_by_chan = density_by_chan(:, :, b);
    E.bands(b).smoothed_density_by_chan = smoothed_density_by_chan(:, :, b);
    E.bands(b).density_mean = density_mean(:, b);
    E.bands(b).smoothed_density_mean = smoothed_density_mean(:, b);
    E.bands(b).total_event_count = total_event_count(:, b);
end

save(save_file, 'E', '-v7.3');
end


function params = apply_default_params(params)
if ~isfield(params, 'bin_sec')
    params.bin_sec = 2;
end

if ~isfield(params, 'smooth_sigma_sec')
    params.smooth_sigma_sec = params.bin_sec;
end

if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end

if ~isscalar(params.bin_sec) || params.bin_sec <= 0
    error('params.bin_sec must be a positive scalar.');
end

if ~isscalar(params.smooth_sigma_sec) || params.smooth_sigma_sec < 0
    error('params.smooth_sigma_sec must be a nonnegative scalar.');
end
end


function t = build_global_time_axis(session_lengths, session_dx)
t = build_global_time_axis_from_sessions(session_lengths, session_dx);
end


function [dx_ref, L_total] = infer_dt_and_length(DetectResults)
% Infer dt and total length from DetectResults when top-level metadata is absent.

dx_ref = [];
L_total = [];

for i = 1:numel(DetectResults)
    S = DetectResults{i};
    if isempty(S)
        continue;
    end

    if isfield(S, 'dx') && ~isempty(S.dx)
        dx_val = double(S.dx);
        if dx_val > 1
            dx_ref = 1 / dx_val;
        else
            dx_ref = dx_val;
        end
    end

    if isfield(S, 'L_total') && ~isempty(S.L_total)
        L_total = double(S.L_total);
    end

    if ~isempty(dx_ref) && ~isempty(L_total)
        return;
    end
end

if isempty(dx_ref) || isempty(L_total)
    error('Could not infer dt and total length from DetectResults.');
end
end


function g = build_gaussian_kernel(sigma_sec, bin_sec)
% Build a normalized Gaussian kernel in density-bin units.

if isempty(sigma_sec) || sigma_sec <= 0
    g = 1;
    return;
end

sig_bins = sigma_sec / bin_sec;
hw = max(1, ceil(4 * sig_bins));
x = (-hw:hw);
g = exp(-(x .^ 2) / (2 * sig_bins ^ 2));
g = g / sum(g);
end


function value = get_optional_field(S, field_name, default_value)
% Read an optional struct field with a default fallback.

if isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function tag = build_save_tag(file_stem, bin_sec)
% Build save tag for event-density output.

sec_tag = strrep(sprintf('%gs', bin_sec), '.', 'p');
tag = sprintf('%s_event_density_%s', file_stem, sec_tag);
end


function tf = can_reuse_saved_result(E, params, source_event_file_signature)
tf = isstruct(E);
if ~tf
    return;
end

required_fields = {'bin_sec', 'smooth_sigma_sec', 'source_event_file_signature'};
for i = 1:numel(required_fields)
    if ~isfield(E, required_fields{i}) || isempty(E.(required_fields{i}))
        tf = false;
        return;
    end
end

tf = double(E.bin_sec) == double(params.bin_sec) && ...
    double(E.smooth_sigma_sec) == double(params.smooth_sigma_sec);
if ~tf
    return;
end

if ~isfield(source_event_file_signature, 'source_ref') || ...
        isempty(source_event_file_signature.source_ref) || ...
        ~isfield(source_event_file_signature, 'exists') || ...
        ~logical(source_event_file_signature.exists)
    tf = false;
    return;
end

tf = file_signature_matches(E.source_event_file_signature, source_event_file_signature);
end
