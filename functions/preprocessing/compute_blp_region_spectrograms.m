function S = compute_blp_region_spectrograms(D, cfg, output_root)
% Compute region-mean spectrograms from loaded BLP data.
%
% Inputs
%   D           Output struct from load_blp_dataset
%   cfg         Dataset config struct
%   output_root Root folder for saving spectrogram results
%
% Output
%   S           Struct containing concatenated region-mean spectrograms
%
% Notes
%   - This function computes spectrograms session by session.
%   - Mirror padding is applied within each session if requested.
%   - Padding columns are cropped by index after spectrogram computation.
%   - Region averaging is done within each session to avoid storing
%     full channel-wise spectrograms.

%% =========================
%  Basic settings
%  =========================
if nargin < 3 || isempty(output_root)
    output_root = get_project_processed_root();
end

freq_range = [0.1, 250];
Nfreqs = 251;
pad_sec = 20;
pad_mode = 'mirror';

if isfield(cfg, 'spectrogram')
    if isfield(cfg.spectrogram, 'freq_range')
        freq_range = cfg.spectrogram.freq_range;
    end
    if isfield(cfg.spectrogram, 'nfreqs')
        Nfreqs = cfg.spectrogram.nfreqs;
    end
    if isfield(cfg.spectrogram, 'pad_sec')
        pad_sec = cfg.spectrogram.pad_sec;
    end
    if isfield(cfg.spectrogram, 'pad_mode')
        pad_mode = cfg.spectrogram.pad_mode;
    end
end

%% =========================
%  Determine region labels
%  =========================
if isfield(cfg.channels, 'selected_labels') && ~isempty(cfg.channels.selected_labels)
    regions = cfg.channels.selected_labels;
else
    regions = unique(cfg.channels.sites(D.selected_channels), 'stable');
end

n_regions = numel(regions);

%% =========================
%  Prepare save paths
%  =========================
if ~isfield(cfg, 'file_stem')
    error('cfg.file_stem is required.');
end

save_dir = fullfile(output_root, cfg.file_stem, 'spectrograms');

if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

if strcmpi(pad_mode, 'mirror')
    pad_tag = sprintf('_mirrorpad_%gs', pad_sec);
else
    pad_tag = '_nopad';
end
pad_tag = strrep(pad_tag, '.', 'p');

file_abs = fullfile(save_dir, [cfg.file_stem, pad_tag, '_regionmean_spectrograms_abs.mat']);
file_complex = fullfile(save_dir, [cfg.file_stem, pad_tag, '_regionmean_spectrograms_complex.mat']);

%% =========================
%  Load cached result if available
%  =========================
if exist(file_abs, 'file') == 2 && exist(file_complex, 'file') == 2
    S_abs = load(file_abs);
    S_complex = load(file_complex, 'tmpall_mean_complex');

    S = S_abs;
    S.tmpall_mean_complex = S_complex.tmpall_mean_complex;
    return;
end

%% =========================
%  Initialize metadata
%  =========================
n_sessions = numel(D.session_ids);

mean_abs_cells = cell(n_sessions, 1);
mean_complex_cells = cell(n_sessions, 1);

session_spec_lengths = zeros(n_sessions, 1);
session_timesout = cell(n_sessions, 1);
session_fs = 1 ./ D.session_dx;
session_pad_samples = zeros(n_sessions, 1);
session_pad_sec = zeros(n_sessions, 1);

%% =========================
%  Map selected channels to regions
%  =========================
channel_sites = cfg.channels.sites(D.selected_channels);
region_idx = zeros(1, numel(D.selected_channels));

for r = 1:n_regions
    region_idx(strcmp(channel_sites, regions{r})) = r;
end

if any(region_idx == 0)
    error('Some selected channels do not match the selected region labels.');
end

freqs = [];

%% =========================
%  Loop over sessions
%  =========================
for k = 1:n_sessions
    sid = D.session_ids(k);

    % Extract the current session from concatenated raw data
    idx1 = D.session_start_idx(k);
    idx2 = D.session_end_idx(k);
    x = read_blp_data_slice(D, idx1:idx2);   % time x channel

    dx = D.session_dx(k);
    Fs_session = 1 / dx;
    raw_n = size(x, 1);

    fprintf('[%d/%d] Spectrogram session %04d ...\n', k, n_sessions, sid);

    %% -------------------------
    %  Build padded signal
    %  -------------------------
    switch lower(pad_mode)
        case 'mirror'
            pad_n = round(pad_sec / dx);
            pad_n = min(pad_n, raw_n);

            if pad_n > 0
                x_pad = [flipud(x(1:pad_n, :)); x; flipud(x(end-pad_n+1:end, :))];
            else
                x_pad = x;
            end

        case {'none', 'no', 'off'}
            pad_n = 0;
            x_pad = x;

        otherwise
            error('Unsupported pad_mode: %s', pad_mode);
    end

    session_pad_samples(k) = pad_n;
    session_pad_sec(k) = pad_n * dx;

    %% -------------------------
    %  Initialize region accumulators
    %  -------------------------
    session_mean_complex_sum = [];
    session_mean_abs_sum = [];
    region_counts = zeros(1, n_regions);

    keep_idx = [];
    n_time_kept = [];

    %% -------------------------
    %  Loop over selected channels
    %  -------------------------
    for d = 1:size(x_pad, 2)
        [spec_full, freqs_this, ~] = timefreqMB( ...
            x_pad(:, d), Fs_session, 'freqs', freq_range, 'nfreqs', Nfreqs);

        % This pipeline assumes one spectrogram column per raw sample
        if size(spec_full, 2) ~= size(x_pad, 1)
            error(['timefreqMB returned %d time bins for an input length of %d. ' ...
                   'This pipeline assumes one spectrogram column per sample.'], ...
                   size(spec_full, 2), size(x_pad, 1));
        end

        % Crop padding columns by direct index
        if pad_n > 0
            keep_idx = (pad_n + 1):(pad_n + raw_n);
        else
            keep_idx = 1:raw_n;
        end

        spec = spec_full(:, keep_idx);

        if d == 1
            if isempty(freqs)
                freqs = freqs_this;
            else
                if numel(freqs) ~= numel(freqs_this) || any(abs(freqs(:) - freqs_this(:)) > 1e-10)
                    error('Frequency grid is inconsistent across sessions.');
                end
            end

            n_time_kept = size(spec, 2);

            session_mean_complex_sum = complex( ...
                zeros(size(spec, 1), n_time_kept, n_regions), ...
                zeros(size(spec, 1), n_time_kept, n_regions));

            session_mean_abs_sum = zeros(size(spec, 1), n_time_kept, n_regions);
        end

        r = region_idx(d);

        session_mean_complex_sum(:, :, r) = session_mean_complex_sum(:, :, r) + spec;
        session_mean_abs_sum(:, :, r) = session_mean_abs_sum(:, :, r) + abs(spec);
        region_counts(r) = region_counts(r) + 1;
    end

    %% -------------------------
    %  Average within each region
    %  -------------------------
    session_mean_complex = complex( ...
        nan(size(session_mean_complex_sum)), ...
        nan(size(session_mean_complex_sum)));

    session_mean_abs = nan(size(session_mean_abs_sum));

    for r = 1:n_regions
        if region_counts(r) > 0
            session_mean_complex(:, :, r) = session_mean_complex_sum(:, :, r) / region_counts(r);
            session_mean_abs(:, :, r) = session_mean_abs_sum(:, :, r) / region_counts(r);
        end
    end

    %% -------------------------
    %  Store session result
    %  -------------------------
    mean_complex_cells{k} = session_mean_complex;
    mean_abs_cells{k} = session_mean_abs;

    session_spec_lengths(k) = size(session_mean_abs, 2);
    session_timesout{k} = (0:raw_n-1) * dx;

    fprintf('         done, raw = [%d, %d], pad = %d samples, spec = [%d, %d, %d], dx = %.12g\n', ...
        size(x, 1), size(x, 2), pad_n, ...
        size(session_mean_abs, 1), size(session_mean_abs, 2), size(session_mean_abs, 3), dx);
end

%% =========================
%  Concatenate sessions
%  =========================
tmpall_mean_complex = cat(2, mean_complex_cells{:});
tmpall_mean_abs = cat(2, mean_abs_cells{:});

session_lengths = session_spec_lengths;
session_end_idx = cumsum(session_lengths);
session_start_idx = [1; session_end_idx(1:end-1) + 1];
border_idx = session_end_idx(1:end-1);

%% =========================
%  Build global time axis
%  =========================
timesout = build_global_time_axis_from_sessions(session_lengths, D.session_dx);

%% =========================
%  Collect metadata
%  =========================
session_ids = D.session_ids;
session_raw_lengths = D.session_lengths;
session_dx = D.session_dx;
selected_channels = D.selected_channels;

%% =========================
%  Save results
%  =========================
save(file_abs, ...
    'tmpall_mean_abs', 'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', 'session_timesout', ...
    'pad_mode', 'pad_sec', 'session_pad_samples', 'session_pad_sec', '-v7.3');

save(file_complex, ...
    'tmpall_mean_complex', 'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', 'session_timesout', ...
    'pad_mode', 'pad_sec', 'session_pad_samples', 'session_pad_sec', '-v7.3');

%% =========================
%  Pack output struct
%  =========================
S = struct();
S.tmpall_mean_abs = tmpall_mean_abs;
S.tmpall_mean_complex = tmpall_mean_complex;
S.freqs = freqs;
S.timesout = timesout;
S.regions = regions;

S.session_ids = session_ids;
S.session_lengths = session_lengths;
S.session_raw_lengths = session_raw_lengths;
S.session_start_idx = session_start_idx;
S.session_end_idx = session_end_idx;
S.border_idx = border_idx;

S.selected_channels = selected_channels;
S.session_dx = session_dx;
S.session_fs = session_fs;
S.session_timesout = session_timesout;

S.pad_mode = pad_mode;
S.pad_sec = pad_sec;
S.session_pad_samples = session_pad_samples;
S.session_pad_sec = session_pad_sec;
end
