function S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, opts)
% Compute region-mean spectrograms and stream large arrays directly to disk.
%
% Inputs
%   D           Output struct from load_blp_dataset
%   cfg         Dataset config struct
%   output_root Root folder for saving spectrogram results
%   opts        Optional struct:
%                 .save_precision  'single' (default) | 'double'
%                 .return_data     false (default) | true
%                 .force_recompute false (default) | true
%
% Output
%   S           Struct containing saved-file metadata and, optionally,
%               the full concatenated spectrogram arrays
%
% Notes
%   - This version is intended for large datasets where concatenating all
%     sessions in memory would exceed MATLAB's available RAM.
%   - The saved MAT files keep the same variable names used by the rest of
%     the pipeline: tmpall_mean_abs and tmpall_mean_complex.

%% =========================
%  Defaults
%  =========================
if nargin < 3 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 4
    opts = struct();
end

freq_range = [0.1, 250];
Nfreqs = 251;
pad_sec = 20;
pad_mode = 'mirror';

if ~isfield(opts, 'save_precision') || isempty(opts.save_precision)
    opts.save_precision = 'single';
end

if ~isfield(opts, 'return_data') || isempty(opts.return_data)
    opts.return_data = false;
end

if ~isfield(opts, 'force_recompute') || isempty(opts.force_recompute)
    opts.force_recompute = false;
end

save_precision = lower(char(string(opts.save_precision)));
return_data = logical(opts.return_data);
force_recompute = logical(opts.force_recompute);

if ~ismember(save_precision, {'single', 'double'})
    error('opts.save_precision must be ''single'' or ''double''.');
end

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
%  Cached result
%  =========================
if ~force_recompute && exist(file_abs, 'file') == 2 && exist(file_complex, 'file') == 2
    if local_has_complete_saved_result(file_abs, file_complex)
        S = local_load_saved_result(file_abs, file_complex, return_data);
        return;
    end

    warning(['Existing streamed spectrogram files are incomplete or corrupted. ' ...
        'Deleting them and recomputing:\n  %s\n  %s'], file_abs, file_complex);
end

if exist(file_abs, 'file') == 2
    delete(file_abs);
end

if exist(file_complex, 'file') == 2
    delete(file_complex);
end

%% =========================
%  Metadata allocation
%  =========================
n_sessions = numel(D.session_ids);
session_spec_lengths = zeros(n_sessions, 1);
session_timesout = cell(n_sessions, 1);
session_fs = 1 ./ D.session_dx;
session_pad_samples = zeros(n_sessions, 1);
session_pad_sec = zeros(n_sessions, 1);

channel_sites = cfg.channels.sites(D.selected_channels);
region_idx = zeros(1, numel(D.selected_channels));

for r = 1:n_regions
    region_idx(strcmp(channel_sites, regions{r})) = r;
end

if any(region_idx == 0)
    error('Some selected channels do not match the selected region labels.');
end

total_spec_length = sum(double(D.session_lengths));
freqs = [];

Mabs = [];
Mcomplex = [];

%% =========================
%  Session loop
%  =========================
for k = 1:n_sessions
    sid = D.session_ids(k);
    idx1 = D.session_start_idx(k);
    idx2 = D.session_end_idx(k);
    x = read_blp_data_slice(D, idx1:idx2);   % time x channel

    dx = D.session_dx(k);
    Fs_session = 1 / dx;
    raw_n = size(x, 1);

    fprintf('[%d/%d] Spectrogram session %04d ...\n', k, n_sessions, sid);

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

    region_counts = zeros(1, n_regions);
    if pad_n > 0
        keep_idx = (pad_n + 1):(pad_n + raw_n);
    else
        keep_idx = 1:raw_n;
    end
    n_time_kept = numel(keep_idx);

    zero_block = zeros(Nfreqs, n_time_kept, n_regions, save_precision);
    session_mean_complex_sum = complex(zero_block, zero_block);
    session_mean_abs_sum = zeros(Nfreqs, n_time_kept, n_regions, save_precision);

    for d = 1:size(x_pad, 2)
        [spec_full, freqs_this, ~] = timefreqMB( ...
            x_pad(:, d), Fs_session, 'freqs', freq_range, 'nfreqs', Nfreqs);

        if size(spec_full, 2) ~= size(x_pad, 1)
            error(['timefreqMB returned %d time bins for an input length of %d. ' ...
                'This pipeline assumes one spectrogram column per sample.'], ...
                size(spec_full, 2), size(x_pad, 1));
        end

        spec = spec_full(:, keep_idx);

        if d == 1
            if isempty(freqs)
                freqs = freqs_this(:);
                [Mabs, Mcomplex] = local_preallocate_matfiles( ...
                    file_abs, file_complex, numel(freqs), total_spec_length, n_regions, save_precision);
            else
                if numel(freqs) ~= numel(freqs_this) || any(abs(freqs(:) - freqs_this(:)) > 1e-10)
                    error('Frequency grid is inconsistent across sessions.');
                end
            end
        end

        r = region_idx(d);
        spec_cast = cast(spec, save_precision);
        session_mean_complex_sum(:, :, r) = session_mean_complex_sum(:, :, r) + spec_cast;
        session_mean_abs_sum(:, :, r) = session_mean_abs_sum(:, :, r) + abs(spec_cast);
        region_counts(r) = region_counts(r) + 1;
    end

    session_mean_complex = complex( ...
        nan(size(session_mean_complex_sum, 1), n_time_kept, n_regions, save_precision), ...
        nan(size(session_mean_complex_sum, 1), n_time_kept, n_regions, save_precision));
    session_mean_abs = nan(size(session_mean_abs_sum, 1), n_time_kept, n_regions, save_precision);

    for r = 1:n_regions
        if region_counts(r) > 0
            denom = cast(region_counts(r), save_precision);
            session_mean_complex(:, :, r) = session_mean_complex_sum(:, :, r) ./ denom;
            session_mean_abs(:, :, r) = session_mean_abs_sum(:, :, r) ./ denom;
        end
    end

    Mabs.tmpall_mean_abs(:, idx1:idx2, :) = session_mean_abs;
    Mcomplex.tmpall_mean_complex(:, idx1:idx2, :) = session_mean_complex;

    session_spec_lengths(k) = size(session_mean_abs, 2);
    session_timesout{k} = (0:raw_n-1) * dx;

    fprintf('         done, raw = [%d, %d], pad = %d samples, spec = [%d, %d, %d], dx = %.12g\n', ...
        size(x, 1), size(x, 2), pad_n, ...
        size(session_mean_abs, 1), size(session_mean_abs, 2), size(session_mean_abs, 3), dx);

    clear x x_pad spec_full spec spec_cast session_mean_complex_sum session_mean_abs_sum
    clear session_mean_complex session_mean_abs zero_block
end

%% =========================
%  Global metadata
%  =========================
session_lengths = session_spec_lengths;
session_end_idx = cumsum(session_lengths);
session_start_idx = [1; session_end_idx(1:end-1) + 1];
border_idx = session_end_idx(1:end-1);

timesout = build_global_time_axis_from_sessions(session_lengths, D.session_dx);

session_ids = D.session_ids;
session_raw_lengths = D.session_lengths;
session_dx = D.session_dx;
selected_channels = D.selected_channels;

tmpall_mean_abs_size = [numel(freqs), sum(session_lengths), n_regions];
tmpall_mean_complex_size = tmpall_mean_abs_size;

save(file_abs, ...
    'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', 'session_timesout', ...
    'pad_mode', 'pad_sec', 'session_pad_samples', 'session_pad_sec', ...
    'tmpall_mean_abs_size', 'save_precision', '-append');

save(file_complex, ...
    'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', 'session_timesout', ...
    'pad_mode', 'pad_sec', 'session_pad_samples', 'session_pad_sec', ...
    'tmpall_mean_complex_size', 'save_precision', '-append');

S = struct();
S.abs_file = file_abs;
S.complex_file = file_complex;
S.tmpall_mean_abs_size = tmpall_mean_abs_size;
S.tmpall_mean_complex_size = tmpall_mean_complex_size;
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
S.save_precision = save_precision;
S.data_in_memory = false;
S.tmpall_mean_abs = [];
S.tmpall_mean_complex = [];

if return_data
    S_abs_full = load(file_abs, 'tmpall_mean_abs');
    S_complex_full = load(file_complex, 'tmpall_mean_complex');
    S.tmpall_mean_abs = S_abs_full.tmpall_mean_abs;
    S.tmpall_mean_complex = S_complex_full.tmpall_mean_complex;
    S.data_in_memory = true;
end
end


function [Mabs, Mcomplex] = local_preallocate_matfiles(file_abs, file_complex, n_freq, n_time, n_regions, save_precision)
% Create writable MAT files and preallocate the large spectrogram arrays.

Mabs = matfile(file_abs, 'Writable', true);
Mcomplex = matfile(file_complex, 'Writable', true);

if strcmp(save_precision, 'single')
    Mabs.tmpall_mean_abs(n_freq, n_time, n_regions) = single(0);
    Mcomplex.tmpall_mean_complex(n_freq, n_time, n_regions) = complex(single(0), single(0));
else
    Mabs.tmpall_mean_abs(n_freq, n_time, n_regions) = double(0);
    Mcomplex.tmpall_mean_complex(n_freq, n_time, n_regions) = complex(double(0), double(0));
end
end


function S = local_load_saved_result(file_abs, file_complex, return_data)
% Load saved streamed spectrogram metadata, optionally including arrays.

meta_fields = { ...
    'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', 'session_timesout', ...
    'pad_mode', 'pad_sec', 'session_pad_samples', 'session_pad_sec', ...
    'tmpall_mean_abs_size', 'save_precision'};

S_abs = load(file_abs, meta_fields{:});
S = S_abs;

S.abs_file = file_abs;
S.complex_file = file_complex;
S.tmpall_mean_complex_size = [];

info_complex = whos('-file', file_complex, 'tmpall_mean_complex');
if ~isempty(info_complex)
    S.tmpall_mean_complex_size = info_complex.size;
end

S.tmpall_mean_abs = [];
S.tmpall_mean_complex = [];
S.data_in_memory = false;

if return_data
    S_abs_full = load(file_abs, 'tmpall_mean_abs');
    S_complex_full = load(file_complex, 'tmpall_mean_complex');
    S.tmpall_mean_abs = S_abs_full.tmpall_mean_abs;
    S.tmpall_mean_complex = S_complex_full.tmpall_mean_complex;
    S.tmpall_mean_abs_size = size(S.tmpall_mean_abs);
    S.tmpall_mean_complex_size = size(S.tmpall_mean_complex);
    S.data_in_memory = true;
end
end


function tf = local_has_complete_saved_result(file_abs, file_complex)
% Check whether both streamed MAT files contain the expected variables.

tf = false;

try
    vars_abs = who('-file', file_abs);
    vars_complex = who('-file', file_complex);
catch
    return;
end

required_abs = { ...
    'tmpall_mean_abs', 'freqs', 'timesout', 'regions', ...
    'session_ids', 'session_lengths', 'session_raw_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'selected_channels', 'session_dx', 'session_fs', 'session_timesout', ...
    'pad_mode', 'pad_sec', 'session_pad_samples', 'session_pad_sec', ...
    'tmpall_mean_abs_size'};

required_complex = {'tmpall_mean_complex'};

tf = all(ismember(required_abs, vars_abs)) && all(ismember(required_complex, vars_complex));
end
