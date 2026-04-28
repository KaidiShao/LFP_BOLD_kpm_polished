function [obs, info, model] = build_bold_slow_band_power_observables(X, D, params)
% Build raw-plus-slow-band BOLD observables session by session.

bands_hz = params.slow_bands_hz;
band_names = cellstr(string(params.slow_band_names(:)));
if size(bands_hz, 2) ~= 2
    error('params.slow_bands_hz must be an nBands x 2 matrix.');
end
if numel(band_names) ~= size(bands_hz, 1)
    error('params.slow_band_names must match the number of rows in params.slow_bands_hz.');
end

X = double(X);
[n_time, n_var] = size(X);
n_bands = size(bands_hz, 1);
band_power = zeros(n_time, n_var, n_bands);

for i_sess = 1:numel(D.session_lengths)
    idx = D.session_start_idx(i_sess):D.session_end_idx(i_sess);
    if numel(idx) < 2
        continue;
    end

    fs = 1 ./ D.session_dx(i_sess);
    X_band = local_filter_slow_bands_one_session(X(idx, :), fs, bands_hz, params.filter_order);
    band_power(idx, :, :) = X_band .^ 2 + params.power_eps;
end

switch lower(params.band_power_transform)
    case 'power'
        band_features = band_power;
    case 'log_power'
        band_features = log(band_power);
    otherwise
        error('params.band_power_transform must be ''power'' or ''log_power''.');
end

band_features_2d = reshape(band_features, n_time, n_var * n_bands);
if params.include_raw
    obs = [X, band_features_2d];
else
    obs = band_features_2d;
end

info = local_build_slow_band_power_info(D, n_var, band_names, params.include_raw);
model = struct();
model.type = 'slow_band_power';
model.bands_hz = bands_hz;
model.band_names = band_names;
model.filter_order = params.filter_order;
model.include_raw = params.include_raw;
model.power_eps = params.power_eps;
model.band_power_transform = params.band_power_transform;
end

function X_band = local_filter_slow_bands_one_session(X, fs, bands_hz, filter_order)
[n_time, n_var] = size(X);
n_bands = size(bands_hz, 1);
X_band = zeros(n_time, n_var, n_bands);
nyq = fs / 2;

X_pad = [flipud(X); X; flipud(X)];
idx_mid = (n_time + 1):(2 * n_time);

for i_band = 1:n_bands
    f1 = bands_hz(i_band, 1);
    f2 = bands_hz(i_band, 2);
    if f2 >= nyq
        f2 = 0.99 * nyq;
    end
    if f1 <= 0 || f2 <= f1
        error('Invalid slow band edges after Nyquist clipping: [%g %g] Hz.', f1, f2);
    end

    [b, a] = butter(filter_order, [f1 f2] / nyq, 'bandpass');
    X_filt_pad = filtfilt(b, a, X_pad);
    X_band(:, :, i_band) = X_filt_pad(idx_mid, :);
end
end

function info = local_build_slow_band_power_info(D, n_var, band_names, include_raw)
base_info = local_build_base_observable_info(D, n_var);

info = table();
if include_raw
    raw_info = base_info;
    raw_info.observable_idx = (1:height(raw_info)).';
    raw_info.source(:) = {'bold_raw'};
    info = raw_info;
end

for i_band = 1:numel(band_names)
    band_info = base_info;
    band_info.observable_idx = (height(info) + (1:n_var)).';
    band_info.source(:) = {['bold_' band_names{i_band} '_power']};
    band_info.observable_label = cellstr(strcat(string(band_info.observable_label), ...
        "_", band_names{i_band}, "_power"));
    if isempty(info)
        info = band_info;
    else
        info = [info; band_info]; %#ok<AGROW>
    end
end
end

function base_info = local_build_base_observable_info(D, n_var)
observable_idx = (1:n_var).';
source = repmat({'bold'}, n_var, 1);

if isfield(D, 'variable_labels') && numel(D.variable_labels) == n_var
    observable_label = cellstr(string(D.variable_labels(:)));
else
    observable_label = cellstr(compose("bold_var%04d", observable_idx));
end

base_info = table(observable_idx, source, observable_label, ...
    'VariableNames', {'observable_idx', 'source', 'observable_label'});

if isfield(D, 'variable_info') && height(D.variable_info) == n_var
    extra = D.variable_info;
    duplicate_names = intersect(base_info.Properties.VariableNames, extra.Properties.VariableNames);
    if ~isempty(duplicate_names)
        extra(:, duplicate_names) = [];
    end
    base_info = [base_info extra];
end
end
