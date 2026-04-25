function O = build_bold_observables(D, params)
% Build BOLD observables from raw, preprocessed, or filtered BOLD data.
%
% params.source:
%   'data'          use D.data
%   'filtered'      use D.filtered_data when available, otherwise D.data
%
% params.mode:
%   'identity'      keep variables as observables
%   'roi_mean'      average variables within each ROI/region
%   'slow_band_power' append canonical BOLD slow-band power observables
%   'svd'           SVD scores across variables

if nargin < 2
    params = struct();
end
params = apply_default_params(params);

X = select_source_data(D, params.source);
if isempty(X)
    error('Selected BOLD source data is empty.');
end

switch lower(params.mode)
    case 'identity'
        obs = double(X);
        observable_info = build_identity_info(D);
        model = struct('type', 'identity');

    case 'roi_mean'
        [obs, observable_info] = build_roi_mean_observables(X, D);
        model = struct('type', 'roi_mean');

    case 'slow_band_power'
        [obs, observable_info, model] = build_slow_band_power_observables(X, D, params);

    case 'svd'
        [obs, observable_info, model] = build_lowdim_observables(X, D, params);

    otherwise
        error('Unsupported params.mode: %s', params.mode);
end

O = struct();
O.data = cast_observables(obs, params.precision);
O.observable_info = observable_info;
O.observable_labels = observable_info.observable_label;
O.params = params;
O.model = model;
O.source_data_kind = params.source;
O.session_ids = D.session_ids;
O.session_lengths = D.session_lengths;
O.session_dx = D.session_dx;
O.session_start_idx = D.session_start_idx;
O.session_end_idx = D.session_end_idx;
O.border_idx = D.border_idx;
O.dx = D.dx;
if isfield(D, 'fs'), O.fs = D.fs; end
if isfield(D, 'dataset_id'), O.dataset_id = D.dataset_id; end
if isfield(D, 'file_stem'), O.file_stem = D.file_stem; end
end

function params = apply_default_params(params)
if ~isfield(params, 'source') || isempty(params.source)
    params.source = 'data';
end
if ~isfield(params, 'mode') || isempty(params.mode)
    params.mode = 'identity';
end
if ~isfield(params, 'n_components') || isempty(params.n_components)
    params.n_components = 50;
end
if ~isfield(params, 'center') || isempty(params.center)
    params.center = true;
end
if ~isfield(params, 'precision') || isempty(params.precision)
    params.precision = 'single';
end
if ~isfield(params, 'slow_bands_hz') || isempty(params.slow_bands_hz)
    params.slow_bands_hz = [0.01 0.027; 0.027 0.073; 0.073 0.198; 0.198 0.25];
end
if ~isfield(params, 'slow_band_names') || isempty(params.slow_band_names)
    params.slow_band_names = {'slow5', 'slow4', 'slow3', 'slow2'};
end
if ~isfield(params, 'filter_order') || isempty(params.filter_order)
    params.filter_order = 4;
end
if ~isfield(params, 'include_raw') || isempty(params.include_raw)
    params.include_raw = true;
end
if ~isfield(params, 'power_eps') || isempty(params.power_eps)
    params.power_eps = 1e-6;
end
if ~isfield(params, 'band_power_transform') || isempty(params.band_power_transform)
    params.band_power_transform = 'power';
end
end

function X = select_source_data(D, source)
switch lower(source)
    case 'data'
        X = D.data;
    case 'filtered'
        if isfield(D, 'filtered_data') && ~isempty(D.filtered_data)
            X = D.filtered_data;
        else
            X = D.data;
            warning('D.filtered_data is not available. Falling back to D.data.');
        end
    otherwise
        error('Unsupported params.source: %s', source);
end
end

function info = build_identity_info(D)
n_obs = size(D.data, 2);
observable_idx = (1:n_obs).';
source = repmat({'bold'}, n_obs, 1);

if isfield(D, 'variable_labels') && numel(D.variable_labels) == n_obs
    observable_label = cellstr(string(D.variable_labels(:)));
else
    observable_label = cellstr(compose("bold_var%04d", observable_idx));
end

info = table(observable_idx, source, observable_label, ...
    'VariableNames', {'observable_idx', 'source', 'observable_label'});

if isfield(D, 'variable_info') && height(D.variable_info) == n_obs
    extra = D.variable_info;
    duplicate_names = intersect(info.Properties.VariableNames, extra.Properties.VariableNames);
    if ~isempty(duplicate_names)
        extra(:, duplicate_names) = [];
    end
    info = [info extra];
end
end

function [obs, info] = build_roi_mean_observables(X, D)
if ~isfield(D, 'variable_info') || ~ismember('region_idx', D.variable_info.Properties.VariableNames)
    obs = double(X);
    info = build_identity_info(D);
    info.source(:) = {'bold_roi_mean'};
    info.observable_label = cellstr(strcat(string(info.observable_label), "_mean"));
    return;
end

region_idx = D.variable_info.region_idx;
region_labels = D.variable_info.region_label;
regions = unique(region_idx(:), 'stable');

obs = zeros(size(X, 1), numel(regions));
observable_idx = (1:numel(regions)).';
source = repmat({'bold_roi_mean'}, numel(regions), 1);
observable_label = cell(numel(regions), 1);
region_id = zeros(numel(regions), 1);
region_label = cell(numel(regions), 1);

for i = 1:numel(regions)
    mask = region_idx == regions(i);
    obs(:, i) = mean(double(X(:, mask)), 2, 'omitnan');
    region_id(i) = regions(i);
    region_label{i} = region_labels{find(mask, 1, 'first')};
    observable_label{i} = sprintf('%s_mean', region_label{i});
end

info = table(observable_idx, source, observable_label, region_id, region_label);
end

function [obs, info, model] = build_slow_band_power_observables(X, D, params)
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
    X_band = filter_slow_bands_one_session(X(idx, :), fs, bands_hz, params.filter_order);
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

info = build_slow_band_power_info(D, n_var, band_names, params.include_raw);
model = struct();
model.type = 'slow_band_power';
model.bands_hz = bands_hz;
model.band_names = band_names;
model.filter_order = params.filter_order;
model.include_raw = params.include_raw;
model.power_eps = params.power_eps;
model.band_power_transform = params.band_power_transform;
end

function X_band = filter_slow_bands_one_session(X, fs, bands_hz, filter_order)
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

function info = build_slow_band_power_info(D, n_var, band_names, include_raw)
base_info = build_identity_info(D);
if height(base_info) ~= n_var
    observable_idx = (1:n_var).';
    source = repmat({'bold'}, n_var, 1);
    observable_label = cellstr(compose("bold_var%04d", observable_idx));
    base_info = table(observable_idx, source, observable_label);
end

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

function [obs, info, model] = build_lowdim_observables(X, D, params)
X = double(X);
if params.center
    mu = mean(X, 1, 'omitnan');
    Xc = X - mu;
else
    mu = zeros(1, size(X, 2));
    Xc = X;
end

n_comp = min([params.n_components, size(Xc, 1) - 1, size(Xc, 2)]);
if n_comp < 1
    error('Not enough samples or variables for low-dimensional observables.');
end

[U, S, V] = svd(Xc, 'econ');
coeff = V(:, 1:n_comp);
score = U(:, 1:n_comp) * S(1:n_comp, 1:n_comp);
latent = diag(S).^2 / max(size(Xc, 1) - 1, 1);
latent = latent(1:n_comp);
explained = 100 * latent / sum(latent);
obs = score;

observable_idx = (1:n_comp).';
source_tag = lower(params.mode);
if isfield(params, 'observable_branch') && ~isempty(params.observable_branch)
    source_tag = lower(char(string(params.observable_branch)));
end
source = repmat({source_tag}, n_comp, 1);
observable_label = cellstr(compose("%s%03d", source_tag, observable_idx));
explained = explained(:);
latent = latent(:);
info = table(observable_idx, source, observable_label, latent, explained);

model = struct();
model.type = lower(params.mode);
model.center = params.center;
model.mu = mu;
model.coeff = coeff;
model.latent = latent;
model.explained = explained;
if isfield(D, 'variable_labels')
    model.input_variable_labels = D.variable_labels;
end
end

function Y = cast_observables(X, precision)
switch lower(precision)
    case 'single'
        Y = single(X);
    case 'double'
        Y = double(X);
    otherwise
        error('params.precision must be ''single'' or ''double''.');
end
end
