function P = preprocess_bold_sessions(D, params)
% Preprocess BOLD data session by session.
%
% Supported params fields:
%   .demean          true by default
%   .zscore          false by default
%   .detrend         false by default
%   .notch_hz        [] by default
%   .notch_q         35 by default
%   .bandpass_hz     [] by default, e.g. [0.01 0.1]
%   .save_filtered   true by default

if nargin < 2
    params = struct();
end

params = apply_default_params(params);

if ~isfield(D, 'data') || isempty(D.data)
    error('D.data is empty. Load BOLD data into memory before preprocessing.');
end

P = D;
P.raw_data = [];
P.data = zeros(size(D.data), 'like', double(D.data));

if params.save_filtered
    P.filtered_data = [];
end

for k = 1:numel(D.session_ids)
    idx = D.session_start_idx(k):D.session_end_idx(k);
    x = double(D.data(idx, :));

    if params.demean
        x = x - mean(x, 1, 'omitnan');
    end

    if params.detrend
        x = detrend(x);
    end

    if ~isempty(params.notch_hz)
        x = apply_notch_filter(x, D.session_dx(k), params.notch_hz, params.notch_q);
    end

    if ~isempty(params.bandpass_hz)
        x = apply_bandpass_filter(x, D.session_dx(k), params.bandpass_hz);
    end

    if params.zscore
        mu = mean(x, 1, 'omitnan');
        sig = std(x, 0, 1, 'omitnan');
        sig(~isfinite(sig) | sig == 0) = 1;
        x = (x - mu) ./ sig;
    end

    P.data(idx, :) = x;
end

P.preprocessing_params = params;
P.data_kind = 'preprocessed_bold';
if params.save_filtered
    P.filtered_data = P.data;
end
end

function params = apply_default_params(params)
if ~isfield(params, 'demean'), params.demean = true; end
if ~isfield(params, 'zscore'), params.zscore = false; end
if ~isfield(params, 'detrend'), params.detrend = false; end
if ~isfield(params, 'notch_hz'), params.notch_hz = []; end
if ~isfield(params, 'notch_q'), params.notch_q = 35; end
if ~isfield(params, 'bandpass_hz'), params.bandpass_hz = []; end
if ~isfield(params, 'save_filtered'), params.save_filtered = true; end
end

function y = apply_notch_filter(x, dx, notch_hz, notch_q)
fs = 1 / dx;
y = x;

for i = 1:numel(notch_hz)
    f0 = notch_hz(i);
    if f0 <= 0 || f0 >= fs / 2
        error('notch_hz must be between 0 and Nyquist frequency.');
    end

    bw = f0 / notch_q;
    if exist('designfilt', 'file') == 2
        filt_obj = designfilt('bandstopiir', 'FilterOrder', 2, ...
            'HalfPowerFrequency1', f0 - bw/2, ...
            'HalfPowerFrequency2', f0 + bw/2, ...
            'DesignMethod', 'butter', 'SampleRate', fs);
        y = filtfilt(filt_obj, y);
    else
        [b, a] = butter(2, [(f0 - bw/2), (f0 + bw/2)] / (fs/2), 'stop');
        y = filtfilt(b, a, y);
    end
end
end

function y = apply_bandpass_filter(x, dx, bandpass_hz)
fs = 1 / dx;
if numel(bandpass_hz) ~= 2 || bandpass_hz(1) <= 0 || ...
        bandpass_hz(2) <= bandpass_hz(1) || bandpass_hz(2) >= fs / 2
    error('bandpass_hz must be [low high] within (0, Nyquist).');
end

[b, a] = butter(2, bandpass_hz / (fs/2), 'bandpass');
y = filtfilt(b, a, x);
end
