function [A, meta] = compute_eigenfunction_activity(X, dt, params, session_info)
%COMPUTE_EIGENFUNCTION_ACTIVITY Convert efun/component series to activity.
%
% Mainline activity transforms are sign/phase invariant.  The envelope
% branch computes a session-aware RMS envelope:
%   A(t,j) = sqrt(movmean(abs(X(t,j)).^2, window_j)).

if nargin < 2
    dt = [];
end
if nargin < 3 || isempty(params)
    params = struct();
end
if nargin < 4
    session_info = struct();
end

params = local_apply_defaults(params);
local_validate_input(X);

[T, K] = size(X);
session_info = local_apply_session_defaults(session_info, T);

activity_transform = lower(char(string(params.lfp_activity_transform)));
envelope_policy = lower(char(string(params.envelope_policy)));
do_envelope = logical(params.envelope_enable) || ...
    contains(activity_transform, 'envelope') || ...
    contains(envelope_policy, 'envelope') || ...
    contains(envelope_policy, 'adaptive');

if do_envelope
    [window_samples, window_sec, tau_sec, window_source, window_status] = ...
        local_resolve_envelope_windows(params, dt, K);
    [A, value_transform_used] = local_compute_windowed_activity( ...
        X, session_info, window_samples, activity_transform, params);
else
    [A, value_transform_used] = local_transform_without_envelope( ...
        X, params.value_transform);
    window_samples = zeros(K, 1);
    window_sec = nan(K, 1);
    tau_sec = nan(K, 1);
    window_source = repmat({'none'}, K, 1);
    window_status = repmat({'not_requested'}, K, 1);
end

meta = struct();
meta.activity_transform = char(params.lfp_activity_transform);
meta.activity_window_policy = char(params.lfp_activity_window_policy);
meta.value_transform_requested = char(params.value_transform);
meta.value_transform_used = char(value_transform_used);
meta.envelope_enable = logical(do_envelope);
meta.envelope_policy = char(params.envelope_policy);
meta.envelope_alpha = double(params.envelope_alpha);
meta.envelope_min_window_sec = double(params.envelope_min_window_sec);
meta.envelope_max_window_sec = double(params.envelope_max_window_sec);
meta.envelope_fallback_window_sec = double(params.envelope_fallback_window_sec);
meta.envelope_window_samples_by_mode = double(window_samples(:));
meta.envelope_window_sec_by_mode = double(window_sec(:));
meta.envelope_tau_sec_by_mode = double(tau_sec(:));
meta.envelope_window_source_by_mode = window_source(:);
meta.envelope_window_status_by_mode = window_status(:);
meta.session_aware = true;
meta.n_sessions = numel(session_info.start_idx);
meta.dt = dt;
end


function params = local_apply_defaults(params)
if ~isfield(params, 'value_transform') || isempty(params.value_transform)
    params.value_transform = 'auto';
end
if ~isfield(params, 'lfp_activity_transform') || isempty(params.lfp_activity_transform)
    params.lfp_activity_transform = local_activity_transform_from_value_transform( ...
        params.value_transform);
end
if ~isfield(params, 'lfp_activity_window_policy') || isempty(params.lfp_activity_window_policy)
    params.lfp_activity_window_policy = 'samplewise_no_envelope';
end
if ~isfield(params, 'envelope_enable') || isempty(params.envelope_enable)
    params.envelope_enable = false;
end
if ~isfield(params, 'envelope_policy') || isempty(params.envelope_policy)
    params.envelope_policy = 'none';
end
if ~isfield(params, 'envelope_alpha') || isempty(params.envelope_alpha)
    params.envelope_alpha = 0.35;
end
if ~isfield(params, 'envelope_min_window_sec') || isempty(params.envelope_min_window_sec)
    params.envelope_min_window_sec = 0.03;
end
if ~isfield(params, 'envelope_max_window_sec') || isempty(params.envelope_max_window_sec)
    params.envelope_max_window_sec = 1.0;
end
if ~isfield(params, 'envelope_fallback_window_sec') || isempty(params.envelope_fallback_window_sec)
    params.envelope_fallback_window_sec = 0.10;
end
if ~isfield(params, 'envelope_window_sec'), params.envelope_window_sec = []; end
if ~isfield(params, 'envelope_window_sec_by_mode'), params.envelope_window_sec_by_mode = []; end
if ~isfield(params, 'envelope_window_samples'), params.envelope_window_samples = []; end
if ~isfield(params, 'envelope_window_samples_by_mode'), params.envelope_window_samples_by_mode = []; end
if ~isfield(params, 'evalues_discrete'), params.evalues_discrete = []; end
if ~isfield(params, 'evalues_continuous'), params.evalues_continuous = []; end
if ~isfield(params, 'evalues_bilinear'), params.evalues_bilinear = []; end
if ~isfield(params, 'activity_output_class') || isempty(params.activity_output_class)
    params.activity_output_class = 'single';
end
end


function local_validate_input(X)
if ~isnumeric(X) || ndims(X) ~= 2 || isempty(X)
    error('X must be a nonempty numeric [time x mode] matrix.');
end
end


function session_info = local_apply_session_defaults(session_info, T)
if ~isstruct(session_info)
    session_info = struct();
end
if ~isfield(session_info, 'start_idx') || isempty(session_info.start_idx) || ...
        ~isfield(session_info, 'end_idx') || isempty(session_info.end_idx)
    session_info.start_idx = 1;
    session_info.end_idx = T;
end
session_info.start_idx = round(double(session_info.start_idx(:)));
session_info.end_idx = round(double(session_info.end_idx(:)));
session_info.start_idx = max(1, min(T, session_info.start_idx));
session_info.end_idx = max(1, min(T, session_info.end_idx));
valid = session_info.end_idx >= session_info.start_idx;
session_info.start_idx = session_info.start_idx(valid);
session_info.end_idx = session_info.end_idx(valid);
if isempty(session_info.start_idx)
    session_info.start_idx = 1;
    session_info.end_idx = T;
end
end


function [A, value_transform_used] = local_transform_without_envelope(X, value_transform)
value_transform = lower(char(value_transform));
switch value_transform
    case 'auto'
        if ~isreal(X)
            A = abs(X);
            value_transform_used = 'abs';
        else
            A = X;
            value_transform_used = 'none';
        end
    case {'none', 'asis', 'as_is'}
        if ~isreal(X)
            error(['value_transform="none" requires real-valued input. ', ...
                'Use "abs", "real", or an envelope transform.']);
        end
        A = X;
        value_transform_used = 'none';
    case {'abs', 'magnitude'}
        A = abs(X);
        value_transform_used = 'abs';
    case 'real'
        A = real(X);
        value_transform_used = 'real';
    otherwise
        error('Unknown value_transform = %s.', value_transform);
end
end


function [A, value_transform_used] = local_compute_windowed_activity( ...
        X, session_info, window_samples, activity_transform, params)
[T, K] = size(X);
output_class = char(string(params.activity_output_class));
switch lower(output_class)
    case 'double'
        A = zeros(T, K);
    otherwise
        A = zeros(T, K, 'single');
end

for k = 1:K
    win_k = max(1, round(double(window_samples(k))));
    y = zeros(T, 1, 'like', A(:, 1));
    for s = 1:numel(session_info.start_idx)
        idx = session_info.start_idx(s):session_info.end_idx(s);
        x_abs = abs(X(idx, k));
        switch activity_transform
            case {'rms_envelope', 'abs_envelope', 'adaptive_rms_envelope', ...
                    'abs_rms_envelope'}
                seg = sqrt(movmean(x_abs .^ 2, win_k, ...
                    'Endpoints', 'shrink'));
                value_transform_used = 'rms_envelope_abs';
            case {'abs_smooth', 'smooth_abs', 'movmean_abs'}
                seg = movmean(x_abs, win_k, 'Endpoints', 'shrink');
                value_transform_used = 'movmean_abs';
            otherwise
                seg = sqrt(movmean(x_abs .^ 2, win_k, ...
                    'Endpoints', 'shrink'));
                value_transform_used = 'rms_envelope_abs';
        end
        y(idx) = cast(seg, 'like', y);
    end
    A(:, k) = y;
end
end


function [window_samples, window_sec, tau_sec, source, status] = ...
        local_resolve_envelope_windows(params, dt, K)
window_samples = [];
window_sec = [];
source = repmat({'unknown'}, K, 1);
status = repmat({'ok'}, K, 1);
tau_sec = local_tau_from_params(params, dt, K);

if ~isempty(params.envelope_window_samples_by_mode)
    window_samples = local_vector_param(params.envelope_window_samples_by_mode, K);
    window_sec = local_samples_to_sec(window_samples, dt);
    source(:) = {'params.envelope_window_samples_by_mode'};
elseif ~isempty(params.envelope_window_samples)
    window_samples = local_vector_param(params.envelope_window_samples, K);
    window_sec = local_samples_to_sec(window_samples, dt);
    source(:) = {'params.envelope_window_samples'};
elseif ~isempty(params.envelope_window_sec_by_mode)
    window_sec = local_vector_param(params.envelope_window_sec_by_mode, K);
    window_samples = local_sec_to_samples(window_sec, dt);
    source(:) = {'params.envelope_window_sec_by_mode'};
elseif ~isempty(params.envelope_window_sec)
    window_sec = local_vector_param(params.envelope_window_sec, K);
    window_samples = local_sec_to_samples(window_sec, dt);
    source(:) = {'params.envelope_window_sec'};
else
    [window_sec, source, status] = local_adaptive_window_sec(params, tau_sec, K);
    window_samples = local_sec_to_samples(window_sec, dt);
end

invalid = ~isfinite(window_samples) | window_samples < 1;
if any(invalid)
    if isempty(dt)
        fallback_samples = max(1, round(double(params.envelope_fallback_window_sec)));
    else
        fallback_samples = max(1, round(double(params.envelope_fallback_window_sec) / dt));
    end
    window_samples(invalid) = fallback_samples;
    if isempty(window_sec)
        window_sec = nan(K, 1);
    end
    if ~isempty(dt)
        window_sec(invalid) = window_samples(invalid) .* dt;
    end
    status(invalid) = {'fallback_samples'};
end

window_samples = max(1, round(double(window_samples(:))));
window_sec = double(window_sec(:));
end


function values = local_vector_param(value_in, K)
values = double(value_in(:));
if isscalar(values)
    values = repmat(values, K, 1);
elseif numel(values) ~= K
    error('Envelope window parameter must be scalar or have one value per mode/component.');
end
end


function samples = local_sec_to_samples(sec, dt)
if isempty(dt)
    samples = nan(size(sec));
else
    samples = max(1, round(double(sec(:)) ./ double(dt)));
end
end


function sec = local_samples_to_sec(samples, dt)
if isempty(dt)
    sec = nan(size(samples));
else
    sec = double(samples(:)) .* double(dt);
end
end


function [window_sec, source, status] = local_adaptive_window_sec(params, tau_sec, K)
alpha = double(params.envelope_alpha);
min_sec = double(params.envelope_min_window_sec);
max_sec = double(params.envelope_max_window_sec);
fallback_sec = double(params.envelope_fallback_window_sec);

window_sec = alpha .* tau_sec(:);
source = repmat({'eigenvalue_timescale_adaptive'}, K, 1);
status = repmat({'ok'}, K, 1);

invalid = ~isfinite(window_sec) | window_sec <= 0;
window_sec(invalid) = fallback_sec;
status(invalid) = {'fallback_invalid_timescale'};

too_short = window_sec < min_sec;
window_sec(too_short) = min_sec;
status(too_short & ~invalid) = {'clamped_min'};

too_long = window_sec > max_sec;
window_sec(too_long) = max_sec;
status(too_long & ~invalid) = {'clamped_max'};
end


function tau_sec = local_tau_from_params(params, dt, K)
tau_sec = nan(K, 1);
lambda_d = local_complex_vector(params.evalues_discrete, K);
if ~isempty(dt) && ~all(isnan(real(lambda_d)))
    mag = abs(lambda_d);
    valid_d = isfinite(mag) & mag > 0 & mag < 1;
    tau_sec(valid_d) = -double(dt) ./ log(mag(valid_d));
    unstable = isfinite(mag) & mag >= 1;
    tau_sec(unstable) = Inf;
end

lambda_c = local_complex_vector(params.evalues_continuous, K);
if all(isnan(real(lambda_c)))
    lambda_c = local_complex_vector(params.evalues_bilinear, K);
end
need_c = ~isfinite(tau_sec) | tau_sec <= 0;
valid_c = isfinite(real(lambda_c)) & real(lambda_c) < 0;
tau_sec(need_c & valid_c) = -1 ./ real(lambda_c(need_c & valid_c));
end


function values = local_complex_vector(values_in, K)
values = complex(nan(K, 1), nan(K, 1));
if isempty(values_in)
    return;
end
candidate = values_in(:);
if numel(candidate) == K
    values = candidate;
end
end


function value = local_activity_transform_from_value_transform(value_transform)
switch lower(char(value_transform))
    case {'abs', 'auto'}
        value = 'abs_magnitude';
    case 'real'
        value = 'signed_real';
    otherwise
        value = 'feature_as_is';
end
end
