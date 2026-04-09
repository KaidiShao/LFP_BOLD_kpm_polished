function fit = fit_empirical_eigenvalues_from_efuns(efuns_raw, lambda_ref, cfg)
%FIT_EMPIRICAL_EIGENVALUES_FROM_EFUNS Estimate empirical eigenvalues from eigenfunction timescales.
%
%   fit = fit_empirical_eigenvalues_from_efuns(efuns_raw, lambda_ref, cfg)
%
% Inputs
%   efuns_raw   [T x K] complex eigenfunctions (time x mode)
%   lambda_ref  [K x 1] reference eigenvalues used to preserve phase if desired
%   cfg         optional struct
%     .dt                    sampling interval (default 1)
%     .lambdaType            'discrete' (default) or 'continuous'
%     .feat                  'abs' (default) or 'real'
%     .maxLag                max lag in samples for ACF (default min(200,T-1))
%     .fit_env_start         default 0.9
%     .fit_env_end           default 0.2
%     .fit_env_floor         default 1e-3
%     .fit_min_points        default 8
%     .fit_fallback          'efold' (default) or 'nan'
%     .lambda_phase_mode     'preserve_edmd_angle' (default) or 'real_positive'
%     .complex_fit_mode      'ls_one_step' (default)
%
% Output
%   fit       struct with fields:
%     .acf_lags
%     .acf_mat
%     .tau_emp
%     .kappa_emp
%     .lambda_ref_discrete
%     .lambda_emp_mag
%     .lambda_emp_discrete
%     .lambda_emp_input
%     .lambda_emp_discrete_safe
%     .lambda_emp_input_safe
%     .lambda_emp_complex_discrete
%     .lambda_emp_complex_input
%     .lambda_emp_complex_discrete_safe
%     .lambda_emp_complex_input_safe
%     .lambda_emp_complex_tau
%     .lambda_emp_complex_kappa
%     .lambda_emp_complex_valid
%     .used_reference_fallback
%     .fit_info

if nargin < 3 || isempty(cfg)
    cfg = struct();
end

cfg = local_set_default(cfg, 'dt', 1);
cfg = local_set_default(cfg, 'lambdaType', 'discrete');
cfg = local_set_default(cfg, 'feat', 'abs');
cfg = local_set_default(cfg, 'maxLag', []);
cfg = local_set_default(cfg, 'fit_env_start', 0.9);
cfg = local_set_default(cfg, 'fit_env_end', 0.2);
cfg = local_set_default(cfg, 'fit_env_floor', 1e-3);
cfg = local_set_default(cfg, 'fit_min_points', 8);
cfg = local_set_default(cfg, 'fit_fallback', 'efold');
cfg = local_set_default(cfg, 'lambda_phase_mode', 'preserve_edmd_angle');
cfg = local_set_default(cfg, 'complex_fit_mode', 'ls_one_step');

[T, K0] = size(efuns_raw);
K = min(K0, numel(lambda_ref));
efuns_raw = efuns_raw(:, 1:K);
lambda_ref = lambda_ref(1:K);

if isempty(cfg.maxLag)
    maxLag = min(200, T - 1);
else
    maxLag = min(cfg.maxLag, T - 1);
end

Phi = local_normalize_efuns(efuns_raw, cfg.feat);
[acf_lags, acf_mat, tau_emp, fit_info] = local_acf_timescale_from_envelope(Phi, maxLag, cfg.dt, cfg);

kappa_emp = 1 ./ tau_emp;
kappa_emp(~isfinite(kappa_emp) | kappa_emp < 0) = 0;

lambda_ref_discrete = local_to_discrete_lambda(lambda_ref, cfg);
lambda_emp_mag = exp(-cfg.dt .* kappa_emp);

switch lower(cfg.lambda_phase_mode)
    case 'preserve_edmd_angle'
        phase = angle(lambda_ref_discrete);
    case 'real_positive'
        phase = zeros(size(lambda_ref_discrete));
    otherwise
        error('Unknown cfg.lambda_phase_mode = %s.', cfg.lambda_phase_mode);
end

lambda_emp_discrete = lambda_emp_mag .* exp(1i * phase);
invalid = ~isfinite(lambda_emp_discrete);
lambda_emp_discrete_safe = lambda_emp_discrete;
lambda_emp_discrete_safe(invalid) = lambda_ref_discrete(invalid);

lambda_emp_input = local_from_discrete_lambda(lambda_emp_discrete, cfg);
lambda_emp_input_safe = local_from_discrete_lambda(lambda_emp_discrete_safe, cfg);
[lambda_emp_complex_discrete, lambda_emp_complex_valid] = ...
    local_fit_complex_lambda(efuns_raw, cfg);
lambda_emp_complex_input = local_from_discrete_lambda(lambda_emp_complex_discrete, cfg);
lambda_emp_complex_discrete_safe = lambda_emp_complex_discrete;
lambda_emp_complex_discrete_safe(~lambda_emp_complex_valid) = 0;
lambda_emp_complex_input_safe = local_from_discrete_lambda(lambda_emp_complex_discrete_safe, cfg);
[lambda_emp_complex_tau, lambda_emp_complex_kappa] = ...
    local_lambda_timescale(lambda_emp_complex_discrete_safe, cfg.dt);

fit = struct();
fit.cfg = cfg;
fit.feat = cfg.feat;
fit.acf_lags = acf_lags;
fit.acf_mat = acf_mat;
fit.tau_emp = tau_emp;
fit.kappa_emp = kappa_emp;
fit.lambda_ref_input = lambda_ref(:);
fit.lambda_ref_discrete = lambda_ref_discrete(:);
fit.lambda_emp_mag = lambda_emp_mag(:);
fit.lambda_emp_discrete = lambda_emp_discrete(:);
fit.lambda_emp_input = lambda_emp_input(:);
fit.lambda_emp_discrete_safe = lambda_emp_discrete_safe(:);
fit.lambda_emp_input_safe = lambda_emp_input_safe(:);
fit.lambda_emp_complex_discrete = lambda_emp_complex_discrete(:);
fit.lambda_emp_complex_input = lambda_emp_complex_input(:);
fit.lambda_emp_complex_discrete_safe = lambda_emp_complex_discrete_safe(:);
fit.lambda_emp_complex_input_safe = lambda_emp_complex_input_safe(:);
fit.lambda_emp_complex_tau = lambda_emp_complex_tau(:);
fit.lambda_emp_complex_kappa = lambda_emp_complex_kappa(:);
fit.lambda_emp_complex_valid = lambda_emp_complex_valid(:);
fit.used_reference_fallback = invalid(:);
fit.fit_info = fit_info;
end

function cfg = local_set_default(cfg, name, value)
if ~isfield(cfg, name) || isempty(cfg.(name))
    cfg.(name) = value;
end
end

function Phi = local_normalize_efuns(efuns_raw, feat)
if exist('normalize_efun', 'file') == 2
    Phi = normalize_efun(efuns_raw, feat);
    return;
end

if strcmpi(feat, 'real')
    X = real(efuns_raw);
else
    X = abs(efuns_raw);
end

Phi = X;
for j = 1:size(X, 2)
    x = X(:, j);
    mx = max(abs(x));
    if mx > 0
        x = x / mx;
    end
    Phi(:, j) = x;
end
end

function lambda_d = local_to_discrete_lambda(lambda, cfg)
if strcmpi(cfg.lambdaType, 'continuous')
    lambda_d = exp(lambda * cfg.dt);
else
    lambda_d = lambda;
end
end

function lambda_out = local_from_discrete_lambda(lambda_d, cfg)
if strcmpi(cfg.lambdaType, 'continuous')
    lambda_out = log(lambda_d) / cfg.dt;
else
    lambda_out = lambda_d;
end
end

function [lambda_hat, is_valid] = local_fit_complex_lambda(efuns_raw, cfg)
% Fit complex lambda from raw eigenfunctions via one-step regression.
    [~, K] = size(efuns_raw);
    lambda_hat = nan(K, 1);
    is_valid = false(K, 1);

    switch lower(cfg.complex_fit_mode)
        case 'ls_one_step'
            for j = 1:K
                phi = efuns_raw(:, j);
                if numel(phi) < 2
                    continue;
                end
                x_prev = phi(1:end-1);
                x_next = phi(2:end);
                denom = sum(abs(x_prev).^2);
                if ~isfinite(denom) || denom <= 0
                    continue;
                end
                lam = sum(conj(x_prev) .* x_next) / denom;
                if isfinite(real(lam)) && isfinite(imag(lam))
                    lambda_hat(j) = lam;
                    is_valid(j) = true;
                end
            end
        otherwise
            error('Unknown cfg.complex_fit_mode = %s.', cfg.complex_fit_mode);
    end
end

function [tau, kappa] = local_lambda_timescale(lambda_d, dt)
% Convert discrete lambda magnitudes to decay timescale and rate.
    alpha = log(lambda_d) / dt;
    kappa = -real(alpha);
    tau = 1 ./ kappa;

    tau(~isfinite(tau)) = inf;
    tau(kappa <= 0) = inf;
    kappa(~isfinite(kappa) | kappa < 0) = 0;
end

function [lags, acf_mat, tau_emp, fit_info] = local_acf_timescale_from_envelope(Phi, maxLag, dt, cfg)
T = size(Phi, 1);
K = size(Phi, 2);

maxLag = min(maxLag, T - 1);
lags = (0:maxLag).';

acf_mat = nan(maxLag + 1, K);
tau_emp = nan(K, 1);
fit_info = repmat(struct('method','', 'idxStart',nan, 'idxEnd',nan, ...
    'slope',nan, 'r2',nan, 'censored',false), K, 1);

for i = 1:K
    x = Phi(:, i);
    [~, acf_i] = local_acf(x, maxLag);
    acf_i = acf_i(:);
    acf_mat(:, i) = acf_i;

    env = abs(acf_i(:));
    valid = isfinite(env) & (env > 0);

    if nnz(valid) < cfg.fit_min_points
        [tau_emp(i), fit_info(i)] = local_fallback_efold(env, valid, lags, dt);
        continue;
    end

    valid(1) = false;
    vidx = find(valid);
    if isempty(vidx)
        [tau_emp(i), fit_info(i)] = local_fallback_efold(env, isfinite(env) & env > 0, lags, dt);
        continue;
    end

    idxStartCand = find(valid & (env <= cfg.fit_env_start), 1, 'first');
    if isempty(idxStartCand)
        idxStartCand = vidx(1);
    end

    idxEndCand = find(valid & ((1:numel(env))' >= idxStartCand) & (env <= cfg.fit_env_end), 1, 'first');
    censored = false;
    if isempty(idxEndCand)
        idxEndCand = vidx(end);
        censored = true;
    end

    fitIdx = vidx(vidx >= idxStartCand & vidx <= idxEndCand);
    fitIdx = fitIdx(env(fitIdx) >= cfg.fit_env_floor);

    if numel(fitIdx) < cfg.fit_min_points
        cand = vidx(vidx <= idxEndCand & env(vidx) >= cfg.fit_env_floor);
        if numel(cand) >= cfg.fit_min_points
            fitIdx = cand(end - cfg.fit_min_points + 1:end);
        else
            if strcmpi(cfg.fit_fallback, 'nan')
                tau_emp(i) = nan;
                fit_info(i).method = 'nan';
                fit_info(i).censored = censored;
            else
                [tau_emp(i), fit_info(i)] = local_fallback_efold(env, valid, lags, dt);
                fit_info(i).censored = censored;
            end
            continue;
        end
    end

    fitIdx = unique(fitIdx(:));
    fitIdx = fitIdx(fitIdx >= 1 & fitIdx <= numel(env));
    if numel(fitIdx) < cfg.fit_min_points
        [tau_emp(i), fit_info(i)] = local_fallback_efold(env, valid, lags, dt);
        fit_info(i).censored = censored;
        continue;
    end

    xfit = lags(fitIdx) * dt;
    yfit = log(env(fitIdx));
    p = polyfit(xfit, yfit, 1);
    slope = p(1);

    if ~isfinite(slope) || slope >= 0
        [tau_emp(i), fit_info(i)] = local_fallback_efold(env, valid, lags, dt);
        fit_info(i).censored = censored;
        continue;
    end

    yhat = polyval(p, xfit);
    ss_res = sum((yfit - yhat).^2);
    ss_tot = sum((yfit - mean(yfit)).^2);
    if ss_tot > 0
        r2 = 1 - ss_res / ss_tot;
    else
        r2 = nan;
    end

    tau_emp(i) = -1 / slope;
    fit_info(i).method = 'loglin';
    fit_info(i).idxStart = fitIdx(1);
    fit_info(i).idxEnd = fitIdx(end);
    fit_info(i).slope = slope;
    fit_info(i).r2 = r2;
    fit_info(i).censored = censored;
end
end

function [lagsPositive, acfPositive] = local_acf(x, maxLag)
x = x(:) - mean(x);
[c, lags] = xcorr(x, maxLag, 'coeff');
idx = lags >= 0;
lagsPositive = lags(idx);
acfPositive = c(idx);
end

function [tau, info] = local_fallback_efold(env, valid, lags, dt)
info = struct('method','efold', 'idxStart',nan, 'idxEnd',nan, ...
    'slope',nan, 'r2',nan, 'censored',false);

env = env(:);
lags = lags(:);

if nargin < 2 || isempty(valid)
    valid = isfinite(env) & (env > 0);
else
    valid = valid(:);
end

valid(1) = false;

thr = exp(-1);
idx = find(valid & (env <= thr), 1, 'first');
if isempty(idx)
    last = find(valid, 1, 'last');
    if isempty(last)
        tau = nan;
        info.method = 'nan';
    else
        tau = lags(last) * dt;
        info.censored = true;
    end
    return;
end

tau = lags(idx) * dt;
info.idxStart = idx;
info.idxEnd = idx;
end
