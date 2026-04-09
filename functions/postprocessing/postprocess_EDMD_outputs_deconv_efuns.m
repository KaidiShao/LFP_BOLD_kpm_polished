function [EDMD_outputs, fig, deconv] = postprocess_EDMD_outputs_deconv_efuns(EDMD_outputs, cfg)
%POSTPROCESS_EDMD_OUTPUTS_DECONV_EFUNS  Estimate internal perturbations from Koopman eigenfunctions.
%
% This function:
%   1) Builds the "ALL" (thresholded & sorted) and "SELECTED" (first K) eigenfunction sets
%      following the same logic used by postprocess_EDMD_outputs.m.
%   2) Estimates an internal perturbation signal u(:,j) for each eigenfunction phi(:,j)
%      using either:
%        - Wiener deconvolution with h_j[k] = lambda_j^(k-1), or
%        - Koopman residual recursion u_t = phi_t - lambda * phi_{t-1}.
%   3) Produces deconvolved time series u(:,j) and two visualization modes:
%        - abs mode  : normalize(abs(u)) column-wise
%        - real mode : normalize(real(u)) column-wise
%      (same style as postprocess_EDMD_outputs.m)
%   4) Optionally saves the deconvolved eigenfunctions to disk.
%   5) Optionally plots a diagnostic figure (2x7 layout) with 4 heatmaps:
%        (ALL abs, SELECTED abs, ALL real, SELECTED real)
%
% Inputs
%   EDMD_outputs : struct (ideally already processed by postprocess_EDMD_outputs)
%   cfg         : struct with options (all optional)
%
% Key cfg fields (defaults)
%   cfg.method          = 'wiener';    % 'wiener' or 'koopman_residual'
%   cfg.lambda_source   = 'edmd';      % 'edmd', 'empirical_abs', 'empirical_real', or 'empirical_complex'
%   cfg.lambdaType      = 'discrete';  % 'discrete' or 'continuous'
%   cfg.dt              = 1;           % used only if lambdaType='continuous'
%   cfg.remove_mean     = false;       % optional preprocessing before deconvolution
%   cfg.first_u_mode    = 'phi1';      % for koopman_residual: 'phi1' or 'zero'
%   cfg.empirical_lambda_phase_mode = 'preserve_edmd_angle'
%   cfg.empirical_maxLag = [];         % default min(200, T-1)
%   cfg.fit_env_start    = 0.9;
%   cfg.fit_env_end      = 0.2;
%   cfg.fit_env_floor    = 1e-3;
%   cfg.fit_min_points   = 8;
%   cfg.fit_fallback     = 'efold';
%
%   cfg.impulseLen      = [];          % fixed impulse length; [] = auto based on tailTol
%   cfg.tailTol         = 1e-3;        % for auto impulseLen: abs(lambda)^(Lh-1) <= tailTol
%   cfg.maxImpulseLen   = 5000;        % cap for auto impulseLen
%
%   cfg.K               = 1e-2;        % Wiener constant (noise-to-signal ratio)
%   cfg.snrDb           = [];          % alternative: K = 10^(-snrDb/10)
%
%   cfg.max_modes_all   = Inf;         % cap number of modes in ALL set (for speed)
%   cfg.max_modes_sel   = Inf;         % cap number of modes in SELECTED set
%
%   cfg.do_plot         = true;
%   cfg.t_plot          = [];          % explicit indices for the plot window
%   cfg.window_idx      = [];          % alias of cfg.t_plot for windowed plotting
%   cfg.max_plot_samples = 2000;       % default window length when no window is supplied
%   cfg.window_start    = 1;           % default window start when no window is supplied
%   cfg.time_vec        = [];          % time axis vector; [] => use indices
%   cfg.colormap_fun    = [];          % e.g. @() flipud(othercolor('Spectral10'))
%   cfg.session_border  = [];          % optional xline positions in current x-axis units
%   cfg.draw_border     = false;
%   cfg.plot_normalize_scope = 'global'; % 'window' or 'global' for heatmap normalisation
%   cfg.normalize_exclude_idx = 1;       % indices excluded when computing global column maxima
%
%   cfg.save_path       = '';          % if non-empty: save to this .mat
%   cfg.save_dir        = '';          % alternative: directory to save
%   cfg.save_prefix     = 'deconv';    % used with save_dir
%   cfg.save_v7_3       = true;        % use -v7.3
%
% Outputs
%   EDMD_outputs : input struct augmented with EDMD_outputs.deconv_efuns
%   fig         : figure handle (empty if cfg.do_plot=false)
%   deconv      : the saved struct EDMD_outputs.deconv_efuns
%
% Note
%   The Wiener option uses frequency-domain deconvolution for speed.
%   The koopman_residual option follows the discrete recursion
%   phi_t = lambda_d * phi_{t-1} + u_t.

% -------------------- defaults --------------------
if nargin < 2, cfg = struct(); end

cfg = local_set_default(cfg, 'method',         'wiener');
cfg = local_set_default(cfg, 'lambda_source',  'edmd');
cfg = local_set_default(cfg, 'lambdaType',     'discrete');
cfg = local_set_default(cfg, 'dt',             1);
cfg = local_set_default(cfg, 'remove_mean',    false);
cfg = local_set_default(cfg, 'first_u_mode',   'phi1');
cfg = local_set_default(cfg, 'empirical_lambda_phase_mode', 'preserve_edmd_angle');
cfg = local_set_default(cfg, 'empirical_maxLag', []);
cfg = local_set_default(cfg, 'fit_env_start', 0.9);
cfg = local_set_default(cfg, 'fit_env_end', 0.2);
cfg = local_set_default(cfg, 'fit_env_floor', 1e-3);
cfg = local_set_default(cfg, 'fit_min_points', 8);
cfg = local_set_default(cfg, 'fit_fallback', 'efold');

cfg = local_set_default(cfg, 'impulseLen',     []);
cfg = local_set_default(cfg, 'tailTol',        1e-3);
cfg = local_set_default(cfg, 'maxImpulseLen',  5000);

cfg = local_set_default(cfg, 'K',              1e-2);
cfg = local_set_default(cfg, 'snrDb',          []);

cfg = local_set_default(cfg, 'max_modes_all',  Inf);
cfg = local_set_default(cfg, 'max_modes_sel',  Inf);

cfg = local_set_default(cfg, 'do_plot',        true);
cfg = local_set_default(cfg, 't_plot',         []);
cfg = local_set_default(cfg, 'window_idx',     []);
cfg = local_set_default(cfg, 'max_plot_samples', 2000);
cfg = local_set_default(cfg, 'window_start',   1);
cfg = local_set_default(cfg, 'time_vec',       []);
cfg = local_set_default(cfg, 'colormap_fun',   []);
cfg = local_set_default(cfg, 'session_border', []);
cfg = local_set_default(cfg, 'draw_border',    false);
cfg = local_set_default(cfg, 'plot_normalize_scope', 'global');
cfg = local_set_default(cfg, 'normalize_exclude_idx', 1);

cfg = local_set_default(cfg, 'save_path',      '');
cfg = local_set_default(cfg, 'save_dir',       '');
cfg = local_set_default(cfg, 'save_prefix',    'deconv');
cfg = local_set_default(cfg, 'save_v7_3',      true);

% -------------------- build ALL and SELECTED sets --------------------
[efuns_all, evalues_all, efuns_sel, evalues_sel, meta] = local_get_sets(EDMD_outputs);

T = size(efuns_all, 1);
K_all = size(efuns_all, 2);
K_sel = size(efuns_sel, 2);

K_all_use = min(K_all, cfg.max_modes_all);
K_sel_use = min(K_sel, cfg.max_modes_sel);
K_deconv  = min(K_all, max(K_all_use, K_sel_use));

% Check whether SELECTED set matches the prefix of ALL set.
% If not, SELECTED will be deconvolved separately so that plots/saved data are correct.
Kcmp = min(K_sel_use, K_all);
tol_prefix = 1e-12;
is_sel_prefix = (Kcmp == 0) || all(abs(evalues_sel(1:Kcmp) - evalues_all(1:Kcmp)) < tol_prefix);

evalues_all_use = evalues_all(1:K_deconv);
efuns_all_use = efuns_all(:, 1:K_deconv);
if is_sel_prefix
    evalues_sel_use = evalues_all_use(1:K_sel_use);
    efuns_sel_use = efuns_all_use(:, 1:K_sel_use);
else
    evalues_sel_use = evalues_sel(1:K_sel_use);
    efuns_sel_use = efuns_sel(:, 1:K_sel_use);
end

lambda_bundle_all = local_resolve_lambda_bundle( ...
    EDMD_outputs, efuns_all_use, evalues_all_use, 'all', cfg);
lambda_bundle_sel = local_resolve_lambda_bundle( ...
    EDMD_outputs, efuns_sel_use, evalues_sel_use, 'selected', cfg);

% -------------------- internal perturbation estimate --------------------

U = zeros(T, K_deconv, 'like', efuns_all);
impulseLen = zeros(K_deconv, 1);
lambda_d_used = zeros(K_deconv, 1);
K_used = zeros(K_deconv, 1);

for j = 1:K_deconv
    phi = efuns_all(:, j);
    lam = lambda_bundle_all.used_input(j);
    if cfg.remove_mean
        phi = phi - mean(phi);
    end

    [u, info] = local_compute_internal_perturbation(phi, lam, cfg);
    U(:, j) = u;
    impulseLen(j)   = info.impulseLen;
    lambda_d_used(j)= info.lambda_d;
    K_used(j)       = info.K;
end

% If SELECTED is not a prefix of ALL, deconvolve SELECTED separately.
U_sel = [];
impulseLen_sel = [];
lambda_d_sel = [];
K_sel_wiener = [];
if ~is_sel_prefix
    U_sel = zeros(T, K_sel_use, 'like', efuns_sel);
    impulseLen_sel = zeros(K_sel_use, 1);
    lambda_d_sel = zeros(K_sel_use, 1);
    K_sel_wiener = zeros(K_sel_use, 1);
    for j = 1:K_sel_use
        phi = efuns_sel(:, j);
        lam = lambda_bundle_sel.used_input(j);
        if cfg.remove_mean
            phi = phi - mean(phi);
        end
        [u, info] = local_compute_internal_perturbation(phi, lam, cfg);
        U_sel(:, j) = u;
        impulseLen_sel(j) = info.impulseLen;
        lambda_d_sel(j) = info.lambda_d;
        K_sel_wiener(j) = info.K;
    end
end

% -------------------- normalize for visualization --------------------
U_abs  = local_normalize_cols_excluding(abs(U), cfg.normalize_exclude_idx);
U_real = local_normalize_cols_excluding(real(U), cfg.normalize_exclude_idx);

% ALL and SELECTED views
deconv = struct();
deconv.cfg = cfg;
deconv.meta = meta;
deconv.method = lower(cfg.method);
deconv.method_label = local_method_label(cfg.method);
deconv.formula = local_method_formula(cfg);
deconv.lambda_source = lower(cfg.lambda_source);
deconv.lambda_source_label = local_lambda_source_label(cfg.lambda_source);
deconv.lambda_source_formula = local_lambda_source_formula(cfg);
deconv.evalues_ref_all = evalues_all_use;
deconv.evalues_all = lambda_bundle_all.used_input;
deconv.evalues_used_all = lambda_bundle_all.used_input;
deconv.lambda_d    = lambda_d_used;
deconv.impulseLen  = impulseLen;
deconv.wienerK     = K_used;
deconv.empirical_fit_all = lambda_bundle_all.empirical_fit;

deconv.u_all  = U(:, 1:K_all_use);
if is_sel_prefix
    deconv.u_sel = U(:, 1:K_sel_use);
    deconv.evalues_ref_sel = evalues_sel_use;
    deconv.evalues_sel = lambda_bundle_sel.used_input;
    deconv.evalues_used_sel = lambda_bundle_sel.used_input;
    deconv.lambda_d_sel = lambda_bundle_all.used_discrete(1:K_sel_use);
    deconv.empirical_fit_sel = local_trim_empirical_fit(lambda_bundle_all.empirical_fit, K_sel_use);
else
    deconv.u_sel = U_sel;
    deconv.evalues_ref_sel = evalues_sel_use;
    deconv.evalues_sel = lambda_bundle_sel.used_input;
    deconv.evalues_used_sel = lambda_bundle_sel.used_input;
    deconv.lambda_d_sel = lambda_d_sel;
    deconv.impulseLen_sel = impulseLen_sel;
    deconv.wienerK_sel = K_sel_wiener;
    deconv.empirical_fit_sel = lambda_bundle_sel.empirical_fit;
end

deconv.norm_efuns.abs_all  = U_abs(:,  1:K_all_use);
deconv.norm_efuns.real_all = U_real(:, 1:K_all_use);
if is_sel_prefix
    deconv.norm_efuns.abs_sel  = U_abs(:,  1:K_sel_use);
    deconv.norm_efuns.real_sel = U_real(:, 1:K_sel_use);
else
    U_sel_abs  = local_normalize_cols_excluding(abs(U_sel), cfg.normalize_exclude_idx);
    U_sel_real = local_normalize_cols_excluding(real(U_sel), cfg.normalize_exclude_idx);
    deconv.norm_efuns.abs_sel  = U_sel_abs;
    deconv.norm_efuns.real_sel = U_sel_real;
end

% Cache the exact matrices used for the heatmaps so they can be inspected later.
t_plot = local_resolve_plot_window(T, cfg);
[A_abs_all, A_abs_sel, A_re_all, A_re_sel, norm_label] = ...
    local_prepare_plot_mats(deconv, t_plot, cfg);

deconv.plot_view = struct();
deconv.plot_view.t_plot = t_plot;
deconv.plot_view.norm_label = norm_label;
deconv.plot_view.normalize_exclude_idx = cfg.normalize_exclude_idx;
deconv.plot_view.abs_all = A_abs_all;
deconv.plot_view.abs_sel = A_abs_sel;
deconv.plot_view.real_all = A_re_all;
deconv.plot_view.real_sel = A_re_sel;

if isempty(cfg.time_vec)
    deconv.plot_view.x = t_plot;
else
    deconv.plot_view.x = cfg.time_vec(t_plot);
end

% Attach to outputs
EDMD_outputs.deconv_efuns = deconv;

% -------------------- save to disk (optional) --------------------
save_path = local_resolve_save_path(cfg);
if ~isempty(save_path)
    try
        if cfg.save_v7_3
            save(save_path, 'deconv', '-v7.3');
        else
            save(save_path, 'deconv');
        end
        EDMD_outputs.deconv_efuns.save_path = save_path;
    catch ME
        warning('Failed to save deconvolved eigenfunctions to "%s": %s', save_path, ME.message);
    end
end

% -------------------- plot (optional) --------------------
fig = [];
if cfg.do_plot
    fig = local_plot_deconv(EDMD_outputs, deconv, efuns_all, evalues_all, evalues_sel, cfg);
end

end

% =====================================================================
% Helpers
% =====================================================================

function cfg = local_set_default(cfg, name, value)
    if ~isfield(cfg, name)
        cfg.(name) = value;
    end
end

function [efuns_all, evalues_all, efuns_sel, evalues_sel, meta] = local_get_sets(EDMD_outputs)
% Build sets consistent with postprocess_EDMD_outputs.m
    meta = struct();

    if isfield(EDMD_outputs, 'original_sorted') && ...
            isfield(EDMD_outputs.original_sorted, 'evalues') && ...
            isfield(EDMD_outputs.original_sorted, 'efuns')

        e_full = EDMD_outputs.original_sorted.evalues(:);
        f_full = EDMD_outputs.original_sorted.efuns;
        if size(f_full, 2) ~= numel(e_full)
            error('Dimension mismatch: original_sorted.efuns columns != numel(original_sorted.evalues).');
        end

        if isfield(EDMD_outputs.original_sorted, 'abs_thresh')
            abs_thresh = EDMD_outputs.original_sorted.abs_thresh;
        elseif isfield(EDMD_outputs, 'masked') && isfield(EDMD_outputs.masked, 'abs_thresh')
            abs_thresh = EDMD_outputs.masked.abs_thresh;
        else
            abs_thresh = 0;
        end

        mask_sorted = abs(e_full) > abs_thresh;
        if ~any(mask_sorted)
            error('All eigenvalues removed by abs_thresh=%g in original_sorted.', abs_thresh);
        end

        efuns_all   = f_full(:, mask_sorted);
        evalues_all = e_full(mask_sorted);
        meta.abs_thresh = abs_thresh;
        meta.source = 'original_sorted + abs_thresh';

    else
        if ~isfield(EDMD_outputs, 'evalues') || ~isfield(EDMD_outputs, 'efuns')
            error('EDMD_outputs must contain fields evalues and efuns (or original_sorted.evalues/efuns).');
        end
        evalues_all = EDMD_outputs.evalues(:);
        efuns_all   = EDMD_outputs.efuns;
        if size(efuns_all, 2) ~= numel(evalues_all)
            error('Dimension mismatch: efuns columns != numel(evalues).');
        end
        meta.source = 'EDMD_outputs.evalues/efuns';
    end

    % SELECTED set: prefer EDMD_outputs.evalues/efuns if present (already truncated)
    if isfield(EDMD_outputs, 'evalues') && isfield(EDMD_outputs, 'efuns')
        evalues_sel = EDMD_outputs.evalues(:);
        efuns_sel   = EDMD_outputs.efuns;
        if size(efuns_sel, 2) ~= numel(evalues_sel)
            error('Dimension mismatch: efuns columns != numel(evalues) for selected set.');
        end
    else
        evalues_sel = evalues_all;
        efuns_sel   = efuns_all;
    end

    if size(efuns_sel, 1) ~= size(efuns_all, 1)
        error('Time dimension mismatch between ALL and SELECTED eigenfunctions.');
    end
end

function bundle = local_resolve_lambda_bundle(EDMD_outputs, efuns_use, evalues_use, set_name, cfg)
% Resolve which eigenvalues should drive the internal perturbation estimate.
    bundle = struct();
    bundle.source = lower(cfg.lambda_source);
    bundle.source_label = local_lambda_source_label(cfg.lambda_source);
    bundle.reference_input = evalues_use(:);
    bundle.reference_discrete = local_to_discrete_lambda(evalues_use(:), cfg);
    bundle.empirical_fit = [];

    switch lower(cfg.lambda_source)
        case 'edmd'
            bundle.used_input = bundle.reference_input;
            bundle.used_discrete = bundle.reference_discrete;

        case {'empirical_abs', 'empirical_real'}
            if contains(lower(cfg.lambda_source), 'real')
                feat = 'real';
            else
                feat = 'abs';
            end

            fit = local_fetch_empirical_lambda_fit(EDMD_outputs, efuns_use, evalues_use, set_name, feat, cfg);
            bundle.empirical_fit = fit;
            bundle.used_input = fit.lambda_emp_input_safe(:);
            bundle.used_discrete = fit.lambda_emp_discrete_safe(:);

        case 'empirical_complex'
            fit = local_fetch_empirical_lambda_fit(EDMD_outputs, efuns_use, evalues_use, set_name, 'abs', cfg);
            bundle.empirical_fit = fit;
            bundle.used_input = fit.lambda_emp_complex_input_safe(:);
            bundle.used_discrete = fit.lambda_emp_complex_discrete_safe(:);

        otherwise
            error('Unknown cfg.lambda_source = %s.', cfg.lambda_source);
    end
end

function fit = local_fetch_empirical_lambda_fit(EDMD_outputs, efuns_use, evalues_use, set_name, feat, cfg)
% Reuse timescale results when available; otherwise fit on the fly.
    if isfield(EDMD_outputs, 'timescale_info') && ~isempty(EDMD_outputs.timescale_info)
        try
            ts_info = EDMD_outputs.timescale_info.(set_name).(feat);
            if isfield(ts_info, 'empirical') && numel(ts_info.empirical.lambda_emp_input_safe) >= numel(evalues_use)
                fit = local_trim_empirical_fit(ts_info.empirical, numel(evalues_use));
                return;
            end
        catch
        end
    end

    fit_cfg = struct();
    fit_cfg.dt = cfg.dt;
    fit_cfg.lambdaType = cfg.lambdaType;
    fit_cfg.feat = feat;
    fit_cfg.maxLag = cfg.empirical_maxLag;
    fit_cfg.fit_env_start = cfg.fit_env_start;
    fit_cfg.fit_env_end = cfg.fit_env_end;
    fit_cfg.fit_env_floor = cfg.fit_env_floor;
    fit_cfg.fit_min_points = cfg.fit_min_points;
    fit_cfg.fit_fallback = cfg.fit_fallback;
    fit_cfg.lambda_phase_mode = cfg.empirical_lambda_phase_mode;
    fit = fit_empirical_eigenvalues_from_efuns(efuns_use, evalues_use, fit_cfg);
end

function fit = local_trim_empirical_fit(fit, K)
% Trim empirical-fit fields to the first K modes.
    if isempty(fit)
        return;
    end

    fields = {'tau_emp', 'kappa_emp', 'lambda_ref_input', 'lambda_ref_discrete', ...
        'lambda_emp_mag', 'lambda_emp_discrete', 'lambda_emp_input', ...
        'lambda_emp_discrete_safe', 'lambda_emp_input_safe', ...
        'lambda_emp_complex_discrete', 'lambda_emp_complex_input', ...
        'lambda_emp_complex_discrete_safe', 'lambda_emp_complex_input_safe', ...
        'lambda_emp_complex_tau', 'lambda_emp_complex_kappa', ...
        'lambda_emp_complex_valid', 'used_reference_fallback'};
    for i = 1:numel(fields)
        name = fields{i};
        if isfield(fit, name) && ~isempty(fit.(name))
            fit.(name) = fit.(name)(1:K);
        end
    end

    if isfield(fit, 'acf_mat') && ~isempty(fit.acf_mat)
        fit.acf_mat = fit.acf_mat(:, 1:K);
    end
    if isfield(fit, 'fit_info') && ~isempty(fit.fit_info)
        fit.fit_info = fit.fit_info(1:K);
    end
end

function [u, info] = local_compute_internal_perturbation(phi, lambda, cfg)
% Dispatch the chosen internal-perturbation estimator.
    switch lower(cfg.method)
        case 'wiener'
            [u, info] = local_wiener_deconv(phi, lambda, cfg);
        case {'koopman_residual', 'residual'}
            [u, info] = local_koopman_residual(phi, lambda, cfg);
        otherwise
            error(['cfg.method="%s" not supported. Use ''wiener'' or ', ...
                '''koopman_residual''.'], cfg.method);
    end
end

function [u, info] = local_wiener_deconv(phi, lambda, cfg)
% Frequency-domain Wiener deconvolution with exponential impulse response.
    phi = phi(:);
    T = numel(phi);

    % Eigenvalue -> discrete
    lambda_d = local_to_discrete_lambda(lambda, cfg);

    % Choose impulse length
    if ~isempty(cfg.impulseLen)
        Lh = min(T, round(cfg.impulseLen));
    else
        if abs(lambda_d) < 1
            Lh_auto = ceil(log(cfg.tailTol) / log(abs(lambda_d))) + 1;
            if ~isfinite(Lh_auto) || Lh_auto < 1
                Lh_auto = min(T, cfg.maxImpulseLen);
            end
            Lh = min(T, min(cfg.maxImpulseLen, Lh_auto));
        else
            Lh = min(T, cfg.maxImpulseLen);
        end
    end
    Lh = max(1, Lh);

    h = (lambda_d .^ (0:Lh-1)).';

    if ~isempty(cfg.snrDb)
        K = 10^(-cfg.snrDb/10);
    else
        K = cfg.K;
    end

    Nfft = 2^nextpow2(T + Lh - 1);
    Hf   = fft(h,   Nfft);
    Phif = fft(phi, Nfft);

    Uf = conj(Hf) .* Phif ./ (abs(Hf).^2 + K);
    u_full = ifft(Uf, Nfft);
    u = u_full(1:T);

    info = struct();
    info.lambda_input = lambda;
    info.lambda_d = lambda_d;
    info.impulseLen = Lh;
    info.K = K;
end

function [u, info] = local_koopman_residual(phi, lambda, cfg)
% Discrete Koopman recursion residual: u_t = phi_t - lambda_d * phi_{t-1}.
    phi = phi(:);
    T = numel(phi);
    lambda_d = local_to_discrete_lambda(lambda, cfg);

    u = zeros(T, 1, 'like', phi);
    if T < 1
        info = struct('lambda_input', lambda, 'lambda_d', lambda_d, ...
            'impulseLen', NaN, 'K', NaN);
        return;
    end

    switch lower(cfg.first_u_mode)
        case 'phi1'
            u(1) = phi(1);
        case 'zero'
            u(1) = zeros(1, 1, 'like', phi);
        otherwise
            error('Unknown cfg.first_u_mode = %s. Use ''phi1'' or ''zero''.', cfg.first_u_mode);
    end

    if T >= 2
        u(2:end) = phi(2:end) - lambda_d .* phi(1:end-1);
    end

    info = struct();
    info.lambda_input = lambda;
    info.lambda_d = lambda_d;
    info.impulseLen = NaN;
    info.K = NaN;
end

function lambda_d = local_to_discrete_lambda(lambda, cfg)
% Convert eigenvalue input to discrete-time form when needed.
    if strcmpi(cfg.lambdaType, 'continuous')
        lambda_d = exp(lambda * cfg.dt);
    else
        lambda_d = lambda;
    end
end

function Xn = local_normalize_cols(X)
% Column-wise: divide by max(abs) only, without de-meaning.
    Xn = X;
    for j = 1:size(X, 2)
        x = X(:, j);
        mx = max(abs(x));
        if mx > 0
            x = x ./ mx;
        end
        Xn(:, j) = x;
    end
end

function Xn = local_normalize_cols_excluding(X, exclude_idx)
% Column-wise normalisation using maxima computed after excluding selected rows.
    Xn = X;
    n_rows = size(X, 1);

    if nargin < 2 || isempty(exclude_idx)
        exclude_mask = false(n_rows, 1);
    else
        exclude_idx = unique(round(exclude_idx(:)));
        exclude_idx = exclude_idx(exclude_idx >= 1 & exclude_idx <= n_rows);
        exclude_mask = false(n_rows, 1);
        exclude_mask(exclude_idx) = true;
    end

    include_mask = ~exclude_mask;
    if ~any(include_mask)
        include_mask(:) = true;
    end

    for j = 1:size(X, 2)
        x = X(:, j);
        mx = max(abs(x(include_mask)));
        if mx > 0
            x = x ./ mx;
        end
        Xn(:, j) = x;
    end
end

function save_path = local_resolve_save_path(cfg)
% Decide where to save, if requested.
    save_path = '';
    if isfield(cfg, 'save_path') && ~isempty(cfg.save_path)
        save_path = cfg.save_path;
        return;
    end
    if isfield(cfg, 'save_dir') && ~isempty(cfg.save_dir)
        if ~exist(cfg.save_dir, 'dir')
            mkdir(cfg.save_dir);
        end
        prefix = cfg.save_prefix;
        if isempty(prefix)
            prefix = 'deconv';
        end
        save_path = fullfile(cfg.save_dir, [prefix '_deconv_efuns.mat']);
    end
end

function fig = local_plot_deconv(~, deconv, efuns_all, ~, ~, cfg)
% Plot deconvolved eigenfunctions in the same layout as postprocess_EDMD_outputs.m

    T = size(efuns_all, 1);
    t_plot = local_resolve_plot_window(T, cfg);
    signal_label = deconv.method_label;
    lambda_label = deconv.lambda_source_label;

    if isempty(cfg.time_vec)
        x = t_plot;
        xlab = 'Time (index)';
    else
        x = cfg.time_vec(t_plot);
        xlab = 'Time';
    end

    % Session borders
    session_border = cfg.session_border;

    % Select how many modes are available in saved deconv
    K_all_use = size(deconv.norm_efuns.abs_all, 2);
    K_sel_use = size(deconv.norm_efuns.abs_sel, 2);

    [A_abs_all, A_abs_sel, A_re_all, A_re_sel, norm_label] = ...
        local_prepare_plot_mats(deconv, t_plot, cfg);

    fig = figure('Color', 'w');

    % (1,1) original eigenvalues (sorted)
    subplot(2,7,1);
    if strcmpi(deconv.lambda_source, 'edmd')
        plot(real(deconv.evalues_used_all(:)), imag(deconv.evalues_used_all(:)), '.', 'MarkerSize', 10);
    else
        plot(real(deconv.evalues_ref_all(:)), imag(deconv.evalues_ref_all(:)), '.', ...
            'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
        hold on;
        plot(real(deconv.evalues_used_all(:)), imag(deconv.evalues_used_all(:)), '.', ...
            'Color', [0 0.4470 0.7410], 'MarkerSize', 10);
    end
    hold on;
    th = linspace(0, 2*pi, 400);
    plot(cos(th), sin(th), 'k--', 'LineWidth', 1);
    hold off;
    grid on; axis equal; axis([-1.1, 1.1, -1.1, 1.1]);
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title(sprintf('Eigenvalues used (%s)', lambda_label), 'Interpreter', 'none');

    % (2,1) selected eigenvalues
    subplot(2,7,8);
    if strcmpi(deconv.lambda_source, 'edmd')
        plot(real(deconv.evalues_used_sel(1:K_sel_use)), imag(deconv.evalues_used_sel(1:K_sel_use)), '.', 'MarkerSize', 10);
    else
        plot(real(deconv.evalues_ref_sel(1:K_sel_use)), imag(deconv.evalues_ref_sel(1:K_sel_use)), '.', ...
            'Color', [0.7 0.7 0.7], 'MarkerSize', 8);
        hold on;
        plot(real(deconv.evalues_used_sel(1:K_sel_use)), imag(deconv.evalues_used_sel(1:K_sel_use)), '.', ...
            'Color', [0 0.4470 0.7410], 'MarkerSize', 10);
    end
    hold on;
    plot(cos(th), sin(th), 'k--', 'LineWidth', 1);
    hold off;
    grid on; axis equal; axis([-1.1, 1.1, -1.1, 1.1]);
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title(sprintf('Selected eigenvalues (%s, K=%d)', lambda_label, K_sel_use), 'Interpreter', 'none');

    % Choose colormap
    if isempty(cfg.colormap_fun)
        if exist('othercolor', 'file') == 2
            use_cmap = @() flipud(othercolor('Spectral10'));
        else
            use_cmap = @() parula(256);
        end
    else
        use_cmap = cfg.colormap_fun;
    end

    % (1,2) ALL abs
    subplot(2,7,2:4);
    imagesc(x, 1:size(A_abs_all,1), A_abs_all);
    if cfg.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir', 'reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title(sprintf('|%s| (ALL, K=%d, %s, %s)', signal_label, K_all_use, lambda_label, norm_label), 'Interpreter', 'none');
    colorbar;
    colormap(use_cmap());

    % (2,2) SELECTED abs
    subplot(2,7,9:11);
    imagesc(x, 1:size(A_abs_sel,1), A_abs_sel);
    if cfg.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir', 'reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title(sprintf('|%s| (SELECTED, K=%d, %s, %s)', signal_label, K_sel_use, lambda_label, norm_label), 'Interpreter', 'none');
    colorbar;
    colormap(use_cmap());

    % (1,3) ALL real
    subplot(2,7,5:7);
    imagesc(x, 1:size(A_re_all,1), A_re_all);
    if cfg.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir', 'reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title(sprintf('Re(%s) (ALL, K=%d, %s, %s)', signal_label, K_all_use, lambda_label, norm_label), 'Interpreter', 'none');
    colorbar;
    colormap(use_cmap());

    % (2,3) SELECTED real
    subplot(2,7,12:14);
    imagesc(x, 1:size(A_re_sel,1), A_re_sel);
    if cfg.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir', 'reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title(sprintf('Re(%s) (SELECTED, K=%d, %s, %s)', signal_label, K_sel_use, lambda_label, norm_label), 'Interpreter', 'none');
    colorbar;
    colormap(use_cmap());

    % Reasonable default size
    set(fig, 'Position', [1440, 654, 1531, 584]);
end


function [A_abs_all, A_abs_sel, A_re_all, A_re_sel, norm_label] = local_prepare_plot_mats(deconv, t_plot, cfg)
% Build heatmap matrices using either full-series or window-only normalization.

switch lower(cfg.plot_normalize_scope)
    case 'global'
        A_abs_all = deconv.norm_efuns.abs_all(t_plot, :).';
        A_abs_sel = deconv.norm_efuns.abs_sel(t_plot, :).';
        A_re_all  = deconv.norm_efuns.real_all(t_plot, :).';
        A_re_sel  = deconv.norm_efuns.real_sel(t_plot, :).';
        if isempty(cfg.normalize_exclude_idx)
            norm_label = 'global-norm';
        else
            norm_label = 'global-norm (skip t=1)';
        end

    case 'window'
        A_abs_all = local_normalize_cols(abs(deconv.u_all(t_plot, :))).';
        A_abs_sel = local_normalize_cols(abs(deconv.u_sel(t_plot, :))).';
        A_re_all  = local_normalize_cols(real(deconv.u_all(t_plot, :))).';
        A_re_sel  = local_normalize_cols(real(deconv.u_sel(t_plot, :))).';
        norm_label = 'window-norm';

    otherwise
        error('Unknown cfg.plot_normalize_scope = %s. Use ''window'' or ''global''.', ...
            cfg.plot_normalize_scope);
end
end


function t_plot = local_resolve_plot_window(T, cfg)
% Resolve a single visualization window while keeping backward compatibility.

if ~isempty(cfg.t_plot)
    idx = cfg.t_plot;
elseif ~isempty(cfg.window_idx)
    idx = cfg.window_idx;
else
    win_len = min(T, max(1, round(cfg.max_plot_samples)));
    start_idx = min(max(1, round(cfg.window_start)), T - win_len + 1);
    idx = start_idx:(start_idx + win_len - 1);
end

t_plot = unique(round(idx(:).'));
t_plot = t_plot(t_plot >= 1 & t_plot <= T);

if isempty(t_plot)
    t_plot = 1:T;
end
end

function label = local_method_label(method)
% Human-readable label for plot titles.
    switch lower(method)
        case 'wiener'
            label = 'deconv efuns';
        case {'koopman_residual', 'residual'}
            label = 'koopman residual';
        otherwise
            label = lower(method);
    end
end

function label = local_lambda_source_label(lambda_source)
% Human-readable label for eigenvalue source.
    switch lower(lambda_source)
        case 'edmd'
            label = 'EDMD-lambda';
        case 'empirical_abs'
            label = 'empirical-lambda(abs)';
        case 'empirical_real'
            label = 'empirical-lambda(real)';
        case 'empirical_complex'
            label = 'empirical-lambda(complex)';
        otherwise
            label = lower(lambda_source);
    end
end

function formula = local_lambda_source_formula(cfg)
% Text description of how lambda was chosen for the perturbation analysis.
    switch lower(cfg.lambda_source)
        case 'edmd'
            formula = 'lambda used directly from EDMD output';
        case 'empirical_abs'
            formula = 'lambda magnitude from tau_emp(|phi|), phase preserved from EDMD';
        case 'empirical_real'
            formula = 'lambda magnitude from tau_emp(Re(phi)), phase preserved from EDMD';
        case 'empirical_complex'
            formula = 'lambda fitted directly from complex phi_t ~ lambda * phi_{t-1}';
        otherwise
            formula = lower(cfg.lambda_source);
    end
end

function formula = local_method_formula(cfg)
% Text description of the perturbation definition used in this run.
    switch lower(cfg.method)
        case 'wiener'
            formula = 'phi = h_lambda * u (Wiener deconvolution)';
        case {'koopman_residual', 'residual'}
            switch lower(cfg.first_u_mode)
                case 'phi1'
                    formula = 'u(1)=phi(1); u(t)=phi(t)-lambda_d*phi(t-1), t>=2';
                case 'zero'
                    formula = 'u(1)=0; u(t)=phi(t)-lambda_d*phi(t-1), t>=2';
                otherwise
                    formula = 'u(t)=phi(t)-lambda_d*phi(t-1)';
            end
        otherwise
            formula = lower(cfg.method);
    end
end
