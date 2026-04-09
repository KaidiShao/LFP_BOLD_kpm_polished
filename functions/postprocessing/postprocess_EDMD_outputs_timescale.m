function [fig, timescale_info] = postprocess_EDMD_outputs_timescale(EDMD_outputs, cfg)
%POSTPROCESS_EDMD_OUTPUTS_TIMESCALE
%   Create a 4x4 diagnostic figure for Koopman eigenfunctions:
%     Row 1: correlation-based similarity matrices (heatmaps)
%     Row 2: theoretical eigenfunction decay implied by eigenvalues (lambda_d^n)
%     Row 3: empirical autocorrelation functions (ACF) of eigenfunctions
%     Row 4: theory-vs-empirical timescale error per eigenfunction
%
%   The 4 columns correspond to:
%     (1) ALL normalized eigenfunctions (abs)
%     (2) ALL normalized eigenfunctions (real)
%     (3) SELECTED normalized eigenfunctions (abs)
%     (4) SELECTED normalized eigenfunctions (real)
%
%   Plot styling and coloring follows plot_acf_and_theoretical.m:
%     - Spectral colormap (othercolor('Spectral10')) when available
%     - Overlaid curves (no legend)
%     - Column-wise normalisation uses only division by max(abs), without de-meaning
%
%   Timescale comparison (Row 4):
%     - Theoretical: kappa_th = -Re(log(lambda_d)/dt), tau_th = 1/kappa_th
%     - Empirical: envelope ACF |rho|, log-linear fit: log(|rho|) ~ -t/tau_emp + c
%     - Error metric (cfg.err_metric):
%         'logtau' (default): |log(tau_emp) - log(tau_th)|
%         'kappa'            : |kappa_emp - kappa_th|, where kappa_emp = 1/tau_emp
%
%   INPUT
%     EDMD_outputs
%       Required:
%         .evalues   [Ksel x 1]  discrete eigenvalues for selected set
%         .efuns     [T x Ksel]  eigenfunctions (time x mode) for selected set
%       Optional (recommended):
%         .original_sorted.evalues [Kall x 1]
%         .original_sorted.efuns   [T x Kall]
%
%     cfg (optional struct)
%       .dt               sampling interval (default 1)
%       .t_plot           indices for similarity computation (default 1:T)
%       .maxLag           max lag in samples for ACF (default min(200, T-1))
%       .xlim_time        x-axis limit in time units for Row 2/3 (default 20)
%       .max_modes_all    cap #modes for the ALL columns (default 80)
%       .max_modes_sel    cap #modes for the SELECTED columns (default Inf)
%       .sim_clim         color limits for similarity heatmaps (default [-1 1])
%       .title_prefix     optional figure title (default '')
%
%       ACF-envelope fit settings:
%       .fit_env_start    start threshold on |rho| (default 0.9)
%       .fit_env_end      end threshold on |rho| (default 0.2)
%       .fit_env_floor    ignore points below this (default 1e-3)
%       .fit_min_points   minimum points for fit (default 8)
%       .fit_fallback     'efold' (default) or 'nan'
%
%       Comparison:
%       .err_metric       'logtau' (default) or 'kappa'
%       .empirical_lambda_phase_mode  'preserve_edmd_angle' (default) or 'real_positive'
%
%   OUTPUT
%     fig  figure handle (results are also stored in fig.UserData.cond)

    if nargin < 2 || isempty(cfg)
        cfg = struct();
    end

    % -------------------- defaults --------------------
    if ~isfield(cfg, 'dt') || isempty(cfg.dt), cfg.dt = 1; end
    if ~isfield(cfg, 't_plot'), cfg.t_plot = []; end
    if ~isfield(cfg, 'maxLag') || isempty(cfg.maxLag), cfg.maxLag = []; end
    if ~isfield(cfg, 'xlim_time') || isempty(cfg.xlim_time), cfg.xlim_time = 20; end
    if ~isfield(cfg, 'max_modes_all') || isempty(cfg.max_modes_all), cfg.max_modes_all = 80; end
    if ~isfield(cfg, 'max_modes_sel') || isempty(cfg.max_modes_sel), cfg.max_modes_sel = inf; end
    if ~isfield(cfg, 'sim_clim') || isempty(cfg.sim_clim), cfg.sim_clim = [-1 1]; end
    if ~isfield(cfg, 'title_prefix') || isempty(cfg.title_prefix), cfg.title_prefix = ''; end

    if ~isfield(cfg, 'fit_env_start') || isempty(cfg.fit_env_start), cfg.fit_env_start = 0.9; end
    if ~isfield(cfg, 'fit_env_end')   || isempty(cfg.fit_env_end),   cfg.fit_env_end   = 0.2; end
    if ~isfield(cfg, 'fit_env_floor') || isempty(cfg.fit_env_floor), cfg.fit_env_floor = 1e-3; end
    if ~isfield(cfg, 'fit_min_points')|| isempty(cfg.fit_min_points),cfg.fit_min_points = 8; end
    if ~isfield(cfg, 'fit_fallback')  || isempty(cfg.fit_fallback),  cfg.fit_fallback  = 'efold'; end

    if ~isfield(cfg, 'err_metric')    || isempty(cfg.err_metric),    cfg.err_metric    = 'logtau'; end
    if ~isfield(cfg, 'empirical_lambda_phase_mode') || isempty(cfg.empirical_lambda_phase_mode)
        cfg.empirical_lambda_phase_mode = 'preserve_edmd_angle';
    end

    dt = cfg.dt;

    % -------------------- fetch ALL vs SELECTED --------------------
    efuns_sel_raw = EDMD_outputs.efuns;
    evals_sel_raw = EDMD_outputs.evalues;

    if isfield(EDMD_outputs, 'original_sorted') && ...
            isfield(EDMD_outputs.original_sorted, 'efuns') && ...
            isfield(EDMD_outputs.original_sorted, 'evalues')
        efuns_all_raw = EDMD_outputs.original_sorted.efuns;
        evals_all_raw = EDMD_outputs.original_sorted.evalues;
    else
        efuns_all_raw = efuns_sel_raw;
        evals_all_raw = evals_sel_raw;
        warning('postprocess_EDMD_outputs_timescale:NoOriginalSorted', ...
            'No EDMD_outputs.original_sorted.* found. Using selected set as ALL.');
    end

    % -------------------- cap modes and align sizes --------------------
    [Tall, Kall0] = size(efuns_all_raw);
    [Tsel, Ksel0] = size(efuns_sel_raw);

    Kall = min([Kall0, numel(evals_all_raw), cfg.max_modes_all]);
    Ksel = min([Ksel0, numel(evals_sel_raw), cfg.max_modes_sel]);

    efuns_all_raw = efuns_all_raw(:, 1:Kall);
    evals_all_raw = evals_all_raw(1:Kall);

    efuns_sel_raw = efuns_sel_raw(:, 1:Ksel);
    evals_sel_raw = evals_sel_raw(1:Ksel);

    % -------------------- time windows --------------------
    if isempty(cfg.t_plot)
        t_plot_all = 1:Tall;
        t_plot_sel = 1:Tsel;
    else
        t_plot_all = cfg.t_plot(cfg.t_plot >= 1 & cfg.t_plot <= Tall);
        t_plot_sel = cfg.t_plot(cfg.t_plot >= 1 & cfg.t_plot <= Tsel);
        if isempty(t_plot_all), t_plot_all = 1:Tall; end
        if isempty(t_plot_sel), t_plot_sel = 1:Tsel; end
    end

    if isempty(cfg.maxLag)
        maxLag_all = min(200, Tall-1);
        maxLag_sel = min(200, Tsel-1);
    else
        maxLag_all = min(cfg.maxLag, Tall-1);
        maxLag_sel = min(cfg.maxLag, Tsel-1);
    end

    % -------------------- prepare the 4 conditions --------------------
    cond(1).name = 'ALL |\phi| (norm)';
    cond(1).efuns_raw = efuns_all_raw;
    cond(1).evals = evals_all_raw;
    cond(1).feat = 'abs';
    cond(1).t_plot = t_plot_all;
    cond(1).maxLag = maxLag_all;

    cond(2).name = 'ALL Re(\phi) (norm)';
    cond(2).efuns_raw = efuns_all_raw;
    cond(2).evals = evals_all_raw;
    cond(2).feat = 'real';
    cond(2).t_plot = t_plot_all;
    cond(2).maxLag = maxLag_all;

    cond(3).name = 'SELECTED |\phi| (norm)';
    cond(3).efuns_raw = efuns_sel_raw;
    cond(3).evals = evals_sel_raw;
    cond(3).feat = 'abs';
    cond(3).t_plot = t_plot_sel;
    cond(3).maxLag = maxLag_sel;

    cond(4).name = 'SELECTED Re(\phi) (norm)';
    cond(4).efuns_raw = efuns_sel_raw;
    cond(4).evals = evals_sel_raw;
    cond(4).feat = 'real';
    cond(4).t_plot = t_plot_sel;
    cond(4).maxLag = maxLag_sel;

    for c = 1:4
        % Normalize eigenfunctions per feature type (abs/real)
        cond(c).Phi = local_normalize_efuns(cond(c).efuns_raw, cond(c).feat); % [T x K]
        cond(c).K = size(cond(c).Phi, 2);
        cond(c).T = size(cond(c).Phi, 1);
        cond(c).t = (0:cond(c).T-1) * dt;

        % Colormap (per-eigenfunction color)
        cond(c).cmap = local_spectral_cmap(cond(c).K);

        % Similarity matrix uses only the specified window
        X = cond(c).Phi(cond(c).t_plot, :);
        cond(c).R = local_corrcoef_safe(X);

        % Theoretical timescales
        [cond(c).tau_th, cond(c).kappa_th] = local_theoretical_timescale(cond(c).evals, dt);

        % Empirical ACF, timescales, and lambda estimates from eigenfunction dynamics
        fit_cfg = cfg;
        fit_cfg.feat = cond(c).feat;
        fit_cfg.maxLag = cond(c).maxLag;
        fit_cfg.lambda_phase_mode = cfg.empirical_lambda_phase_mode;
        cond(c).empirical = fit_empirical_eigenvalues_from_efuns(cond(c).efuns_raw, cond(c).evals, fit_cfg);
        cond(c).acf_lags = cond(c).empirical.acf_lags;
        cond(c).acf_mat = cond(c).empirical.acf_mat;
        cond(c).tau_emp = cond(c).empirical.tau_emp;
        cond(c).kappa_emp = cond(c).empirical.kappa_emp;
        cond(c).lambda_emp = cond(c).empirical.lambda_emp_discrete;
        cond(c).lambda_emp_safe = cond(c).empirical.lambda_emp_discrete_safe;
        cond(c).lambda_emp_complex = cond(c).empirical.lambda_emp_complex_discrete;
        cond(c).lambda_emp_complex_safe = cond(c).empirical.lambda_emp_complex_discrete_safe;
        cond(c).lambda_emp_complex_tau = cond(c).empirical.lambda_emp_complex_tau;
        cond(c).lambda_emp_complex_kappa = cond(c).empirical.lambda_emp_complex_kappa;
        cond(c).fit_info = cond(c).empirical.fit_info;

        % Error metric
        cond(c).err = local_timescale_error(cond(c).tau_th, cond(c).kappa_th, ...
                                            cond(c).tau_emp, cfg.err_metric);
    end

    % -------------------- plot --------------------
    fig = figure('Color', 'w', 'Name', 'EDMD timescales & similarity');
    tl = tiledlayout(fig, 4, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

    if ~isempty(cfg.title_prefix)
        title(tl, cfg.title_prefix, 'Interpreter', 'none');
    end

    % Row 1: similarity matrices
    for c = 1:4
        nexttile(tl, c);
        imagesc(cond(c).R);
        axis image;
        set(gca, 'YDir', 'normal');
        clim(cfg.sim_clim);
        colormap(gca, local_diverging_cmap());
        colorbar;
        title({cond(c).name, 'Similarity: corrcoef(\phi)'}, 'FontWeight', 'normal');
        xlabel('mode #');
        ylabel('mode #');
    end

    % Row 2: theoretical eigenfunctions (lambda_d^n), normalized over full T (as in plot_acf_and_theoretical)
    for c = 1:4
        nexttile(tl, 4 + c);
        hold on;

        if cond(c).K == 0
            title({'Theoretical: \lambda_d^n', cond(c).name}, 'FontWeight', 'normal');
            hold off;
            continue;
        end

        n_plot = min(cond(c).T, floor(cfg.xlim_time / dt) + 1);
        t_plot = cond(c).t(1:n_plot);

        for i = 1:cond(c).K
            ld = cond(c).evals(i);

            % Compute normalized theoretical curve; normalization uses full T.
            x_plot = local_theoretical_curve_plot(ld, cond(c).feat, cond(c).T, n_plot);

            plot(t_plot, x_plot, 'Color', cond(c).cmap(i,:), 'LineWidth', 0.8);
        end

        grid on;
        xlim([0 cfg.xlim_time]);
        xlabel('time');
        ylabel('amplitude (normalised)');
        title({'Theoretical: \lambda_d^n', cond(c).name}, 'FontWeight', 'normal');
        hold off;
    end

    % Row 3: ACF curves
    for c = 1:4
        nexttile(tl, 8 + c);
        hold on;

        if cond(c).K == 0
            title({'Empirical: ACF(\phi)', cond(c).name}, 'FontWeight', 'normal');
            hold off;
            continue;
        end

        for i = 1:cond(c).K
            plot(cond(c).acf_lags * dt, cond(c).acf_mat(:, i), ...
                'Color', cond(c).cmap(i,:), 'LineWidth', 0.8);
        end

        grid on;
        xlim([0 cfg.xlim_time]);
        xlabel('lag (time)');
        ylabel('ACF');
        title({'Empirical: ACF(\phi)', cond(c).name}, 'FontWeight', 'normal');
        hold off;
    end

    % Row 4: timescale error per eigenfunction index
    for c = 1:4
        nexttile(tl, 12 + c);
        hold on;

        if cond(c).K == 0
            title({'Timescale error', cond(c).name}, 'FontWeight', 'normal');
            hold off;
            continue;
        end

        if strcmpi(cond(c).feat, 'abs')
            mk = '*';
        else
            mk = '+';
        end

        for i = 1:cond(c).K
            if ~isfinite(cond(c).err(i))
                continue;
            end
            plot(i, cond(c).err(i), mk, ...
                'Color', cond(c).cmap(i,:), 'LineWidth', 1.0, 'MarkerSize', 6);
        end

        grid on;
        xlim([0 cond(c).K + 1]);
        xlabel('mode #');
        if strcmpi(cfg.err_metric, 'kappa')
            ylabel('| \kappa_{emp} - \kappa_{th} |');
            title({'Timescale error (rate)', cond(c).name}, 'FontWeight', 'normal');
        else
            ylabel('| log(\tau_{emp}) - log(\tau_{th}) |');
            title({'Timescale error (log-\tau)', cond(c).name}, 'FontWeight', 'normal');
        end
        hold off;
    end

    % Store results for programmatic access
    timescale_info = local_pack_timescale_info(cond, cfg);
    fig.UserData.cond = cond;
    fig.UserData.cfg = cfg;
    fig.UserData.timescale_info = timescale_info;
end

% =====================================================================
% Helpers
% =====================================================================

function Phi = local_normalize_efuns(efuns_raw, feat)
%LOCAL_NORMALIZE_EFUNS Convert (abs/real) and normalize each column.
%   Uses the shared project helper when available so normalization stays
%   consistent across EDMD postprocessing scripts.
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

function R = local_corrcoef_safe(X)
%LOCAL_CORRCOEF_SAFE corrcoef with NaN handling.
    if isempty(X)
        R = 0;
        return;
    end
    try
        R = corrcoef(X);
    catch
        % Fallback: compute with normalization
        X = X - mean(X, 1, 'omitnan');
        denom = sqrt(sum(X.^2, 1, 'omitnan'));
        denom(denom == 0) = 1;
        Xn = X ./ denom;
        R = Xn' * Xn;
    end
    R(~isfinite(R)) = 0;
    R(1:size(R,1)+1:end) = 1;
end

function cmap = local_spectral_cmap(n)
%LOCAL_SPECTRAL_CMAP Use othercolor('Spectral10') if available.
    if n <= 0
        cmap = zeros(0, 3);
        return;
    end

    if exist('othercolor', 'file') == 2
        base = othercolor('Spectral10');
    else
        base = lines(10);
    end

    if size(base, 1) == n
        cmap = base;
        return;
    end

    % Interpolate to desired size
    x0 = linspace(0, 1, size(base, 1));
    x1 = linspace(0, 1, n);
    cmap = interp1(x0, base, x1);
end

function cmap = local_diverging_cmap()
%LOCAL_DIVERGING_CMAP Diverging colormap for correlation heatmaps.
    if exist('othercolor', 'file') == 2
        cmap = othercolor('RdBu11');
        cmap = flipud(cmap);
    else
        cmap = jet(256);
    end
end

function [tau_th, kappa_th] = local_theoretical_timescale(lambda_d, dt)
%LOCAL_THEORETICAL_TIMESCALE kappa_th = -Re(log(lambda_d)/dt), tau_th = 1/kappa_th.
    lambda_d = lambda_d(:);
    alpha = log(lambda_d) / dt;
    kappa_th = -real(alpha);
    tau_th = 1 ./ kappa_th;

    tau_th(~isfinite(tau_th)) = inf;
    tau_th(kappa_th <= 0) = inf;
    kappa_th(~isfinite(kappa_th) | kappa_th < 0) = 0;
end

function timescale_info = local_pack_timescale_info(cond, cfg)
%LOCAL_PACK_TIMESCALE_INFO Create a convenient structure for downstream use.
    timescale_info = struct();
    timescale_info.cfg = cfg;
    timescale_info.all = struct();
    timescale_info.all.abs = local_make_timescale_summary(cond(1));
    timescale_info.all.real = local_make_timescale_summary(cond(2));

    timescale_info.selected = struct();
    timescale_info.selected.abs = local_make_timescale_summary(cond(3));
    timescale_info.selected.real = local_make_timescale_summary(cond(4));
end

function summary = local_make_timescale_summary(cond_entry)
%LOCAL_MAKE_TIMESCALE_SUMMARY Keep only compact downstream-relevant fields.
    summary = struct();
    summary.name = cond_entry.name;
    summary.feat = cond_entry.feat;
    summary.evals = cond_entry.evals;
    summary.tau_th = cond_entry.tau_th;
    summary.kappa_th = cond_entry.kappa_th;
    summary.tau_emp = cond_entry.tau_emp;
    summary.kappa_emp = cond_entry.kappa_emp;
    summary.lambda_emp = cond_entry.lambda_emp;
    summary.lambda_emp_safe = cond_entry.lambda_emp_safe;
    summary.lambda_emp_complex = cond_entry.lambda_emp_complex;
    summary.lambda_emp_complex_safe = cond_entry.lambda_emp_complex_safe;
    summary.lambda_emp_complex_tau = cond_entry.lambda_emp_complex_tau;
    summary.lambda_emp_complex_kappa = cond_entry.lambda_emp_complex_kappa;
    summary.err = cond_entry.err;
    summary.fit_info = cond_entry.fit_info;
    summary.empirical = cond_entry.empirical;
end

function err = local_timescale_error(tau_th, kappa_th, tau_emp, err_metric)
%LOCAL_TIMESCALE_ERROR Compute per-mode error.
    tau_th  = tau_th(:);
    kappa_th = kappa_th(:);
    tau_emp = tau_emp(:);

    err = nan(size(tau_th));
    if strcmpi(err_metric, 'kappa')
        kappa_emp = 1 ./ tau_emp;
        good = isfinite(kappa_emp) & isfinite(kappa_th) & (kappa_emp >= 0) & (kappa_th >= 0);
        err(good) = abs(kappa_emp(good) - kappa_th(good));
    else
        good = isfinite(tau_emp) & isfinite(tau_th) & (tau_emp > 0) & (tau_th > 0);
        err(good) = abs(log(tau_emp(good)) - log(tau_th(good)));
    end
end

function x_plot = local_theoretical_curve_plot(ld, feat, T, n_plot)
%LOCAL_THEORETICAL_CURVE_PLOT
%   Compute normalized theoretical curve for plotting:
%     x_th(n) = abs(ld^n) or real(ld^n), n=0..T-1
%   Normalization matches the project EDMD visualisation helpers:
%     x_th <- x_th / max(abs(x_th)) over full length T
%   Returns only the first n_plot samples (for plotting), but normalization
%   is computed using the full length T.

    n_plot = min(max(1, n_plot), T);

    % First pass: max abs over full T (no large allocation)
    z = 1;
    mx = 0;
    for k = 1:T
        if strcmpi(feat, 'real')
            raw = real(z);
        else
            raw = abs(z);
        end
        mx = max(mx, abs(raw));
        z = z * ld;
    end

    % Second pass: cache first n_plot samples
    z = 1;
    x_cache = zeros(n_plot, 1);
    for k = 1:T
        if strcmpi(feat, 'real')
            raw = real(z);
        else
            raw = abs(z);
        end
        if k <= n_plot
            x_cache(k) = raw;
        end
        z = z * ld;
    end

    if mx > 0
        x_plot = x_cache / mx;
    else
        x_plot = x_cache;
    end
end
