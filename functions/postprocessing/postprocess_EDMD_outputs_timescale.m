function [fig, timescale_info] = postprocess_EDMD_outputs_timescale(EDMD_outputs, cfg)
%POSTPROCESS_EDMD_OUTPUTS_TIMESCALE
%   Create a 6x4 diagnostic figure for Koopman eigenfunctions:
%     Row 1: correlation-based similarity matrices (heatmaps)
%     Row 2: theoretical eigenfunction decay implied by eigenvalues (lambda_d^n)
%     Row 3: empirical eigenfunction decay implied by fitted eigenvalues
%     Row 4: empirical autocorrelation functions (ACF) of eigenfunctions
%     Row 5: theory-vs-empirical lambda-fit timescale error
%     Row 6: theory-vs-empirical ACF-envelope timescale error
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
%   Timescale comparison:
%     - Theoretical: kappa_th = -Re(log(lambda_d)/dt), tau_th = 1/kappa_th
%     - Empirical lambda fit: fit lambda_emp from phi(t+1) ~= lambda_emp * phi(t)
%     - Empirical ACF envelope: log(|rho|) ~ -t/tau_emp + c
%     - Error metric (cfg.err_metric):
%         'logtau' (default): |log(tau_emp) - log(tau_th)|
%         'kappa'            : |kappa_emp - kappa_th|
%
%   INPUT
%     EDMD_outputs
%       Required, either:
%         .evalues   [Ksel x 1]  discrete eigenvalues for selected set
%         .efuns     [T x Ksel]  eigenfunctions (time x mode)
%       Or, for large runs:
%         .evalues
%         .source_efuns [T x Ksource]
%         .efun_col_indices [Ksel x 1]
%       Optional (recommended): .original_sorted with the same matrix or
%       source_efuns/efun_col_indices convention for the ALL columns.
%
%     cfg (optional struct)
%       .dt               sampling interval (default 1)
%       .t_plot           indices for similarity computation (default 1:T)
%       .maxLag           max lag in samples for ACF (default covers xlim_time)
%       .xlim_time        x-axis limit in time units for Row 2/3/4 (default 20)
%       .max_modes_all    cap #modes for the ALL columns (default Inf)
%       .max_modes_sel    cap #modes for the SELECTED columns (default Inf)
%       .sim_clim         color limits for similarity heatmaps (default [-1 1])
%       .title_prefix     optional figure title (default '')
%       .match_empirical_scale_to_theoretical  default true
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
    if ~isfield(cfg, 'max_modes_all') || isempty(cfg.max_modes_all), cfg.max_modes_all = inf; end
    if ~isfield(cfg, 'max_modes_sel') || isempty(cfg.max_modes_sel), cfg.max_modes_sel = inf; end
    if ~isfield(cfg, 'sim_clim') || isempty(cfg.sim_clim), cfg.sim_clim = [-1 1]; end
    if ~isfield(cfg, 'title_prefix') || isempty(cfg.title_prefix), cfg.title_prefix = ''; end
    if ~isfield(cfg, 'match_empirical_scale_to_theoretical') || ...
            isempty(cfg.match_empirical_scale_to_theoretical)
        cfg.match_empirical_scale_to_theoretical = true;
    end

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
    sel_set = local_resolve_selected_mode_set(EDMD_outputs);
    if isfield(EDMD_outputs, 'original_sorted')
        all_set = local_resolve_original_mode_set( ...
            EDMD_outputs.original_sorted, sel_set, EDMD_outputs);
    else
        all_set = sel_set;
        warning('postprocess_EDMD_outputs_timescale:NoOriginalSorted', ...
            'No EDMD_outputs.original_sorted.* found. Using selected set as ALL.');
    end

    % -------------------- cap modes and align sizes --------------------
    all_set = local_cap_mode_set(all_set, cfg.max_modes_all);
    sel_set = local_cap_mode_set(sel_set, cfg.max_modes_sel);
    all_label = local_all_mode_label(EDMD_outputs);

    Tall = all_set.T;
    Tsel = sel_set.T;

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

    xlim_lag = max(1, ceil(cfg.xlim_time / dt));
    if isempty(cfg.maxLag) || ~isfinite(cfg.maxLag)
        if cfg.match_empirical_scale_to_theoretical
            maxLag_requested = xlim_lag;
        else
            maxLag_requested = 200;
        end
    else
        maxLag_requested = cfg.maxLag;
        if cfg.match_empirical_scale_to_theoretical
            maxLag_requested = max(maxLag_requested, xlim_lag);
        end
    end

    maxLag_all = min(maxLag_requested, Tall-1);
    maxLag_sel = min(maxLag_requested, Tsel-1);

    % -------------------- prepare the 4 conditions --------------------
    cond(1).name = [all_label ' |\phi| (norm)'];
    cond(1).source_efuns = all_set.source_efuns;
    cond(1).col_indices = all_set.col_indices;
    cond(1).evals = all_set.evals;
    cond(1).feat = 'abs';
    cond(1).t_plot = t_plot_all;
    cond(1).maxLag = maxLag_all;

    cond(2).name = [all_label ' Re(\phi) (norm)'];
    cond(2).source_efuns = all_set.source_efuns;
    cond(2).col_indices = all_set.col_indices;
    cond(2).evals = all_set.evals;
    cond(2).feat = 'real';
    cond(2).t_plot = t_plot_all;
    cond(2).maxLag = maxLag_all;

    cond(3).name = 'SELECTED |\phi| (norm)';
    cond(3).source_efuns = sel_set.source_efuns;
    cond(3).col_indices = sel_set.col_indices;
    cond(3).evals = sel_set.evals;
    cond(3).feat = 'abs';
    cond(3).t_plot = t_plot_sel;
    cond(3).maxLag = maxLag_sel;

    cond(4).name = 'SELECTED Re(\phi) (norm)';
    cond(4).source_efuns = sel_set.source_efuns;
    cond(4).col_indices = sel_set.col_indices;
    cond(4).evals = sel_set.evals;
    cond(4).feat = 'real';
    cond(4).t_plot = t_plot_sel;
    cond(4).maxLag = maxLag_sel;

    for c = 1:4
        cond(c).K = numel(cond(c).col_indices);
        cond(c).T = size(cond(c).source_efuns, 1);
        cond(c).t = (0:cond(c).T-1) * dt;

        % Colormap (per-eigenfunction color)
        cond(c).cmap = local_spectral_cmap(cond(c).K);

        % Similarity matrix uses only the specified window
        X = local_normalize_efuns_from_source( ...
            cond(c).source_efuns, cond(c).col_indices, ...
            cond(c).t_plot, cond(c).feat);
        cond(c).R = local_corrcoef_safe(X);
        clear X;

        % Theoretical timescales
        [cond(c).tau_th, cond(c).kappa_th] = local_theoretical_timescale(cond(c).evals, dt);

        % Empirical ACF, timescales, and lambda estimates from eigenfunction dynamics
        fit_cfg = cfg;
        fit_cfg.feat = cond(c).feat;
        fit_cfg.maxLag = cond(c).maxLag;
        fit_cfg.lambda_phase_mode = cfg.empirical_lambda_phase_mode;
        fit_input = struct();
        fit_input.source_efuns = cond(c).source_efuns;
        fit_input.col_indices = cond(c).col_indices;
        cond(c).empirical = fit_empirical_eigenvalues_from_efuns( ...
            fit_input, cond(c).evals, fit_cfg);
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
        cond(c).lambda_emp_complex_valid = cond(c).empirical.lambda_emp_complex_valid;
        cond(c).fit_info = cond(c).empirical.fit_info;

        cond(c).tau_emp_lambda_fit = cond(c).lambda_emp_complex_tau;
        cond(c).kappa_emp_lambda_fit = cond(c).lambda_emp_complex_kappa;
        cond(c).tau_emp_lambda_fit(~cond(c).lambda_emp_complex_valid) = nan;
        cond(c).kappa_emp_lambda_fit(~cond(c).lambda_emp_complex_valid) = nan;

        cond(c).tau_emp_acf_envelope = cond(c).tau_emp;
        cond(c).kappa_emp_acf_envelope = cond(c).kappa_emp;

        cond(c).err_lambda_fit = local_timescale_error( ...
            cond(c).tau_th, cond(c).kappa_th, ...
            cond(c).tau_emp_lambda_fit, cond(c).kappa_emp_lambda_fit, cfg.err_metric);
        cond(c).err_acf_envelope = local_timescale_error( ...
            cond(c).tau_th, cond(c).kappa_th, ...
            cond(c).tau_emp_acf_envelope, cond(c).kappa_emp_acf_envelope, cfg.err_metric);
    end

    % -------------------- plot --------------------
    fig = figure('Color', 'w', 'Name', 'EDMD timescales & similarity', ...
        'Position', [100, 100, 2400, 3000]);
    tl = tiledlayout(fig, 6, 4, 'TileSpacing', 'compact', 'Padding', 'compact');

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

    % Row 3: empirical lambda-fit curves
    for c = 1:4
        nexttile(tl, 8 + c);
        hold on;

        if cond(c).K == 0
            title({'Empirical fit: \lambda_{emp}^n', cond(c).name}, 'FontWeight', 'normal');
            hold off;
            continue;
        end

        n_plot = min(cond(c).T, floor(cfg.xlim_time / dt) + 1);
        t_plot = cond(c).t(1:n_plot);
        for i = 1:cond(c).K
            if ~cond(c).lambda_emp_complex_valid(i)
                continue;
            end
            ld = cond(c).lambda_emp_complex(i);
            x_plot = local_theoretical_curve_plot(ld, cond(c).feat, cond(c).T, n_plot);
            plot(t_plot, x_plot, 'Color', cond(c).cmap(i,:), 'LineWidth', 0.8);
        end

        grid on;
        xlim([0 cfg.xlim_time]);
        xlabel('time');
        ylabel('amplitude (normalised)');
        title({'Empirical fit: \lambda_{emp}^n', cond(c).name}, 'FontWeight', 'normal');
        hold off;
    end

    % Row 4: empirical ACF curves
    for c = 1:4
        nexttile(tl, 12 + c);
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

    % Row 5: lambda-fit timescale error per eigenfunction index
    for c = 1:4
        nexttile(tl, 16 + c);
        local_plot_error_points(cond(c), cond(c).err_lambda_fit, cfg, ...
            'Timescale error: lambda-fit');
    end

    % Row 6: ACF-envelope timescale error per eigenfunction index
    for c = 1:4
        nexttile(tl, 20 + c);
        local_plot_error_points(cond(c), cond(c).err_acf_envelope, cfg, ...
            'Timescale error: ACF-envelope');
    end

    % Store results for programmatic access
    timescale_info = local_pack_timescale_info(cond, cfg);
    fig.UserData.cond = local_strip_source_efuns_from_cond(cond);
    fig.UserData.cfg = cfg;
    fig.UserData.timescale_info = timescale_info;
end

% =====================================================================
% Helpers
% =====================================================================

function mode_set = local_resolve_selected_mode_set(EDMD_outputs)
%LOCAL_RESOLVE_SELECTED_MODE_SET Resolve selected modes without copying columns.
    if isfield(EDMD_outputs, 'source_efuns') && ...
            isfield(EDMD_outputs, 'efun_col_indices') && ...
            ~isempty(EDMD_outputs.source_efuns)
        source_efuns = EDMD_outputs.source_efuns;
        col_indices = EDMD_outputs.efun_col_indices(:);
    elseif isfield(EDMD_outputs, 'efuns') && ~isempty(EDMD_outputs.efuns)
        source_efuns = EDMD_outputs.efuns;
        col_indices = (1:size(source_efuns, 2)).';
    else
        error('EDMD_outputs must contain efuns or source_efuns + efun_col_indices.');
    end

    mode_set = local_make_mode_set(source_efuns, col_indices, EDMD_outputs.evalues);
end

function mode_set = local_resolve_original_mode_set(original_sorted, fallback_set, EDMD_outputs)
%LOCAL_RESOLVE_ORIGINAL_MODE_SET Resolve ALL modes without copying columns.
    if isfield(original_sorted, 'source_efuns') && ...
            isfield(original_sorted, 'efun_col_indices') && ...
            ~isempty(original_sorted.source_efuns)
        source_efuns = original_sorted.source_efuns;
        col_indices = original_sorted.efun_col_indices(:);
        mode_set = local_make_mode_set(source_efuns, col_indices, original_sorted.evalues);
    elseif isfield(EDMD_outputs, 'source_efuns') && ...
            isfield(original_sorted, 'efun_col_indices') && ...
            ~isempty(EDMD_outputs.source_efuns)
        source_efuns = EDMD_outputs.source_efuns;
        col_indices = original_sorted.efun_col_indices(:);
        mode_set = local_make_mode_set(source_efuns, col_indices, original_sorted.evalues);
    elseif isfield(original_sorted, 'efuns') && ~isempty(original_sorted.efuns)
        source_efuns = original_sorted.efuns;
        col_indices = (1:size(source_efuns, 2)).';
        mode_set = local_make_mode_set(source_efuns, col_indices, original_sorted.evalues);
    else
        mode_set = fallback_set;
        warning('postprocess_EDMD_outputs_timescale:NoOriginalSortedEfuns', ...
            'No original_sorted eigenfunction source found. Using selected set as ALL.');
    end
end

function mode_set = local_make_mode_set(source_efuns, col_indices, evals)
    evals = evals(:);
    col_indices = col_indices(:);
    n_source_cols = size(source_efuns, 2);
    keep = isfinite(col_indices) & col_indices >= 1 & col_indices <= n_source_cols;
    col_indices = col_indices(keep);
    n = min(numel(col_indices), numel(evals));

    mode_set = struct();
    mode_set.source_efuns = source_efuns;
    mode_set.col_indices = col_indices(1:n);
    mode_set.evals = evals(1:n);
    mode_set.T = size(source_efuns, 1);
    mode_set.K = n;
end

function mode_set = local_cap_mode_set(mode_set, max_modes)
    if isfinite(max_modes)
        n = min(mode_set.K, max(0, floor(max_modes)));
    else
        n = mode_set.K;
    end

    mode_set.col_indices = mode_set.col_indices(1:n);
    mode_set.evals = mode_set.evals(1:n);
    mode_set.K = n;
end

function label = local_all_mode_label(EDMD_outputs)
label = 'ALL';
if isfield(EDMD_outputs, 'original_sorted') && ...
        isfield(EDMD_outputs.original_sorted, 'abs_thresh') && ...
        ~isempty(EDMD_outputs.original_sorted.abs_thresh)
    label = sprintf('THRESHOLDED ALL |\\lambda|>%g', ...
        EDMD_outputs.original_sorted.abs_thresh);
end
end

function X = local_normalize_efuns_from_source(source_efuns, col_indices, rows, feat)
%LOCAL_NORMALIZE_EFUNS_FROM_SOURCE Normalize one column at a time.
    rows = rows(:);
    col_indices = col_indices(:);
    X = zeros(numel(rows), numel(col_indices));
    for j = 1:numel(col_indices)
        x = source_efuns(rows, col_indices(j));
        if strcmpi(feat, 'real')
            x = real(x);
        else
            x = abs(x);
        end
        mx = max(abs(x));
        if mx > 0
            x = x / mx;
        end
        X(:, j) = x;
    end
end

function cond_user = local_strip_source_efuns_from_cond(cond)
%LOCAL_STRIP_SOURCE_EFUNS_FROM_COND Keep figure UserData compact.
    cond_user = cond;
    for i = 1:numel(cond_user)
        cond_user(i).source_size = size(cond_user(i).source_efuns);
        if isfield(cond_user, 'source_efuns')
            cond_user(i).source_efuns = [];
        end
    end
end

function local_plot_error_points(cond_entry, err, cfg, title_prefix)
%LOCAL_PLOT_ERROR_POINTS Plot per-mode timescale error for one empirical definition.
    hold on;

    if cond_entry.K == 0
        title({title_prefix, cond_entry.name}, 'FontWeight', 'normal');
        hold off;
        return;
    end

    if strcmpi(cond_entry.feat, 'abs')
        mk = '*';
    else
        mk = '+';
    end

    for i = 1:cond_entry.K
        if ~isfinite(err(i))
            continue;
        end
        plot(i, err(i), mk, ...
            'Color', cond_entry.cmap(i,:), 'LineWidth', 1.0, 'MarkerSize', 6);
    end

    grid on;
    xlim([0 cond_entry.K + 1]);
    xlabel('mode #');
    if strcmpi(cfg.err_metric, 'kappa')
        ylabel('| \kappa_{emp} - \kappa_{th} |');
        title({[title_prefix ' (rate)'], cond_entry.name}, 'FontWeight', 'normal');
    else
        ylabel('| log(\tau_{emp}) - log(\tau_{th}) |');
        title({[title_prefix ' (log-\tau)'], cond_entry.name}, 'FontWeight', 'normal');
    end
    hold off;
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
    summary.lambda_emp_complex_valid = cond_entry.lambda_emp_complex_valid;
    summary.tau_emp_lambda_fit = cond_entry.tau_emp_lambda_fit;
    summary.kappa_emp_lambda_fit = cond_entry.kappa_emp_lambda_fit;
    summary.tau_emp_acf_envelope = cond_entry.tau_emp_acf_envelope;
    summary.kappa_emp_acf_envelope = cond_entry.kappa_emp_acf_envelope;
    summary.err_lambda_fit = cond_entry.err_lambda_fit;
    summary.err_acf_envelope = cond_entry.err_acf_envelope;
    summary.err = cond_entry.err_lambda_fit;
    summary.fit_info = cond_entry.fit_info;
    summary.empirical = cond_entry.empirical;
end

function err = local_timescale_error(tau_th, kappa_th, tau_emp, kappa_emp, err_metric)
%LOCAL_TIMESCALE_ERROR Compute per-mode error.
    tau_th  = tau_th(:);
    kappa_th = kappa_th(:);
    tau_emp = tau_emp(:);
    kappa_emp = kappa_emp(:);

    err = nan(size(tau_th));
    if strcmpi(err_metric, 'kappa')
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
    n = (0:n_plot-1).';

    if strcmpi(feat, 'real')
        x_cache = real(ld .^ n);
        if abs(ld) <= 1
            mx = 1;
        else
            mx = max(abs(x_cache));
            log_mx = (T - 1) * log(abs(ld));
            if isfinite(log_mx) && log_mx < log(realmax)
                mx = max(mx, exp(log_mx));
            end
        end
    else
        abs_ld = abs(ld);
        if abs_ld == 0
            x_cache = zeros(n_plot, 1);
            x_cache(1) = 1;
            mx = 1;
        elseif abs_ld <= 1
            x_cache = abs_ld .^ n;
            mx = 1;
        else
            log_mx = (T - 1) * log(abs_ld);
            log_vals = n * log(abs_ld);
            if isfinite(log_mx) && log_mx < log(realmax)
                mx = exp(log_mx);
                x_cache = exp(log_vals);
            else
                mx = 1;
                x_cache = exp(log_vals - log_mx);
            end
        end
    end

    if mx > 0
        x_plot = x_cache / mx;
    else
        x_plot = x_cache;
    end
end
