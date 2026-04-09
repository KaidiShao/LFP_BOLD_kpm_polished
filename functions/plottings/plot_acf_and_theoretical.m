function fig = plot_acf_and_theoretical(rep, EDMD_outputs, cfg)
% PLOT_ACF_AND_THEORETICAL
%   1×4 layout:
%     (1) ACFs of all Koopman eigenfunctions (normalized)
%     (2) Theoretical eigenfunctions: exp(lambda_c * t) for all eigenvalues
%     (3) ACFs of all reduced components (Z_time)
%     (4) ACFs of all smoothed reduced components (Z_time_smooth)
%
%   INPUT:
%     rep.Z_time         [T × K]          low-d trajectories
%     rep.Z_time_smooth  [T × K]          smoothed trajectories (optional)
%     rep.feature        'efuns_abs' / 'efuns_real' etc. (to choose Phi)
%     EDMD_outputs.evalues      [N_eig × 1] discrete eigenvalues
%     EDMD_outputs.efuns / norm_efuns.abs / norm_efuns.real
%     cfg.dt             scalar sampling interval (optional)
%
%   Normalisation convention:
%     - Eigenfunctions and theoretical curves are scaled only by max(abs)
%     - No mean subtraction is applied during visualisation normalisation

    % ----------- slow variables (low-d components) -----------
    Z  = rep.Z_time;
    Zs = rep.Z_time_smooth;
    if isempty(Z)
        error('plot_acf_and_theoretical: rep.Z_time is empty.');
    end
    if isempty(Zs)
        Zs = Z;
    end
    [Tz, Kslow] = size(Z);

    % ----------- sampling interval & time vector ------------
    if isfield(cfg, 'dt') && ~isempty(cfg.dt)
        dt = cfg.dt;
    else
        dt = 1;
    end

    % ----------- Koopman eigenfunctions Phi(t,i) ------------
    % choose real or abs normalisation based on feature name
    if contains(rep.feature, 'real')
        if isfield(EDMD_outputs, 'norm_efuns') && isfield(EDMD_outputs.norm_efuns, 'real')
            Phi = EDMD_outputs.norm_efuns.real;    % [T × N]
        else
            Phi = normalize_efun(EDMD_outputs.efuns, 'real');
        end
    else
        if isfield(EDMD_outputs, 'norm_efuns') && isfield(EDMD_outputs.norm_efuns, 'abs')
            Phi = EDMD_outputs.norm_efuns.abs;     % [T × N]
        else
            Phi = normalize_efun(EDMD_outputs.efuns, 'abs');
        end
    end

    [Tphi, Neig] = size(Phi);

    % align time length across Phi and Z
    T = min(Tphi, Tz);
    Phi = Phi(1:T, :);
    Z   = Z(1:T, :);
    Zs  = Zs(1:T, :);
    t   = (0:T-1) * dt;

    % ----------- eigenvalues -----------
    if ~isfield(EDMD_outputs, 'evalues') || isempty(EDMD_outputs.evalues)
        error('plot_acf_and_theoretical: EDMD_outputs.evalues is missing.');
    end
    lambda_d = EDMD_outputs.evalues(:);     % [N_eig × 1]

    % continuous eigenvalues via bilinear transform
    lambda_c = (2/dt) * (lambda_d - 1) ./ (lambda_d + 1);   % [N_eig × 1]

    max(abs(lambda_c));
    % ----------- max lag for ACF -----------
    maxLag = min(200, T-1);

    % create figure
    fig = figure('Name', ['ACF & theoretical: ' rep.rep_id], ...
                 'Color','w', ...
                 'Position',[2378, 559, 1021, 301]);

    % =======================================================
    % (1) ACFs of all Koopman eigenfunctions
    % =======================================================
    subplot(1,4,1); hold on;
    if exist('othercolor', 'file') == 2
        cmap_eig = othercolor('Spectral10');
    else
        cmap_eig = lines(max(Neig, 10));
    end
    % if fewer colors than eigenfunctions, tile them
    if size(cmap_eig,1) < Neig
        cmap_eig = interp1( ...
            linspace(0,1,size(cmap_eig,1)), cmap_eig, ...
            linspace(0,1,Neig));
    end

    for i = 1:Neig
        x = Phi(:,i);
        [lags_i, acf_i] = local_acf(x, maxLag);
        plot(lags_i * dt, acf_i, 'Color', cmap_eig(i,:), 'LineWidth', 0.8);
    end
    xlabel('lag (time)');
    ylabel('ACF');
    title('ACFs of all eigenfunctions');
    grid on;
    xlim([0,20]);
    % do not clutter with legend

    % =======================================================
    % (2) Theoretical eigenfunctions: exp(lambda_c * t)
    % =======================================================
    subplot(1,4,2); hold on;
    % use same colormap for consistency
    % discrete-time index n = 0,1,...,T-1
    n = (0:T-1).';           % [T×1]
    
    for i = 1:Neig
        ld = lambda_d(i);    % discrete eigenvalue λ_d

        % theoretical discrete-time eigenfunction: λ_d^n
        if contains(rep.feature, 'real')
            x_th = real(ld.^n);      % real part of λ_d^n
        else
            x_th = abs(ld.^n);       % magnitude of λ_d^n
        end

        % normalise to avoid crazy scales
        max_abs = max(abs(x_th));
        if max_abs > 0
            x_th = x_th / max_abs;
        end

        plot(t, x_th, 'Color', cmap_eig(i,:), 'LineWidth', 0.8); hold on
    end

    xlabel('time');
    ylabel('amplitude (normalised)');
    title('Theoretical eigenfunctions: \lambda_d^n');
    grid on;
    xlim([0,20])

    % =======================================================
    % (3) ACFs of reduced components (Z_time)
    % =======================================================
    subplot(1,4,3); hold on;
    cmap_red = lines(max(Kslow,1));

    for k = 1:Kslow
        z_k = Z(:,k);
        [lags_k, acf_k] = local_acf(z_k, maxLag);
        plot(lags_k * dt, acf_k, 'Color', cmap_red(k,:), 'LineWidth', 1.0);
    end
    xlabel('lag (time)');
    ylabel('ACF');
    title('ACFs of reduced components (Z\_time)');
    grid on;
    legend(arrayfun(@(k) sprintf('dim %d',k), 1:Kslow, 'UniformOutput', false), ...
           'Location','best');

    xlim([0,20])
    % =======================================================
    % (4) ACFs of smoothed reduced components (Z_time_smooth)
    % =======================================================
    subplot(1,4,4); hold on;

    for k = 1:Kslow
        z_ks = Zs(:,k);
        [lags_k, acf_k] = local_acf(z_ks, maxLag);
        plot(lags_k * dt, acf_k, 'Color', cmap_red(k,:), 'LineWidth', 1.0);
    end
    xlabel('lag (time)');
    ylabel('ACF');
    title('ACFs of smoothed reduced components');
    grid on;
    legend(arrayfun(@(k) sprintf('dim %d',k), 1:Kslow, 'UniformOutput', false), ...
           'Location','best');

end

% -------- helper: ACF with non-negative lags only --------
function [lagsPositive, acfPositive] = local_acf(x, maxLag)
    x = x(:) - mean(x);
    [c, lags] = xcorr(x, maxLag, 'coeff');
    idx = lags >= 0;
    lagsPositive = lags(idx);
    acfPositive  = c(idx);
end
