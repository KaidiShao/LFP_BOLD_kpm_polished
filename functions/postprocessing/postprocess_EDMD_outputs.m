% https://chatgpt.com/c/6969d5fc-33a4-8329-a9e5-7e8b9a0946b7

function [EDMD_outputs, fig] = postprocess_EDMD_outputs(EDMD_outputs, opts)
%POSTPROCESS_EDMD_OUTPUTS Threshold + sort + subset + optional plotting.
%
% Required fields in EDMD_outputs:
%   EDMD_outputs.evalues      [K x 1] (or [1 x K])
%   EDMD_outputs.efuns        [T x K]
%   EDMD_outputs.kpm_modes    [K x N]  (or compatible with eigen-basis dimension = K)
%
% opts fields (with defaults):
%   opts.abs_thresh   = 0.01;         % pre-mask: keep abs(evalues) > abs_thresh
%   opts.sort_by      = 'modulus';    % 'modulus' | 'real'
%   opts.sort_dir     = 'descend';    % 'ascend' | 'descend'
%   opts.dt           = [];           % sampling interval for bilinear transform (optional)
%   opts.do_plot      = false;        % whether to plot the 2x3 diagnostic figure
%   opts.t_plot       = [];           % explicit indices for the plot window
%   opts.window_idx   = [];           % alias of opts.t_plot for windowed plotting
%   opts.max_plot_samples = 2000;     % default window length when no window is supplied
%   opts.window_start = 1;            % default window start when no window is supplied
%   opts.time_vec     = [];           % x-axis vector for time; if empty uses t_plot index
%   opts.session_border = [];         % optional xline positions in current x-axis units
%   opts.draw_border  = false;        % draw session_border if provided
%   opts.colormap_fun = [];           % function handle for colormap (optional), e.g. @() flipud(othercolor('Spectral10'))
%
% Output (EDMD_outputs updated):
%   EDMD_outputs.original                 (backup of original if not present)
%   EDMD_outputs.original_sorted          (sorted full set)
%   EDMD_outputs.mask                     (logical mask after threshold in sorted coordinates)
%   EDMD_outputs.masked                   (masked subset, already sorted)
%   EDMD_outputs.norm                     (normalized masked efuns)
%   EDMD_outputs.evalues_bilinear         (masked, if dt provided)

% -------------------- defaults --------------------
if nargin < 2, opts = struct(); end
if ~isfield(opts, 'abs_thresh'),   opts.abs_thresh = 0.01; end
if ~isfield(opts, 'sort_by'),      opts.sort_by    = 'modulus'; end
if ~isfield(opts, 'sort_dir'),     opts.sort_dir   = 'descend'; end
if ~isfield(opts, 'max_basis'), opts.max_basis = Inf; end   % take first K after sorting (Inf = keep all)
if ~isfield(opts, 'dt'),           opts.dt         = []; end
if ~isfield(opts, 'do_plot'),      opts.do_plot    = false; end
if ~isfield(opts, 't_plot'),       opts.t_plot     = []; end
if ~isfield(opts, 'window_idx'),   opts.window_idx = []; end
if ~isfield(opts, 'max_plot_samples'), opts.max_plot_samples = 2000; end
if ~isfield(opts, 'window_start'), opts.window_start = 1; end
if ~isfield(opts, 'time_vec'),     opts.time_vec   = []; end
if ~isfield(opts, 'session_border'), opts.session_border = []; end
if ~isfield(opts, 'draw_border'),  opts.draw_border = false; end
if ~isfield(opts, 'colormap_fun'), opts.colormap_fun = []; end

fig = [];

% -------------------- backup originals (once) --------------------
if ~isfield(EDMD_outputs, 'original')
    EDMD_outputs.original.evalues   = EDMD_outputs.evalues;
    EDMD_outputs.original.efuns     = EDMD_outputs.efuns;
    EDMD_outputs.original.kpm_modes = EDMD_outputs.kpm_modes;
end

% -------------------- normalization on MASKED efuns --------------------


% Ensure column vector eigenvalues
evalues0 = EDMD_outputs.original.evalues(:);
efuns0   = EDMD_outputs.original.efuns;
modes0   = EDMD_outputs.original.kpm_modes;

K = numel(evalues0);
if size(efuns0, 2) ~= K
    error('Dimension mismatch: size(efuns,2)=%d but numel(evalues)=%d.', size(efuns0,2), K);
end
if size(modes0, 1) ~= K
    error('Dimension mismatch: size(kpm_modes,1)=%d but numel(evalues)=%d.', size(modes0,1), K);
end


% -------------------- threshold pre-mask in ORIGINAL coordinates --------------------
mask_thresh0 = abs(evalues0) > opts.abs_thresh;
if ~any(mask_thresh0)
    error('Threshold removed all eigen-basis: abs_thresh=%g. Nothing left after masking.', opts.abs_thresh);
end
% -------------------- sort the FULL set (all K) --------------------
switch lower(opts.sort_by)
    case {'modulus','abs'}
        key_full = abs(evalues0);
    case {'real','realpart'}
        key_full = real(evalues0);
    otherwise
        error('Unknown opts.sort_by = %s. Use ''modulus'' or ''real''.', opts.sort_by);
end

[~, ord_full] = sort(key_full, opts.sort_dir);

evalues_sorted_full = evalues0(ord_full);
efuns_sorted_full   = efuns0(:, ord_full);
modes_sorted_full   = modes0(ord_full, :);
%%

%%
% Save "original_sorted" as the FULL sorted set (all dimensions)
EDMD_outputs.original_sorted.evalues    = evalues_sorted_full;
EDMD_outputs.original_sorted.efuns      = efuns_sorted_full;
EDMD_outputs.original_sorted.kpm_modes  = modes_sorted_full;
EDMD_outputs.original_sorted.sort_by    = opts.sort_by;
EDMD_outputs.original_sorted.sort_dir   = opts.sort_dir;
EDMD_outputs.original_sorted.abs_thresh = opts.abs_thresh;
EDMD_outputs.original_sorted.max_basis  = opts.max_basis;

% -------------------- threshold mask on the FULL sorted set --------------------
mask_sorted = abs(evalues_sorted_full) > opts.abs_thresh;
if ~any(mask_sorted)
    error('Threshold removed all eigen-basis: abs_thresh=%g. Nothing left after masking.', opts.abs_thresh);
end

% Store the threshold mask in ORIGINAL coordinates (optional but useful)
mask_thresh0 = false(size(evalues0));
mask_thresh0(ord_full(mask_sorted)) = true;
EDMD_outputs.mask = mask_thresh0;


% -------------------- masked set: thresholded + sorted (then take first K) --------------------
evalues_masked = evalues_sorted_full(mask_sorted);
efuns_masked   = efuns_sorted_full(:, mask_sorted);
modes_masked   = modes_sorted_full(mask_sorted, :);



EDMD_outputs.norm_efuns.real = local_normalize_efun(efuns_masked, 'real');
EDMD_outputs.norm_efuns.abs  = local_normalize_efun(efuns_masked, 'abs');


% Apply max_basis ONLY to masked
if isfinite(opts.max_basis)
    Kkeep = min(numel(evalues_masked), opts.max_basis);
else
    Kkeep = numel(evalues_masked);
end

evalues_masked = evalues_masked(1:Kkeep);
efuns_masked   = efuns_masked(:, 1:Kkeep);
modes_masked   = modes_masked(1:Kkeep, :);

EDMD_outputs.evalues   = evalues_masked;
EDMD_outputs.efuns     = efuns_masked;
EDMD_outputs.kpm_modes = modes_masked;

% Indices of the final selected basis in ORIGINAL coordinates
idx_sorted_in_original = ord_full(mask_sorted);      % original indices after sorting+threshold
EDMD_outputs.idx_final_in_original = idx_sorted_in_original(1:Kkeep);

% Bookkeeping for masked
EDMD_outputs.masked.evalues    = EDMD_outputs.evalues;
EDMD_outputs.masked.efuns      = EDMD_outputs.efuns;
EDMD_outputs.masked.kpm_modes  = EDMD_outputs.kpm_modes;
EDMD_outputs.masked.abs_thresh = opts.abs_thresh;
EDMD_outputs.masked.sort_by    = opts.sort_by;
EDMD_outputs.masked.sort_dir   = opts.sort_dir;
EDMD_outputs.masked.max_basis  = opts.max_basis;




% -------------------- bilinear transform (optional) --------------------
if ~isempty(opts.dt)
    dt = opts.dt;
    EDMD_outputs.evalues_bilinear = (2/dt) * (EDMD_outputs.evalues - 1) ./ (EDMD_outputs.evalues + 1);
end

% -------------------- plotting (optional) --------------------
if opts.do_plot
    T = size(efuns0, 1);
    t_plot = local_resolve_plot_window(T, opts);

    if isempty(opts.time_vec)
        x = t_plot;
        xlab = 'Time (index)';
    else
        x = opts.time_vec(t_plot);
        xlab = 'Time';
    end

    % Choose colormap
    if isempty(opts.colormap_fun)
        use_cmap = [];
    else
        use_cmap = opts.colormap_fun;
    end

    session_border = opts.session_border;
    % Prepare matrices for heatmaps
    A_abs_orig = EDMD_outputs.norm_efuns.abs(t_plot, :).';   % [K_all x |t_plot|]
    A_abs_mask = EDMD_outputs.norm_efuns.abs(t_plot, 1:Kkeep).';                   % [K_mask x |t_plot|]
    A_re_orig  = EDMD_outputs.norm_efuns.real(t_plot, :).';
    A_re_mask  = EDMD_outputs.norm_efuns.real(t_plot, 1:Kkeep).';

    % Layout: use subplot(2,7,...) to make the first column narrow
    fig = figure('Color','w');

    % (1,1) original eigenvalues (complex plane) - narrow
    subplot(2,7,1);
    plot(real(EDMD_outputs.original_sorted.evalues), imag(EDMD_outputs.original_sorted.evalues), '.', 'MarkerSize', 10);
    hold on;
    theta = linspace(0, 2*pi, 400);
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);  % dashed unit circle
    hold off;
    grid on; axis equal; axis([-1.1,1.1,-1.1,1.1])
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title(sprintf('Original eigenvalues (sorted by %s)', opts.sort_by), 'Interpreter','none');

    % (2,1) masked eigenvalues (complex plane) - narrow
    subplot(2,7,8);
    plot(real(EDMD_outputs.evalues), imag(EDMD_outputs.evalues), '.', 'MarkerSize', 10);
    hold on;
    theta = linspace(0, 2*pi, 400);
    plot(cos(theta), sin(theta), 'k--', 'LineWidth', 1);  % dashed unit circle
    hold off;
    grid on; axis equal; axis([-1.1,1.1,-1.1,1.1])
    xlabel('Re(\lambda)'); ylabel('Im(\lambda)');
    title(sprintf('Masked eigenvalues (|\\lambda|>%g)', opts.abs_thresh));

    % (1,2) original efuns abs heatmap - wide (columns 2:4)
    subplot(2,7,2:4);
    imagesc(x, 1:size(A_abs_orig,1), A_abs_orig);
    if opts.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir','reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title('|eigenfunctions| (original, sorted)');
    colorbar;
    if ~isempty(use_cmap), colormap(use_cmap()); end

    % (2,2) masked efuns abs heatmap - wide (columns 2:4)
    subplot(2,7,9:11);
    imagesc(x, 1:size(A_abs_mask,1), A_abs_mask);
    if opts.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir','reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title('|eigenfunctions| (masked)');
    colorbar;
    if ~isempty(use_cmap), colormap(use_cmap()); end

    % (1,3) original efuns real heatmap - wide (columns 5:7)
    subplot(2,7,5:7);
    imagesc(x, 1:size(A_re_orig,1), A_re_orig);
    if opts.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir','reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title('Re(eigenfunctions) (original, sorted)');
    colorbar;
    if ~isempty(use_cmap), colormap(use_cmap()); end

    % (2,3) masked efuns real heatmap - wide (columns 5:7)
    subplot(2,7,12:14);
    imagesc(x, 1:size(A_re_mask,1), A_re_mask);
    if opts.draw_border && ~isempty(session_border), xline(session_border, 'k--'); end
    set(gca, 'YDir','reverse');
    xlabel(xlab); ylabel('Eigenfunction #');
    title('Re(eigenfunctions) (masked)');
    colorbar;
    if ~isempty(use_cmap), colormap(use_cmap()); end
    set(fig, 'Position', [1440         654        1531         584]);
end
end


function efuns_norm = local_normalize_efun(efuns, mode)
% Normalize each column only by its maximum magnitude, without de-meaning.

if strcmpi(mode, 'real')
    X = real(efuns);
else
    X = abs(efuns);
end

efuns_norm = zeros(size(X), 'like', X);

for n_efun = 1:size(X, 2)
    x = X(:, n_efun);
    denom = max(abs(x));
    if denom > 0
        efuns_norm(:, n_efun) = x ./ denom;
    else
        efuns_norm(:, n_efun) = x;
    end
end
end


function t_plot = local_resolve_plot_window(T, opts)
% Resolve a single visualization window while keeping backward compatibility.

if ~isempty(opts.t_plot)
    idx = opts.t_plot;
elseif ~isempty(opts.window_idx)
    idx = opts.window_idx;
else
    win_len = min(T, max(1, round(opts.max_plot_samples)));
    start_idx = min(max(1, round(opts.window_start)), T - win_len + 1);
    idx = start_idx:(start_idx + win_len - 1);
end

t_plot = unique(round(idx(:).'));
t_plot = t_plot(t_plot >= 1 & t_plot <= T);

if isempty(t_plot)
    t_plot = 1:T;
end
end
