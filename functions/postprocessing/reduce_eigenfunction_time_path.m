function reduction = reduce_eigenfunction_time_path(F_time_by_mode, cfg)
%REDUCE_EIGENFUNCTION_TIME_PATH Reduce eigenfunction features directly in time.

cfg = local_apply_defaults(cfg);

[T, N] = size(F_time_by_mode);
if T < 2 || N < 1
    error('F_time_by_mode must be T x N with T >= 2 and N >= 1.');
end

switch lower(cfg.method)
    case 'svd'
        X = F_time_by_mode;
        baseline = zeros(1, N);
        if cfg.center_modes
            baseline = mean(X, 1);
            X = X - baseline;
        end

        [U, S, V] = svd(X, 'econ');
        k = min([cfg.n_components, size(U, 2), size(V, 2)]);

        C = U(:, 1:k) * S(1:k, 1:k);
        W = V(:, 1:k);
        Xhat = C * W.';
        recon = Xhat + baseline;

        singular_values = diag(S);
        energy = singular_values .^ 2;
        explained_ratio = energy(1:k) ./ max(sum(energy), eps);

        aux = struct();
        aux.time = struct();
        aux.time.method = 'SVD';
        aux.time.center_modes = cfg.center_modes;
        aux.time.baseline_mode = baseline(:);
        aux.time.singular_values = singular_values;

        quality = struct();
        quality.explained_variance_ratio = explained_ratio(:);
        quality.reconstruction_error_fro = norm(F_time_by_mode - recon, 'fro');
        quality.reconstruction_error_domain = 'feature';

        domain_name = 'feature';

    case 'logsvd'
        X = log(abs(F_time_by_mode) + cfg.log_epsilon);
        baseline = zeros(1, N);
        if cfg.center_modes
            baseline = mean(X, 1);
            X = X - baseline;
        end

        [U, S, V] = svd(X, 'econ');
        k = min([cfg.n_components, size(U, 2), size(V, 2)]);

        C = U(:, 1:k) * S(1:k, 1:k);
        W = V(:, 1:k);
        Xhat = C * W.';
        recon = Xhat + baseline;

        singular_values = diag(S);
        energy = singular_values .^ 2;
        explained_ratio = energy(1:k) ./ max(sum(energy), eps);

        aux = struct();
        aux.time = struct();
        aux.time.method = 'logSVD';
        aux.time.center_modes = cfg.center_modes;
        aux.time.baseline_mode = baseline(:);
        aux.time.log_epsilon = cfg.log_epsilon;
        aux.time.singular_values = singular_values;

        quality = struct();
        quality.explained_variance_ratio = explained_ratio(:);
        quality.reconstruction_error_fro = norm(log(abs(F_time_by_mode) + cfg.log_epsilon) - recon, 'fro');
        quality.reconstruction_error_domain = 'logabs';

        domain_name = 'logabs';

    case 'nmf'
        [X_model, offset] = local_make_nonnegative(F_time_by_mode, cfg.nmf_nonnegative_strategy);
        opts = statset('MaxIter', cfg.nmf_max_iter);

        k = min([cfg.n_components, size(X_model, 1), size(X_model, 2)]);
        [C, H] = nnmf(X_model, k, 'algorithm', 'als', 'options', opts);
        W = H.';
        Xhat_model = C * W.';
        recon = Xhat_model - offset;

        aux = struct();
        aux.time = struct();
        aux.time.method = 'NMF';
        aux.time.nmf_nonnegative_strategy = cfg.nmf_nonnegative_strategy;
        aux.time.input_offset = offset;

        quality = struct();
        quality.explained_variance_ratio = [];
        quality.reconstruction_error_fro = norm(F_time_by_mode - recon, 'fro');
        quality.reconstruction_error_domain = 'feature';

        domain_name = 'feature';

    case 'nmf2'
        if exist('meta_nmf', 'file') ~= 2
            error(['meta_nmf was not found on the MATLAB path. ', ...
                'Add it before using method ''NMF2''.']);
        end

        [X_model, offset] = local_make_nonnegative(F_time_by_mode, cfg.nmf_nonnegative_strategy);
        k = min([cfg.n_components, size(X_model, 1), size(X_model, 2)]);

        nmfpars = struct();
        nmfpars.init = 'random';
        res = meta_nmf(X_model.', k, 1, 'ismufast', nmfpars);

        W = res.w;
        C = res.h.';
        Xhat_model = C * W.';
        recon = Xhat_model - offset;

        aux = struct();
        aux.time = struct();
        aux.time.method = 'NMF2';
        aux.time.nmf_nonnegative_strategy = cfg.nmf_nonnegative_strategy;
        aux.time.input_offset = offset;

        quality = struct();
        quality.explained_variance_ratio = [];
        quality.reconstruction_error_fro = norm(F_time_by_mode - recon, 'fro');
        quality.reconstruction_error_domain = 'feature';

        domain_name = 'feature';

    otherwise
        error('Unknown cfg.method = %s.', cfg.method);
end

core = struct();
core.n_components = size(C, 2);
core.temporal_components_time_by_comp = C;
core.mode_weights_mode_by_comp = W;
core.mode_weights_norm_mode_by_comp = local_normalize_columns_absmax(W);
core.reconstruction_time_by_mode = recon;
core.component_domain = domain_name;
core.reconstruction_domain = quality.reconstruction_error_domain;
core.weight_semantics = 'loading';

reduction = struct();
reduction.core = core;
reduction.quality = quality;
reduction.aux = aux;
end


function cfg = local_apply_defaults(cfg)
if ~isfield(cfg, 'method') || isempty(cfg.method)
    cfg.method = 'SVD';
end

if ~isfield(cfg, 'n_components') || isempty(cfg.n_components)
    cfg.n_components = 4;
end

if ~isfield(cfg, 'center_modes') || isempty(cfg.center_modes)
    cfg.center_modes = true;
end

if ~isfield(cfg, 'log_epsilon') || isempty(cfg.log_epsilon)
    cfg.log_epsilon = 1e-10;
end

if ~isfield(cfg, 'nmf_max_iter') || isempty(cfg.nmf_max_iter)
    cfg.nmf_max_iter = 1000;
end

if ~isfield(cfg, 'nmf_nonnegative_strategy') || isempty(cfg.nmf_nonnegative_strategy)
    cfg.nmf_nonnegative_strategy = 'shift_global';
end
end


function [X_model, offset] = local_make_nonnegative(X, strategy)
offset = 0;

if all(X(:) >= 0)
    X_model = X;
    return;
end

switch lower(strategy)
    case 'shift_global'
        offset = -min(X(:));
        X_model = X + offset;

    case 'clip_zero'
        X_model = max(X, 0);

    otherwise
        error('Unknown nmf_nonnegative_strategy = %s.', strategy);
end
end


function W_norm = local_normalize_columns_absmax(W)
W_norm = W;

for i = 1:size(W, 2)
    denom = max(abs(W(:, i)));
    if denom > 0
        W_norm(:, i) = W(:, i) ./ denom;
    end
end
end
