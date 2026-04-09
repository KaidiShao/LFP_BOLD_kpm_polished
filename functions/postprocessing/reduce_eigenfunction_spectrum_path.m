function reduction = reduce_eigenfunction_spectrum_path(F_time_by_mode, cfg)
%REDUCE_EIGENFUNCTION_SPECTRUM_PATH Reduce eigenfunction features through mode geometry.

cfg = local_apply_defaults(cfg);

[T, N] = size(F_time_by_mode);
if T < 2 || N < 2
    error('F_time_by_mode must be T x N with T >= 2 and N >= 2.');
end

D = local_compute_mode_distance(F_time_by_mode, cfg.distance);
Z = local_embed_modes(F_time_by_mode, D, cfg);
[A, labels, embedding_info] = local_cluster_embedding(Z, cfg);

sign_alignment = ones(N, 1);
F_aligned = F_time_by_mode;
if strcmpi(cfg.feature_variant, 'real') && ~strcmpi(cfg.sign_alignment, 'none')
    sign_alignment = local_align_mode_signs(F_time_by_mode, A);
    F_aligned = F_time_by_mode .* reshape(sign_alignment, 1, []);
end

P = local_normalize_assignment_columns(A);
C = F_aligned * P;
W_recon = local_solve_reconstruction_weights(C, F_aligned, cfg.reconstruction_weight_method);
Fhat = C * W_recon.';

core = struct();
core.n_components = size(C, 2);
core.temporal_components_time_by_comp = C;
core.mode_weights_mode_by_comp = W_recon;
core.mode_weights_norm_mode_by_comp = local_normalize_columns_absmax(W_recon);
core.reconstruction_time_by_mode = Fhat;
core.component_domain = 'feature';
core.reconstruction_domain = 'feature';
core.weight_semantics = 'reconstruction_weight';

quality = struct();
quality.explained_variance_ratio = [];
quality.reconstruction_error_fro = norm(F_aligned - Fhat, 'fro');
quality.reconstruction_error_domain = 'feature';
quality.assignment_reconstruction_error_fro = norm(F_aligned - C * P.', 'fro');

aux = struct();
aux.spectrum = struct();
aux.spectrum.distance_name = cfg.distance;
aux.spectrum.distance_mode_by_mode = D;
aux.spectrum.embedding_method = cfg.method;
aux.spectrum.embedding_mode_by_dim = Z;
aux.spectrum.embedding_info = embedding_info;
aux.spectrum.cluster_method = cfg.cluster_method;
aux.spectrum.cluster_labels_mode = labels;
aux.spectrum.assignment_weights_mode_by_comp = A;
aux.spectrum.assignment_weights_norm_mode_by_comp = P;
aux.spectrum.sign_alignment_mode = sign_alignment;

reduction = struct();
reduction.core = core;
reduction.quality = quality;
reduction.aux = aux;
end


function cfg = local_apply_defaults(cfg)
if ~isfield(cfg, 'method') || isempty(cfg.method)
    cfg.method = 'MDS';
end

if ~isfield(cfg, 'output_dim') || isempty(cfg.output_dim)
    cfg.output_dim = 3;
end

if ~isfield(cfg, 'n_components') || isempty(cfg.n_components)
    cfg.n_components = 4;
end

if ~isfield(cfg, 'distance') || isempty(cfg.distance)
    cfg.distance = 'corr_abs';
end

if ~isfield(cfg, 'downsample_step') || isempty(cfg.downsample_step)
    cfg.downsample_step = 10;
end

if ~isfield(cfg, 'cluster_method') || isempty(cfg.cluster_method)
    cfg.cluster_method = 'kmeans';
end

if ~isfield(cfg, 'cluster_replicates') || isempty(cfg.cluster_replicates)
    cfg.cluster_replicates = 10;
end

if ~isfield(cfg, 'gmm_replicates') || isempty(cfg.gmm_replicates)
    cfg.gmm_replicates = 3;
end

if ~isfield(cfg, 'gmm_regularization') || isempty(cfg.gmm_regularization)
    cfg.gmm_regularization = 1e-5;
end

if ~isfield(cfg, 'mds_max_iter') || isempty(cfg.mds_max_iter)
    cfg.mds_max_iter = 1000;
end

if ~isfield(cfg, 'tsne_perplexity') || isempty(cfg.tsne_perplexity)
    cfg.tsne_perplexity = 30;
end

if ~isfield(cfg, 'umap_dir')
    cfg.umap_dir = '';
end

if ~isfield(cfg, 'feature_variant') || isempty(cfg.feature_variant)
    cfg.feature_variant = 'abs';
end

if ~isfield(cfg, 'sign_alignment') || isempty(cfg.sign_alignment)
    cfg.sign_alignment = 'dominant_cluster';
end

if ~isfield(cfg, 'reconstruction_weight_method') || isempty(cfg.reconstruction_weight_method)
    cfg.reconstruction_weight_method = 'least_squares';
end
end


function D = local_compute_mode_distance(F, distance_name)
switch lower(distance_name)
    case 'corr_abs'
        C = corr(F);
        C(~isfinite(C)) = 0;
        D = 1 - abs(C);
        D(1:size(D, 1)+1:end) = 0;

    otherwise
        error('Unsupported distance_name = %s.', distance_name);
end
end


function Z = local_embed_modes(F_time_by_mode, D, cfg)
output_dim = min(cfg.output_dim, size(F_time_by_mode, 2) - 1);
output_dim = max(output_dim, 1);

switch lower(cfg.method)
    case 'mds'
        opts = statset('MaxIter', cfg.mds_max_iter);
        Z = mdscale(D, output_dim, ...
            'Criterion', 'stress', ...
            'Start', 'random', ...
            'Options', opts);

    case 'diffusion_map'
        epsilon = median(D(:));
        if ~isfinite(epsilon) || epsilon <= 0
            epsilon = 1;
        end
        K = exp(-(D .^ 2) / (epsilon ^ 2));
        P = K ./ max(sum(K, 2), eps);
        [V, L] = eig(P);
        [~, idx] = sort(diag(L), 'descend');
        V = V(:, idx);
        max_dim = min(output_dim + 1, size(V, 2));
        Z = V(:, 2:max_dim);

    case 'umap'
        if ~isempty(cfg.umap_dir)
            addpath(cfg.umap_dir);
        end
        if exist('run_umap', 'file') ~= 2
            error(['run_umap was not found on the MATLAB path. ', ...
                'Set cfg.path.spectrum.umap_dir before using UMAP.']);
        end
        step = max(1, cfg.downsample_step);
        X_mode_by_time = F_time_by_mode(1:step:end, :).';
        Z = run_umap(X_mode_by_time, ...
            'metric', 'correlation', ...
            'n_components', output_dim);

    case 'tsne'
        step = max(1, cfg.downsample_step);
        X_mode_by_time = F_time_by_mode(1:step:end, :).';
        Z = tsne(X_mode_by_time, ...
            'NumDimensions', output_dim, ...
            'Distance', 'correlation', ...
            'Perplexity', cfg.tsne_perplexity);

    otherwise
        error('Unknown cfg.method = %s.', cfg.method);
end
end


function [A, labels, info] = local_cluster_embedding(Z, cfg)
k = min(cfg.n_components, size(Z, 1));
info = struct();

switch lower(cfg.cluster_method)
    case 'kmeans'
        labels = kmeans(Z, k, ...
            'Replicates', cfg.cluster_replicates, ...
            'Distance', 'sqeuclidean');
        A = zeros(size(Z, 1), k);
        for i = 1:size(Z, 1)
            A(i, labels(i)) = 1;
        end
        info.gm = [];

    case 'gmm'
        gm = fitgmdist(Z, k, ...
            'Start', 'plus', ...
            'Replicates', cfg.gmm_replicates, ...
            'RegularizationValue', cfg.gmm_regularization, ...
            'Options', statset('MaxIter', 1000));
        labels = cluster(gm, Z);
        A = posterior(gm, Z);
        info.gm = gm;

    otherwise
        error('Unknown cfg.cluster_method = %s.', cfg.cluster_method);
end
end


function sign_alignment = local_align_mode_signs(F_time_by_mode, A)
N = size(F_time_by_mode, 2);
sign_alignment = ones(N, 1);
[~, dominant_label] = max(A, [], 2);

for k = 1:size(A, 2)
    idx = find(dominant_label == k);
    if isempty(idx)
        continue;
    end

    [~, best_local] = max(A(idx, k));
    ref_idx = idx(best_local);
    ref_trace = F_time_by_mode(:, ref_idx);

    if norm(ref_trace) <= eps
        continue;
    end

    for j = 1:numel(idx)
        mode_idx = idx(j);
        score = dot(F_time_by_mode(:, mode_idx), ref_trace);
        if score < 0
            sign_alignment(mode_idx) = -1;
        end
    end
end
end


function P = local_normalize_assignment_columns(A)
P = A;
for i = 1:size(A, 2)
    denom = sum(A(:, i));
    if denom > 0
        P(:, i) = A(:, i) ./ denom;
    else
        P(:, i) = 0;
    end
end
end


function W = local_solve_reconstruction_weights(C, F, method)
K = size(C, 2);
N = size(F, 2);

switch lower(method)
    case 'least_squares'
        W = (C \ F).';

    case 'nonnegative'
        W = zeros(N, K);
        for i = 1:N
            W(i, :) = lsqnonneg(C, F(:, i)).';
        end

    otherwise
        error('Unknown reconstruction_weight_method = %s.', method);
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
