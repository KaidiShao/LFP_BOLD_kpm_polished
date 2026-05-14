function mode_info = build_blp_residual_mode_info(evalues_in, residual_cfg)
%BUILD_BLP_RESIDUAL_MODE_INFO Select and order modes for residual analysis.

evalues0 = evalues_in(:);

switch lower(residual_cfg.sort_by)
    case {'modulus', 'abs'}
        key_full = abs(evalues0);
    case {'real', 'realpart'}
        key_full = real(evalues0);
    otherwise
        error('Unknown residual sort_by = %s.', residual_cfg.sort_by);
end

[~, ord_full] = sort(key_full, residual_cfg.sort_dir);
evalues_sorted = evalues0(ord_full);
mask_sorted = abs(evalues_sorted) > residual_cfg.abs_thresh;

if ~any(mask_sorted)
    error('Residual threshold removed all modes.');
end

idx_sorted_in_original = ord_full(mask_sorted);
Kkeep = min(numel(idx_sorted_in_original), residual_cfg.max_basis);
Ksel = min(Kkeep, residual_cfg.max_modes);

selected_idx = idx_sorted_in_original(1:Ksel);
selected_evalues = evalues0(selected_idx);
lambda_d = local_to_discrete_lambda(selected_evalues, residual_cfg);

mode_info = struct();
mode_info.selected_idx = selected_idx(:).';
mode_info.selected_evalues = selected_evalues(:);
mode_info.lambda_d = lambda_d(:);
mode_info.lambda_d_row = reshape(lambda_d(:).', 1, []);
mode_info.n_modes = numel(selected_idx);
mode_info.abs_thresh = residual_cfg.abs_thresh;
mode_info.sort_by = residual_cfg.sort_by;
mode_info.sort_dir = residual_cfg.sort_dir;
mode_info.max_basis = residual_cfg.max_basis;
mode_info.max_modes = residual_cfg.max_modes;
end


function lambda_d = local_to_discrete_lambda(lambda_in, residual_cfg)
switch lower(residual_cfg.lambdaType)
    case 'discrete'
        lambda_d = lambda_in;
    case 'continuous'
        lambda_d = exp(lambda_in .* residual_cfg.dt);
    otherwise
        error('Unknown residual lambdaType = %s.', residual_cfg.lambdaType);
end
end
