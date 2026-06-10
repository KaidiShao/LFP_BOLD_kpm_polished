function u_full = compute_blp_residual_matrix(phi_full, lambda_d_row, first_u_mode)
%COMPUTE_BLP_RESIDUAL_MATRIX Compute residuals for one full run workspace matrix.

[T, K] = size(phi_full);
u_full = zeros(T, K, 'like', phi_full);

if T < 1
    return;
end

switch lower(first_u_mode)
    case 'phi1'
        u_full(1, :) = phi_full(1, :);
    case 'zero'
        u_full(1, :) = 0;
    otherwise
        error('Unknown first_u_mode = %s.', first_u_mode);
end

if T >= 2
    u_full(2:end, :) = phi_full(2:end, :) - phi_full(1:end-1, :) .* lambda_d_row;
end
end
