function u_chunk = compute_blp_chunk_residual(phi_chunk, lambda_d_row, prev_phi_last, first_u_mode, is_first_chunk)
%COMPUTE_BLP_CHUNK_RESIDUAL Compute one streamed residual chunk.

[T, K] = size(phi_chunk);
u_chunk = zeros(T, K, 'like', phi_chunk);

if T < 1
    return;
end

if is_first_chunk
    switch lower(first_u_mode)
        case 'phi1'
            u_chunk(1, :) = phi_chunk(1, :);
        case 'zero'
            u_chunk(1, :) = 0;
        otherwise
            error('Unknown first_u_mode = %s.', first_u_mode);
    end
else
    u_chunk(1, :) = phi_chunk(1, :) - prev_phi_last .* lambda_d_row;
end

if T >= 2
    u_chunk(2:end, :) = phi_chunk(2:end, :) - phi_chunk(1:end-1, :) .* lambda_d_row;
end
end
