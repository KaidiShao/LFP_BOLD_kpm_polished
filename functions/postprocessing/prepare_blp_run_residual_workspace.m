function residual_bundle = prepare_blp_run_residual_workspace(EDMD_full, params)
%PREPARE_BLP_RUN_RESIDUAL_WORKSPACE Compute one reusable residual bundle for a run.

if nargin < 1 || isempty(EDMD_full)
    error('EDMD_full must be provided.');
end
if nargin < 2
    params = struct();
end

residual_cfg = local_build_residual_cfg(params);
mode_info = build_blp_residual_mode_info(EDMD_full.evalues, residual_cfg);

phi_full = EDMD_full.efuns(:, mode_info.selected_idx);
u_full = compute_blp_residual_matrix(phi_full, mode_info.lambda_d_row, residual_cfg.first_u_mode);

residual_bundle = struct();
residual_bundle.mode_info = mode_info;
residual_bundle.phi_full = phi_full;
residual_bundle.u_full = u_full;
residual_bundle.n_samples = size(u_full, 1);
residual_bundle.max_modes = mode_info.n_modes;
end


function residual_cfg = local_build_residual_cfg(params)
residual_cfg = struct();
residual_cfg.abs_thresh = params.abs_thresh;
residual_cfg.sort_by = params.sort_by;
residual_cfg.sort_dir = params.sort_dir;
residual_cfg.max_basis = params.max_basis;
residual_cfg.max_modes = max([1, params.spkt_cross_max_modes, params.mua_cross_max_modes]);
residual_cfg.lambdaType = 'discrete';
residual_cfg.first_u_mode = 'phi1';
end
