function residual_cfg = build_blp_residual_cfg_from_cross_params(params)
%BUILD_BLP_RESIDUAL_CFG_FROM_CROSS_PARAMS Build shared residual mode-selection config from cross-correlation params.

residual_cfg = struct();
residual_cfg.abs_thresh = params.post.abs_thresh;
residual_cfg.sort_by = params.post.sort_by;
residual_cfg.sort_dir = params.post.sort_dir;
residual_cfg.max_basis = params.post.max_basis;
residual_cfg.max_modes = params.residual.max_modes;
residual_cfg.lambdaType = params.residual.lambdaType;
residual_cfg.dt = params.residual.dt;
end
