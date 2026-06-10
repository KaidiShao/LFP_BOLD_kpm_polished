function cmp_params = build_mua_residual_cross_correlation_params(run_info, residual_bundle, save_dir, params)
%BUILD_MUA_RESIDUAL_CROSS_CORRELATION_PARAMS Build MUA xcorr params for one run.

cmp_params = struct();
cmp_params.source_cfg = struct();
cmp_params.source_cfg.mode = 'residual_workspace';
cmp_params.source_cfg.residual_bundle = residual_bundle;
cmp_params.source_cfg.data_dir = run_info.output_dir;
cmp_params.source = struct();
cmp_params.source.base_dir = fileparts(run_info.output_dir);
cmp_params.source.name_contains = run_info.run_name;
cmp_params.source.prefer_non_smoke = true;
cmp_params.filename_pattern = params.filename_pattern;
cmp_params.variable_name = params.variable_name;
cmp_params.post = struct();
cmp_params.post.abs_thresh = params.abs_thresh;
cmp_params.post.sort_by = params.sort_by;
cmp_params.post.sort_dir = params.sort_dir;
cmp_params.post.max_basis = params.max_basis;
cmp_params.residual = struct();
cmp_params.residual.max_modes = params.mua_cross_max_modes;
cmp_params.blp_channels = params.mua_cross_channels;
cmp_params.blp_band = params.mua_cross_band;
cmp_params.pairings = params.mua_cross_pairings;
cmp_params.output_root = params.processed_root;
cmp_params.save_dir = save_dir;
cmp_params.verbose = params.mua_cross_verbose;
cmp_params.progress_every = params.mua_cross_progress_every;
cmp_params.top_n_rows = params.mua_cross_top_n_rows;
end
