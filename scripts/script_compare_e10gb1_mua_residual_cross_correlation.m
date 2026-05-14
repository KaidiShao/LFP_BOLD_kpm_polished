this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg = cfg_E10gb1();

params = struct();

% Leave empty to auto-pick the latest non-smoke EDMD output directory that
% matches params.source.name_contains under E:\autodl_results\e10gb1\mlp\outputs.
params.source_cfg = struct();
params.source_cfg.mode = 'chunk_dir';
params.source_cfg.data_dir = '';

params.source = struct();
params.source.base_dir = fullfile('E:\autodl_results', cfg.file_stem, 'mlp', 'outputs');
params.source.name_contains = 'projected_kv_abs';
params.source.prefer_non_smoke = true;

params.post = struct();
params.post.abs_thresh = 0.01;
params.post.sort_by = 'modulus';
params.post.sort_dir = 'descend';
params.post.max_basis = 30;

params.residual = struct();
params.residual.method = 'koopman_residual';
params.residual.lambda_source = 'edmd';
params.residual.lambdaType = 'discrete';
params.residual.first_u_mode = 'phi1';
params.residual.max_modes = 20;

params.blp_channels = 'selected';
params.blp_band = 'last';
params.pairings = {'abs_abs', 'raw_real', 'raw_imag'};

% Optional:
% params.session_filter_ids = [6];
% params.dictionary_file = 'E:\DataPons_processed\e10gb1\pipeline1_reskoopnet_dictionary\e10gb1_low50_high250_g2_abs_single.mat';
% params.blp_channels = 'all';
% params.blp_band = 3;

params.output_root = io_project.get_project_processed_root();
params.save_dir = io_project.get_pipeline_stage_dir( ...
    params.output_root, cfg, 6, 'mua_residual_cross_correlation');
params.verbose = true;
params.progress_every = 50;
params.top_n_rows = 100;

result = compute_mua_residual_cross_correlation(cfg, params);
result = save_mua_residual_cross_correlation_result(result, params.save_dir);

fprintf('\nSaved MUA-residual cross-correlation outputs:\n');
fprintf('  MAT  : %s\n', result.save_paths.main_mat);
fprintf('  CSV  : %s\n', result.save_paths.pooled_corr_csv);
fprintf('  Top  : %s\n', result.save_paths.pooled_top_csv);
fprintf('  Band : %s (index %d of %d)\n', ...
    result.source.blp_band_selection_label, ...
    result.source.blp_band_index, ...
    result.source.n_blp_bands);

if ~isempty(result.pooled_top_corr)
    disp(result.pooled_top_corr(1:min(10, height(result.pooled_top_corr)), :));
end
