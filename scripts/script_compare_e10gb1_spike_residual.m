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
params.spike_channels = 'all';

% Optional:
% params.session_filter_ids = [6];
% params.dictionary_file = 'E:\DataPons_processed\e10gb1\reskoopnet_dictionary\e10gb1_low50_high250_g2_abs_single.mat';

params.output_root = io_project.get_project_processed_root();
params.save_dir = fullfile(params.output_root, cfg.file_stem, 'spike_residual_comparison');
params.save_results = true;
params.verbose = true;
params.progress_every = 50;
params.top_n_rows = 100;

result = compute_spike_residual_comparison(cfg, params);

fprintf('\nSaved spike-residual comparison outputs:\n');
fprintf('  MAT  : %s\n', result.save_paths.main_mat);
fprintf('  CSV  : %s\n', result.save_paths.pooled_corr_csv);
fprintf('  Top  : %s\n', result.save_paths.pooled_top_csv);

if ~isempty(result.pooled_top_corr)
    disp(result.pooled_top_corr(1:min(10, height(result.pooled_top_corr)), :));
end
