this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

log_root = fullfile(repo_root, 'tmp', 'pipeline6_logs');
if exist(log_root, 'dir') ~= 7
    mkdir(log_root);
end

log_file = fullfile(log_root, 'e10gb1_mua_residual_cross_correlation_run_log.txt');
fid = fopen(log_file, 'w');
if fid < 0
    error('Could not open log file for writing: %s', log_file);
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, 'START %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')));
fprintf(fid, 'PWD %s\n', pwd);
fprintf(fid, 'REPO %s\n', repo_root);

try
    cfg = struct();
    cfg.dataset_id = 'E10.gb1';
    cfg.raw_data_root = 'E:\DataPons\E10.gb1\';
    cfg.data_subfolder = 'blp';
    cfg.file_stem = 'e10gb1';
    cfg.channels = struct();
    cfg.channels.sites = {'lgn' 'lgn' 'st' 'lgn' 'lgn' 'lgn' 'st' 'st' 'hp' 'hp' 'hp' 'hp' 'pl' 'pl' 'pl' 'pl'};
    cfg.channels.selected_labels = {'hp', 'pl'};
    cfg.channels.selected_site1 = 9:12;
    cfg.channels.selected_site2 = 13:16;
    cfg.channels.selected_all = 9:16;
    cfg.sessions = struct([]);
    cfg.sessions(1).session_id = 1:5;
    cfg.sessions(1).include = true;
    cfg.sessions(1).selected_channels = cfg.channels.selected_all;
    cfg.sessions(1).notes = 'polar';
    cfg.sessions(2).session_id = [6:7 9:13 20:25];
    cfg.sessions(2).include = true;
    cfg.sessions(2).selected_channels = cfg.channels.selected_all;
    cfg.sessions(2).notes = 'spont';
    cfg.sessions(3).session_id = 31:35;
    cfg.sessions(3).include = false;
    cfg.sessions(3).selected_channels = cfg.channels.selected_all;
    cfg.sessions(3).notes = 'vspont';
    cfg.bold = struct();
    cfg.bold.data_subfolder = 'roits';
    cfg.bold.input_var = 'roiTs';
    cfg.bold.role_map = struct();
    cfg.bold.role_map.elehp = 'eleHP';
    cfg.bold.role_map.hp = 'HP';

    fprintf(fid, 'CHECKPOINT cfg_built\n');

    params = struct();
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
    if exist('session_filter_ids_override', 'var') && ~isempty(session_filter_ids_override)
        params.session_filter_ids = session_filter_ids_override;
        fprintf(fid, 'SESSION_FILTER_OVERRIDE %s\n', mat2str(session_filter_ids_override));
    end

params.output_root = io_project.get_project_processed_root();
params.save_dir = io_project.get_pipeline_stage_dir( ...
    params.output_root, cfg, 6, 'mua_residual_cross_correlation');
    params.debug_log_file = fullfile(log_root, 'e10gb1_mua_residual_cross_correlation_inner_log.txt');
    params.verbose = false;
    params.progress_every = 50;
    params.top_n_rows = 100;

    fprintf(fid, 'OUTPUT_ROOT %s\n', params.output_root);
    fprintf(fid, 'SAVE_DIR %s\n', params.save_dir);
    fprintf(fid, 'CHECKPOINT before_which\n');
    fprintf(fid, 'WHICH %s\n', which('compute_mua_residual_cross_correlation'));
    fprintf(fid, 'CHECKPOINT before_call\n');

    result = compute_mua_residual_cross_correlation(cfg, params);
    result = save_mua_residual_cross_correlation_result(result, params.save_dir);

    fprintf(fid, 'CHECKPOINT after_call\n');
    fprintf(fid, 'RESULT_MAT %s\n', result.save_paths.main_mat);
    fprintf(fid, 'RESULT_POOLED %s\n', result.save_paths.pooled_corr_csv);
    fprintf(fid, 'RESULT_TOP %s\n', result.save_paths.pooled_top_csv);
    fprintf(fid, 'BLP_BAND_LABEL %s\n', result.source.blp_band_selection_label);
    fprintf(fid, 'BLP_BAND_INDEX %d\n', result.source.blp_band_index);
    fprintf(fid, 'DONE OK\n');
catch ME
    fprintf(fid, 'DONE ERROR\n');
    fprintf(fid, '%s\n', getReport(ME, 'extended', 'hyperlinks', 'off'));
    rethrow(ME);
end
