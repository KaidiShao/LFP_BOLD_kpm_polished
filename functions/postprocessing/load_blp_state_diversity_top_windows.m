function cache = load_blp_state_diversity_top_windows(cfg, params)
%LOAD_BLP_STATE_DIVERSITY_TOP_WINDOWS Load or compute saved state-diversity windows.

window_file = fullfile( ...
    io_project.get_pipeline_stage_dir(params.processed_root, cfg, 2, 'consensus_state_diversity_windows'), ...
    sprintf('%s_consensus_state_diversity_windows_%dsamp_globalwin.mat', ...
    cfg.file_stem, params.window_length_samples));

if exist(window_file, 'file') ~= 2
    fprintf('State-diversity window file missing; computing it:\n  %s\n', window_file);
    loader_cfg = struct('file_stem', cfg.file_stem);
    [C, source_consensus_file] = io_results.load_consensus_state_results(loader_cfg, params.processed_root, []);
    win_params = struct();
    win_params.window_length_samples = params.window_length_samples;
    win_params.window_mode = params.window_mode;
    win_params.keep_partial_window = false;
    win_params.top_k = params.n_top_windows;
    win_params.save_csv = true;
    win_params.force_recompute = false;
    W = analyze_blp_consensus_state_diversity_windows( ...
        cfg, params.processed_root, C, win_params, source_consensus_file);
    window_file = W.save_file;
else
    S = load(window_file, 'W');
    if ~isfield(S, 'W')
        error('State-diversity file does not contain variable W: %s', window_file);
    end
    W = S.W;
end

local_validate_top_window_result(W, window_file, params);
cache = struct('W', W, 'window_file', window_file);
end


function local_validate_top_window_result(W, window_file, params)
if ~isfield(W, 'top_windows_table') || isempty(W.top_windows_table)
    error('Top-window table is missing or empty: %s', window_file);
end
if ~isfield(W, 'window_mode') || ~strcmpi(char(string(W.window_mode)), params.window_mode)
    error('Top-window file is not a %s-window result: %s', params.window_mode, window_file);
end
if ~isfield(W, 'window_length_samples') || ...
        double(W.window_length_samples) ~= double(params.window_length_samples)
    error('Top-window file does not use %d-sample windows: %s', ...
        params.window_length_samples, window_file);
end

required_vars = {'state_diversity_rank', 'global_start_idx', 'global_end_idx'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table is missing required column "%s": %s', ...
            required_vars{i}, window_file);
    end
end
end
