this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg = cfg_eigenfunction_reduction_minimal();
[result, EDMD_outputs, concat_info, source_info] = run_eigenfunction_reduction_pipeline(cfg);

fprintf('Loaded EDMD source mode: %s\n', source_info.mode);
fprintf('Input size: T=%d, N=%d\n', ...
    size(result.data.efun_feature_time_by_mode, 1), ...
    size(result.data.efun_feature_time_by_mode, 2));
fprintf('Path kind: %s\n', result.meta.path_kind);
fprintf('Method: %s\n', result.meta.method);
fprintf('Temporal components size: [%s]\n', ...
    num2str(size(result.core.temporal_components_time_by_comp)));
fprintf('Mode weights size: [%s]\n', ...
    num2str(size(result.core.mode_weights_mode_by_comp)));

if isfield(result.artifacts, 'result_mat_file') && ~isempty(result.artifacts.result_mat_file)
    fprintf('Saved result to:\n  %s\n', result.artifacts.result_mat_file);
end

if ~isempty(concat_info)
    fprintf('Concatenated %d chunks into %d samples.\n', ...
        concat_info.n_chunks, concat_info.total_length);
end

clear EDMD_outputs concat_info source_info
