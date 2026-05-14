function result_file = find_latest_blp_eigenfunction_reduction_result(processed_root, dataset_cfg)
%FIND_LATEST_BLP_EIGENFUNCTION_REDUCTION_RESULT Find the newest pipeline 5 MAT.

stage_root = io_project.get_pipeline_stage_dir( ...
    processed_root, dataset_cfg, 5, 'eigenfunction_reduction');

if exist(stage_root, 'dir') ~= 7
    error('Pipeline 5 reduction directory does not exist:\n  %s', stage_root);
end

L = dir(fullfile(stage_root, '**', '*.mat'));
L = L(~[L.isdir]);
if isempty(L)
    error('No pipeline 5 reduction MAT files were found under:\n  %s', stage_root);
end

[~, idx] = max([L.datenum]);
result_file = fullfile(L(idx).folder, L(idx).name);
end
