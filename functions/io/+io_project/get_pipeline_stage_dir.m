function stage_dir = get_pipeline_stage_dir(output_root, cfg_or_stem, pipeline_idx, stage_key)
%GET_PIPELINE_STAGE_DIR Resolve the canonical processed directory for one pipeline stage.

if nargin < 1 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

dataset_dir = io_project.get_processed_dataset_dir(output_root, cfg_or_stem);
stage_name = io_project.get_pipeline_stage_name(pipeline_idx, stage_key);
stage_dir = fullfile(dataset_dir, stage_name);
end
