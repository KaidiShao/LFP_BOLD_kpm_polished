function summary_dir = get_pipeline_summary_dir(output_root, pipeline_idx, stage_key)
%GET_PIPELINE_SUMMARY_DIR Resolve the canonical shared summary directory.

if nargin < 1 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

summary_root = io_project.get_summary_figures_root(output_root);
stage_name = io_project.get_pipeline_stage_name(pipeline_idx, stage_key);
summary_dir = fullfile(summary_root, stage_name);
end
