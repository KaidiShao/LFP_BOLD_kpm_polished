function sync_one_pipeline7_postprocessing_summary_figures(run_info, result, output_root)
%SYNC_ONE_PIPELINE7_POSTPROCESSING_SUMMARY_FIGURES Copy pipeline 7 PNGs into shared summaries.

if nargin < 2 || ~isstruct(result)
    return;
end
if nargin < 3 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', local_get_field(result, 'dataset_stem', ''))));
run_name = char(string(local_get_field(run_info, 'run_name', local_get_field(result, 'run_name', ''))));
if isempty(dataset_stem) || isempty(run_name)
    return;
end

copy_one(result, 'main_png', output_root, dataset_stem, ...
    io_project.get_pipeline_summary_dir(output_root, 7, 'figures_bold_reskoopnet_main_summary'));
copy_one(result, 'deconv_png', output_root, dataset_stem, ...
    io_project.get_pipeline_summary_dir(output_root, 7, 'figures_bold_reskoopnet_deconv_summary'));
copy_one(result, 'timescale_png', output_root, dataset_stem, ...
    io_project.get_pipeline_summary_dir(output_root, 7, 'figures_bold_reskoopnet_timescale_summary'));
end


function copy_one(result, field_name, output_root, dataset_stem, summary_root) %#ok<INUSD>
file_path = char(string(local_get_field(result, field_name, '')));
if isempty(file_path) || exist(file_path, 'file') ~= 2
    return;
end

dataset_dir = fullfile(summary_root, dataset_stem);
if exist(dataset_dir, 'dir') ~= 7
    mkdir(dataset_dir);
end

[~, name, ext] = fileparts(file_path);
copyfile(file_path, fullfile(dataset_dir, [name ext]), 'f');
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
