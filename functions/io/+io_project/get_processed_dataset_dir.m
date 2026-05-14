function dataset_dir = get_processed_dataset_dir(output_root, cfg_or_stem)
%GET_PROCESSED_DATASET_DIR Resolve the canonical processed directory for one dataset.

if nargin < 1 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 2 || isempty(cfg_or_stem)
    error('cfg_or_stem is required.');
end

if isstruct(cfg_or_stem)
    if ~isfield(cfg_or_stem, 'file_stem') || isempty(cfg_or_stem.file_stem)
        error('cfg_or_stem.file_stem is required.');
    end
    file_stem = cfg_or_stem.file_stem;
else
    file_stem = cfg_or_stem;
end

file_stem = char(string(file_stem));
dataset_dir = fullfile(output_root, file_stem);
end
