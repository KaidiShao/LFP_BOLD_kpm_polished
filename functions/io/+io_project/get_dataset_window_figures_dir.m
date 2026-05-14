function dataset_dir = get_dataset_window_figures_dir(output_root, cfg_or_stem)
%GET_DATASET_WINDOW_FIGURES_DIR Resolve the shared top-window figure dir for one dataset.

window_root = io_project.get_window_figures_root(output_root);

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

dataset_dir = fullfile(window_root, char(string(file_stem)));
end
