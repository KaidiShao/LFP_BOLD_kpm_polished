function window_root = get_window_figures_root(output_root)
%GET_WINDOW_FIGURES_ROOT Resolve the shared processed top-window figure root.

if nargin < 1 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

window_root = fullfile(output_root, 'window_figures');
end
