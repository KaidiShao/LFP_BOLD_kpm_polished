function summary_root = get_summary_figures_root(output_root)
%GET_SUMMARY_FIGURES_ROOT Resolve the shared processed summary-figure root.

if nargin < 1 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

summary_root = fullfile(output_root, 'summary_figures');
end
