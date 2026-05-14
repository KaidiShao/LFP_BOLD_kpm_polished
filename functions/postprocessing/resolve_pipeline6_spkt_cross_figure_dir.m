function figure_dir = resolve_pipeline6_spkt_cross_figure_dir(params, run_info, run_output_root)
%RESOLVE_PIPELINE6_SPKT_CROSS_FIGURE_DIR Resolve where SPKT xcorr overview figures should go.

if params.spkt_cross_save_under_run
    figure_dir = fullfile(run_output_root, ...
        io_project.get_pipeline_stage_name(6, 'figures_spkt_residual_cross_correlation'));
else
    figure_dir = io_project.get_pipeline_stage_dir( ...
        params.processed_root, run_info.dataset, 6, 'figures_spkt_residual_cross_correlation');
end
end
