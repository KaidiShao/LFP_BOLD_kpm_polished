function save_dir = resolve_pipeline6_spkt_cross_save_dir(params, run_info, run_output_root)
%RESOLVE_PIPELINE6_SPKT_CROSS_SAVE_DIR Resolve where SPKT xcorr data should go.

if params.spkt_cross_save_under_run
    save_dir = fullfile(run_output_root, params.spkt_cross_save_dir_name);
else
    save_dir = io_project.get_pipeline_stage_dir( ...
        params.processed_root, run_info.dataset, 6, 'spkt_residual_cross_correlation');
end
end
