function stage_name = get_pipeline_stage_name(pipeline_idx, stage_key)
%GET_PIPELINE_STAGE_NAME Resolve the canonical folder name for one pipeline stage.

stage_key = char(string(stage_key));

switch double(pipeline_idx)
    case 1
        switch lower(stage_key)
            case 'spectrograms'
                stage_name = 'pipeline1_spectrograms';
            case 'dictionary'
                stage_name = 'pipeline1_reskoopnet_dictionary';
            otherwise
                error('Unsupported pipeline 1 stage_key: %s', stage_key);
        end

    case 2
        switch lower(stage_key)
            case 'event_detection'
                stage_name = 'pipeline2_event_detection';
            case 'event_density'
                stage_name = 'pipeline2_event_density';
            case 'legacy_event_density'
                stage_name = 'pipeline2_legacy_event_density';
            case 'figures_legacy_event_density'
                stage_name = 'pipeline2_figures_legacy_event_density';
            case 'consensus_states'
                stage_name = 'pipeline2_consensus_states';
            case 'consensus_state_summary'
                stage_name = 'pipeline2_consensus_state_summary';
            case 'event_diversity_windows'
                stage_name = 'pipeline2_event_diversity_windows';
            case 'consensus_state_diversity_windows'
                stage_name = 'pipeline2_consensus_state_diversity_windows';
            case 'figures_event_diversity_top_window_plots'
                stage_name = 'pipeline2_figures_event_diversity_top_window_plots';
            case 'figures_consensus_state_top_window_plots'
                stage_name = 'pipeline2_figures_consensus_state_top_window_plots';
            otherwise
                error('Unsupported pipeline 2 stage_key: %s', stage_key);
        end

    case 3
        switch lower(stage_key)
            case 'bold_observables'
                stage_name = 'pipeline3_bold_observables';
            case 'figures_bold_pre_reskoopnet_qc'
                stage_name = 'pipeline3_figures_bold_pre_reskoopnet_qc';
            case 'figures_bold_pre_reskoopnet_qc_summary'
                stage_name = 'pipeline3_figures_bold_pre_reskoopnet_qc_summary';
            otherwise
                error('Unsupported pipeline 3 stage_key: %s', stage_key);
        end

    case 5
        switch lower(stage_key)
            case 'eigenfunction_reduction'
                stage_name = 'pipeline5_eigenfunction_reduction';
            case 'efun_dimred_top30'
                stage_name = 'pipeline5_efun_dimred_top30';
            case 'raw_thresholded_density'
                stage_name = 'pipeline5_raw_thresholded_density';
            case 'raw_thresholded_density_scan'
                stage_name = 'pipeline5_raw_thresholded_density_scan';
            case 'raw_thresholded_events'
                stage_name = 'pipeline5_raw_thresholded_events';
            case 'dimred_thresholded_density'
                stage_name = 'pipeline5_dimred_thresholded_density';
            case 'dimred_thresholded_events'
                stage_name = 'pipeline5_dimred_thresholded_events';
            case 'eigenfunction_peaks_by_state'
                stage_name = 'pipeline5_eigenfunction_peaks_by_state';
            case 'figures_eigenfunction_reduction_ssc_summary'
                stage_name = 'pipeline5_figures_eigenfunction_reduction_ssc_summary';
            otherwise
                error('Unsupported pipeline 5 stage_key: %s', stage_key);
        end

    case 6
        switch lower(stage_key)
            case 'top_state_diversity_postprocessing'
                stage_name = 'pipeline6_top_state_diversity_postprocessing';
            case 'spkt_residual_cross_correlation'
                stage_name = 'pipeline6_spkt_residual_cross_correlation';
            case 'spkt_residual_lagged_cross_correlation'
                stage_name = 'pipeline6_spkt_residual_lagged_cross_correlation';
            case 'mua_residual_cross_correlation'
                stage_name = 'pipeline6_mua_residual_cross_correlation';
            case 'figures_timescale_diagnostics'
                stage_name = 'pipeline6_figures_timescale_diagnostics';
            case 'figures_spkt_residual_cross_correlation'
                stage_name = 'pipeline6_figures_spkt_residual_cross_correlation';
            case 'figures_spkt_residual_lagged_cross_correlation'
                stage_name = 'pipeline6_figures_spkt_residual_lagged_cross_correlation';
            case 'figures_mua_residual_cross_correlation'
                stage_name = 'pipeline6_figures_mua_residual_cross_correlation';
            otherwise
                error('Unsupported pipeline 6 stage_key: %s', stage_key);
        end

    case 7
        switch lower(stage_key)
            case 'bold_postprocessing'
                stage_name = 'pipeline7_bold_reskoopnet_postprocessing';
            case 'figures_bold_reskoopnet_main_summary'
                stage_name = 'pipeline7_figures_bold_reskoopnet_main_summary';
            case 'figures_bold_reskoopnet_deconv_summary'
                stage_name = 'pipeline7_figures_bold_reskoopnet_deconv_summary';
            case 'figures_bold_reskoopnet_timescale_summary'
                stage_name = 'pipeline7_figures_bold_reskoopnet_timescale_summary';
            otherwise
                error('Unsupported pipeline 7 stage_key: %s', stage_key);
        end

    case 8
        switch lower(stage_key)
            case 'efun_density_cross_correlation'
                stage_name = 'pipeline8_xcorr';
            case 'figures_bold_top_xcorr_activation_maps'
                stage_name = 'pipeline8_top_maps';
            case 'figures_bold_xcorr_summary'
                stage_name = 'pipeline8_xcorr_summary';
            case 'figures_top_xcorr_activation_map'
                stage_name = 'pipeline8_top_maps_summary';
            otherwise
                error('Unsupported pipeline 8 stage_key: %s', stage_key);
        end

    case 9
        switch lower(stage_key)
            case 'bold_eigenfunction_reduction'
                stage_name = 'pipeline9_bold_eigenfunction_reduction';
            case 'figures_bold_eigenfunction_reduction'
                stage_name = 'pipeline9_figures_bold_eigenfunction_reduction';
            otherwise
                error('Unsupported pipeline 9 stage_key: %s', stage_key);
        end

    case 10
        switch lower(stage_key)
            case 'bold_dimred_density_cross_correlation'
                stage_name = 'pipeline10_dimred_xcorr';
            case 'figures_bold_dimred_top_xcorr_activation_maps'
                stage_name = 'pipeline10_dimred_top_maps';
            case 'figures_bold_dimred_xcorr_summary'
                stage_name = 'pipeline10_dimred_xcorr_summary';
            case 'figures_dimred_top_xcorr_activation_map'
                stage_name = 'pipeline10_dimred_top_maps_summary';
            otherwise
                error('Unsupported pipeline 10 stage_key: %s', stage_key);
        end

    otherwise
        error('Unsupported pipeline index: %d', pipeline_idx);
end
end
