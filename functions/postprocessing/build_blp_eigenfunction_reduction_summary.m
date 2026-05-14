function summary = build_blp_eigenfunction_reduction_summary(core, summary_cfg)
%BUILD_BLP_EIGENFUNCTION_REDUCTION_SUMMARY Build lightweight summary readouts.

summary = struct();
summary.temporal_components_smooth_time_by_comp = [];

if ~summary_cfg.smooth.enable || isempty(core.temporal_components_time_by_comp)
    return;
end

switch lower(summary_cfg.smooth.method)
    case 'movmean'
        summary.temporal_components_smooth_time_by_comp = movmean( ...
            core.temporal_components_time_by_comp, ...
            summary_cfg.smooth.window, 1, ...
            'Endpoints', 'shrink');
    otherwise
        error('Unsupported summary smoothing method %s.', summary_cfg.smooth.method);
end
end
