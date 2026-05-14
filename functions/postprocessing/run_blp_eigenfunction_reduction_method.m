function [reduction, method_name] = run_blp_eigenfunction_reduction_method(prep, cfg)
%RUN_BLP_EIGENFUNCTION_REDUCTION_METHOD Run one configured reduction method.

switch lower(cfg.path.kind)
    case 'time'
        method_name = cfg.path.time.method;
        reduction = reduce_eigenfunction_time_path( ...
            prep.efun_feature_time_by_mode, cfg.path.time);

    case 'spectrum'
        spectrum_cfg = cfg.path.spectrum;
        spectrum_cfg.evalues_discrete = prep.evalues_discrete;
        spectrum_cfg.evalues_bilinear = prep.evalues_bilinear;
        spectrum_cfg.feature_variant = cfg.feature.variant;
        method_name = cfg.path.spectrum.method;
        reduction = reduce_eigenfunction_spectrum_path( ...
            prep.efun_feature_time_by_mode, spectrum_cfg);

    otherwise
        error('Unknown cfg.path.kind = %s. Use ''time'' or ''spectrum''.', cfg.path.kind);
end
end
