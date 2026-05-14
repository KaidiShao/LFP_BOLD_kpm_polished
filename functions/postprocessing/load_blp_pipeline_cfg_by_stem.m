function cfg = load_blp_pipeline_cfg_by_stem(dataset_stem)
%LOAD_BLP_PIPELINE_CFG_BY_STEM Resolve cfg_*.m from canonical dataset stem.

switch lower(char(string(dataset_stem)))
    case 'e10fv1'
        cfg = cfg_E10fV1();
    case 'e10gb1'
        cfg = cfg_E10gb1();
    case 'e10gh1'
        cfg = cfg_E10gH1();
    case 'f12m01'
        cfg = cfg_F12m01();
    otherwise
        error('No dataset config mapping is defined for %s.', char(string(dataset_stem)));
end

if ~isfield(cfg, 'spectrogram')
    cfg.spectrogram = struct();
end
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
end
