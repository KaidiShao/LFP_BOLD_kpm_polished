function cfg = load_blp_pipeline_cfg_by_stem(dataset_stem)
%LOAD_BLP_PIPELINE_CFG_BY_STEM Resolve cfg_*.m from canonical dataset stem.

switch lower(char(string(dataset_stem)))
    case 'e10fv1'
        cfg = cfg_E10fV1();
    case 'e10gb1'
        cfg = cfg_E10gb1();
    case 'e10gh1'
        cfg = cfg_E10gH1();
    case 'e10gw1'
        cfg = cfg_E10gW1();
    case 'f12m01'
        cfg = cfg_F12m01();
    case 'f12m02'
        cfg = cfg_F12m02();
    case 'f12m03'
        cfg = cfg_F12m03();
    case 'f12m04'
        cfg = cfg_F12m04();
    case 'f12m05'
        cfg = cfg_F12m05();
    case 'k13m17'
        cfg = cfg_K13m17();
    case 'k13m18'
        cfg = cfg_K13m18();
    case 'k13m19'
        cfg = cfg_K13m19();
    case 'k13m20'
        cfg = cfg_K13m20();
    case 'k13m21'
        cfg = cfg_K13m21();
    case 'k13m23'
        cfg = cfg_K13m23();
    otherwise
        error('No dataset config mapping is defined for %s.', char(string(dataset_stem)));
end

if ~isfield(cfg, 'spectrogram')
    cfg.spectrogram = struct();
end
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
end
