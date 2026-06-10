function cfg = strip_blp_preloaded_eigenfunction_source_payload(cfg)
%STRIP_BLP_PRELOADED_EIGENFUNCTION_SOURCE_PAYLOAD Remove large preloaded source fields.

if isfield(cfg, 'source') && isstruct(cfg.source)
    remove_fields = {'preloaded_EDMD_outputs', 'preloaded_concat_info', ...
        'preloaded_source_info'};
    for i = 1:numel(remove_fields)
        if isfield(cfg.source, remove_fields{i})
            cfg.source = rmfield(cfg.source, remove_fields{i});
        end
    end
end
end
