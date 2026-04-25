function [EDMD_outputs, concat_info, source_info] = load_edmd_source(source_cfg)
%LOAD_EDMD_SOURCE Load EDMD outputs from chunks, one MAT file, or memory.

if nargin < 1 || isempty(source_cfg)
    error('source_cfg must be provided.');
end

if ~isfield(source_cfg, 'mode') || isempty(source_cfg.mode)
    error('source_cfg.mode must be provided.');
end

switch lower(source_cfg.mode)
    case 'preloaded'
        if ~isfield(source_cfg, 'preloaded_EDMD_outputs') || ...
                isempty(source_cfg.preloaded_EDMD_outputs)
            error('source_cfg.preloaded_EDMD_outputs must be provided for preloaded mode.');
        end

        EDMD_outputs = source_cfg.preloaded_EDMD_outputs;
        if isfield(source_cfg, 'preloaded_concat_info')
            concat_info = source_cfg.preloaded_concat_info;
        else
            concat_info = [];
        end

        if isfield(source_cfg, 'preloaded_source_info') && ...
                ~isempty(source_cfg.preloaded_source_info)
            source_info = source_cfg.preloaded_source_info;
        else
            source_info = struct();
        end
        source_info.cache_mode = 'preloaded';

    case 'chunk_dir'
        if ~isfield(source_cfg, 'data_dir') || isempty(source_cfg.data_dir)
            error('source_cfg.data_dir must be provided for chunk_dir mode.');
        end
        if ~isfield(source_cfg, 'concat') || isempty(source_cfg.concat)
            source_cfg.concat = struct();
        end

        [EDMD_outputs, concat_info] = load_and_concat_edmd_output_chunks( ...
            source_cfg.data_dir, source_cfg.concat);

        source_info = struct();
        source_info.mode = 'chunk_dir';
        source_info.path = source_cfg.data_dir;

    case 'mat_file'
        edmd_file = '';
        if isfield(source_cfg, 'edmd_file') && ~isempty(source_cfg.edmd_file)
            edmd_file = source_cfg.edmd_file;
        elseif isfield(source_cfg, 'data_dir') && ~isempty(source_cfg.data_dir)
            edmd_file = local_find_default_edmd_file(source_cfg.data_dir);
        else
            error('Provide source_cfg.edmd_file or source_cfg.data_dir for mat_file mode.');
        end

        S = load(edmd_file);
        if ~isfield(S, 'EDMD_outputs')
            error('File %s does not contain variable EDMD_outputs.', edmd_file);
        end

        EDMD_outputs = S.EDMD_outputs;
        if isfield(S, 'concat_info')
            concat_info = S.concat_info;
        else
            concat_info = [];
        end

        source_info = struct();
        source_info.mode = 'mat_file';
        source_info.path = edmd_file;

    otherwise
        error(['Unknown source_cfg.mode = %s. Use ''chunk_dir'', ', ...
            '''mat_file'', or ''preloaded''.'], source_cfg.mode);
end
end


function edmd_file = local_find_default_edmd_file(data_dir)
patterns = {'*_concat.mat', '*_outputs_*_to_*_concat.mat'};

for i = 1:numel(patterns)
    L = dir(fullfile(data_dir, patterns{i}));
    if isempty(L)
        continue;
    end

    [~, order] = sort([L.datenum], 'descend');
    L = L(order);
    edmd_file = fullfile(L(1).folder, L(1).name);
    return;
end

error(['No concatenated EDMD MAT file was found in %s. ', ...
    'Set source_cfg.edmd_file explicitly or use source_cfg.mode = ''chunk_dir'' instead.'], data_dir);
end
