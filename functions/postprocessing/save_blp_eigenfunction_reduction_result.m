function result = save_blp_eigenfunction_reduction_result(result, cfg)
%SAVE_BLP_EIGENFUNCTION_REDUCTION_RESULT Save compact or full reduction result artifact.

if exist(cfg.save.dir, 'dir') ~= 7
    mkdir(cfg.save.dir);
end

save_path = local_build_save_path(cfg, result);
result.artifacts.result_mat_file = save_path;
result_to_save = local_build_save_result(result, cfg.save.payload);
save_vars = struct('result', result_to_save);

if cfg.save.v7_3
    save(save_path, '-struct', 'save_vars', 'result', '-v7.3');
else
    save(save_path, '-struct', 'save_vars', 'result');
end
end


function result_to_save = local_build_save_result(result, payload)
if isstring(payload)
    payload = char(payload);
end
payload = lower(strtrim(payload));

switch payload
    case {'full', 'all'}
        result_to_save = result;

    case {'compact', 'minimal'}
        result_to_save = result;
        result_to_save.meta.save_payload = 'compact';

        if isfield(result_to_save, 'data') && isstruct(result_to_save.data)
            result_to_save.data = local_compact_data(result_to_save.data);
        end

        if isfield(result_to_save, 'core') && isstruct(result_to_save.core)
            result_to_save.core = local_compact_core(result_to_save.core);
        end

        if isfield(result_to_save, 'aux') && isstruct(result_to_save.aux)
            result_to_save.aux = local_compact_aux(result_to_save.aux);
        end

        if isfield(result_to_save, 'concat') && isstruct(result_to_save.concat)
            result_to_save.concat = local_compact_concat(result_to_save.concat);
        end

        if isfield(result_to_save, 'cfg') && isstruct(result_to_save.cfg)
            result_to_save.cfg = local_compact_cfg(result_to_save.cfg);
        end

    otherwise
        error('Unknown cfg.save.payload = %s. Use ''compact'' or ''full''.', payload);
end
end


function data = local_compact_data(data)
remove_fields = {'efun_raw_time_by_mode', 'efun_feature_time_by_mode', ...
    'kpm_modes_mode_by_dict'};
omitted = local_existing_fields(data, remove_fields);
data = local_rmfield_if_present(data, remove_fields);
if ~isempty(omitted)
    data.omitted_compact_fields = omitted;
end
end


function core = local_compact_core(core)
remove_fields = {'reconstruction_time_by_mode'};
omitted = local_existing_fields(core, remove_fields);
core = local_rmfield_if_present(core, remove_fields);
if ~isempty(omitted)
    core.omitted_compact_fields = omitted;
end
end


function aux = local_compact_aux(aux)
if isfield(aux, 'spectrum') && isstruct(aux.spectrum)
    remove_fields = {'distance_mode_by_mode'};
    omitted = local_existing_fields(aux.spectrum, remove_fields);
    aux.spectrum = local_rmfield_if_present(aux.spectrum, remove_fields);

    if isfield(aux.spectrum, 'embedding_info') && ...
            isstruct(aux.spectrum.embedding_info)
        aux.spectrum.embedding_info = local_rmfield_if_present( ...
            aux.spectrum.embedding_info, {'gm'});
    end

    if ~isempty(omitted)
        aux.spectrum.omitted_compact_fields = omitted;
    end
end
end


function concat = local_compact_concat(concat)
if isfield(concat, 'chunk_ids') && ~isempty(concat.chunk_ids)
    concat.chunk_id_range = [concat.chunk_ids(1), concat.chunk_ids(end)];
end

remove_fields = {'files', 'chunk_ids', 'chunk_lengths', ...
    'chunk_start_idx', 'chunk_end_idx'};
omitted = local_existing_fields(concat, remove_fields);
concat = local_rmfield_if_present(concat, remove_fields);
if ~isempty(omitted)
    concat.omitted_compact_fields = omitted;
end
end


function cfg = local_compact_cfg(cfg)
cfg = strip_blp_preloaded_eigenfunction_source_payload(cfg);

if isfield(cfg, 'viz') && isstruct(cfg.viz)
    if isfield(cfg.viz, 'spectrum') && isstruct(cfg.viz.spectrum)
        cfg.viz.spectrum = local_rmfield_if_present(cfg.viz.spectrum, ...
            {'distance_colormap', 'embedding_index_colormap', ...
            'cluster_colormap'});
    end

    if isfield(cfg.viz, 'state_space') && isstruct(cfg.viz.state_space)
        cfg.viz.state_space = local_rmfield_if_present(cfg.viz.state_space, ...
            {'time_colormap', 'value_colormap'});
    end
end
end


function save_path = local_build_save_path(cfg, result)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {local_default_file_stem(cfg), cfg.path.kind, lower(cfg.feature.variant), lower(result.meta.method)};

if isfield(cfg.save, 'tag') && ~isempty(cfg.save.tag)
    pieces{end+1} = cfg.save.tag;
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(local_default_save_dir(cfg), sprintf('%s__%s.mat', filename, timestamp));
end


function save_dir = local_default_save_dir(cfg)
if ~isfield(cfg, 'save') || ~isfield(cfg.save, 'dir') || isempty(cfg.save.dir)
    source_run_name = local_source_run_name(cfg);
    dataset_name = local_dataset_name(cfg);
    stage_root = io_project.get_pipeline_stage_dir( ...
        io_project.get_project_processed_root(), dataset_name, 5, 'eigenfunction_reduction');
    save_dir = fullfile(stage_root, source_run_name, 'mat');
else
    save_dir = cfg.save.dir;
end
end


function file_stem = local_default_file_stem(cfg)
if isfield(cfg, 'save') && isfield(cfg.save, 'file_stem') && ~isempty(cfg.save.file_stem)
    file_stem = cfg.save.file_stem;
elseif isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    file_stem = sprintf('%s_efun', cfg.dataset.name);
else
    file_stem = 'efun';
end
end


function source_run_name = local_source_run_name(cfg)
source_run_name = 'unspecified_source';

if isfield(cfg, 'output') && isfield(cfg.output, 'output_run_name') && ...
        ~isempty(cfg.output.output_run_name)
    source_run_name = cfg.output.output_run_name;
elseif isfield(cfg, 'output') && isfield(cfg.output, 'source_run_name') && ...
        ~isempty(cfg.output.source_run_name)
    source_run_name = cfg.output.source_run_name;
elseif isfield(cfg, 'source') && isfield(cfg.source, 'data_dir') && ...
        ~isempty(cfg.source.data_dir)
    [~, source_run_name] = fileparts(cfg.source.data_dir);
elseif isfield(cfg, 'source') && isfield(cfg.source, 'edmd_file') && ...
        ~isempty(cfg.source.edmd_file)
    [~, source_run_name] = fileparts(cfg.source.edmd_file);
end

source_run_name = regexprep(char(source_run_name), '[^\w\-]+', '_');
end


function dataset_name = local_dataset_name(cfg)
dataset_name = 'efun';
if isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    dataset_name = char(cfg.dataset.name);
end
dataset_name = regexprep(dataset_name, '[^\w\-]+', '_');
end


function fields = local_existing_fields(S, names)
fields = {};
for i = 1:numel(names)
    if isfield(S, names{i})
        fields{end+1, 1} = names{i}; %#ok<AGROW>
    end
end
end


function S = local_rmfield_if_present(S, names)
for i = 1:numel(names)
    if isfield(S, names{i})
        S = rmfield(S, names{i});
    end
end
end
