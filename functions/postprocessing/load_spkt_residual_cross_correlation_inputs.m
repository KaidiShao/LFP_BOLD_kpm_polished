function [ctx, params] = load_spkt_residual_cross_correlation_inputs(cfg, params)
%LOAD_SPKT_RESIDUAL_CROSS_CORRELATION_INPUTS Load SPKT xcorr inputs and context.

params = local_apply_defaults(cfg, params);

ctx = struct();
ctx.cfg = cfg;
ctx.params = params;
ctx.source_mode = lower(strtrim(char(string(params.source_cfg.mode))));

[ctx.source_dir, ctx.chunk_files, ctx.source_select_info] = local_resolve_source_dir(cfg, params);
if params.verbose
    fprintf('EDMD chunk source:\n  %s\n', ctx.source_dir);
    fprintf('Found %d chunk files.\n', numel(ctx.chunk_files));
end

ctx.metadata_file = local_resolve_dictionary_file(cfg, params, ctx.chunk_files);
if params.verbose
    fprintf('Dictionary metadata file:\n  %s\n', ctx.metadata_file);
end

ctx.meta = local_load_dictionary_metadata(ctx.metadata_file);
ctx.Spk = local_load_spike_dataset_for_comparison(cfg, ctx.meta, params);
local_validate_metadata_alignment(ctx.meta, ctx.Spk);

ctx.session_mask = true(numel(ctx.meta.session_ids), 1);
if ~isempty(params.session_filter_ids)
    ctx.session_mask = ismember(ctx.meta.session_ids, params.session_filter_ids(:));
    if ~any(ctx.session_mask)
        error('None of params.session_filter_ids were found in the metadata session list.');
    end
end

ctx.residual_bundle = [];
ctx.first_outputs = [];
if strcmp(ctx.source_mode, 'residual_workspace')
    if ~isfield(params.source_cfg, 'residual_bundle') || isempty(params.source_cfg.residual_bundle)
        error('params.source_cfg.residual_bundle must be provided for source_cfg.mode = ''residual_workspace''.');
    end
    ctx.residual_bundle = params.source_cfg.residual_bundle;
    ctx.mode_info = subset_blp_residual_bundle(ctx.residual_bundle, params.residual.max_modes);
    ctx.ref_info = build_blp_residual_ref_info(ctx.residual_bundle);
else
    S0 = load(ctx.chunk_files(1).fullpath, params.variable_name);
    if ~isfield(S0, params.variable_name)
        error('File %s does not contain variable %s.', ctx.chunk_files(1).fullpath, params.variable_name);
    end
    ctx.first_outputs = S0.(params.variable_name);
    ctx.mode_info = build_blp_residual_mode_info( ...
        ctx.first_outputs.evalues, build_blp_residual_cfg_from_cross_params(params));
    ctx.ref_info = local_extract_ref_chunk_info(ctx.first_outputs);
end

if params.verbose
    fprintf('Selected %d residual modes after threshold/sort.\n', ctx.mode_info.n_modes);
end
end


function params = local_apply_defaults(cfg, params)
if ~isfield(params, 'source_cfg') || isempty(params.source_cfg)
    params.source_cfg = struct();
end
if ~isfield(params.source_cfg, 'mode') || isempty(params.source_cfg.mode)
    params.source_cfg.mode = 'chunk_dir';
end
if ~strcmpi(params.source_cfg.mode, 'chunk_dir') && ~strcmpi(params.source_cfg.mode, 'residual_workspace')
    error(['Unsupported params.source_cfg.mode = %s. ', ...
        'Use ''chunk_dir'' or ''residual_workspace''.'], params.source_cfg.mode);
end
if ~isfield(params.source_cfg, 'data_dir')
    params.source_cfg.data_dir = '';
end
if ~isfield(params.source_cfg, 'residual_bundle')
    params.source_cfg.residual_bundle = [];
end

if ~isfield(params, 'source') || isempty(params.source)
    params.source = struct();
end
if ~isfield(params.source, 'base_dir') || isempty(params.source.base_dir)
    params.source.base_dir = fullfile('E:\autodl_results', cfg.file_stem, 'mlp', 'outputs');
end
if ~isfield(params.source, 'name_contains') || isempty(params.source.name_contains)
    params.source.name_contains = 'projected_kv_abs';
end
if ~isfield(params.source, 'prefer_non_smoke')
    params.source.prefer_non_smoke = true;
end

if ~isfield(params, 'filename_pattern') || isempty(params.filename_pattern)
    params.filename_pattern = '*_outputs_*.mat';
end
if ~isfield(params, 'variable_name') || isempty(params.variable_name)
    params.variable_name = 'EDMD_outputs';
end

if ~isfield(params, 'post') || isempty(params.post)
    params.post = struct();
end
if ~isfield(params.post, 'abs_thresh') || isempty(params.post.abs_thresh)
    params.post.abs_thresh = 0.01;
end
if ~isfield(params.post, 'sort_by') || isempty(params.post.sort_by)
    params.post.sort_by = 'modulus';
end
if ~isfield(params.post, 'sort_dir') || isempty(params.post.sort_dir)
    params.post.sort_dir = 'descend';
end
if ~isfield(params.post, 'max_basis') || isempty(params.post.max_basis)
    params.post.max_basis = 30;
end

if ~isfield(params, 'residual') || isempty(params.residual)
    params.residual = struct();
end
if ~isfield(params.residual, 'method') || isempty(params.residual.method)
    params.residual.method = 'koopman_residual';
end
if ~strcmpi(params.residual.method, 'koopman_residual')
    error('Only params.residual.method = ''koopman_residual'' is supported here.');
end
if ~isfield(params.residual, 'lambda_source') || isempty(params.residual.lambda_source)
    params.residual.lambda_source = 'edmd';
end
if ~strcmpi(params.residual.lambda_source, 'edmd')
    error('Only params.residual.lambda_source = ''edmd'' is supported in streamed mode.');
end
if ~isfield(params.residual, 'lambdaType') || isempty(params.residual.lambdaType)
    params.residual.lambdaType = 'discrete';
end
if ~isfield(params.residual, 'dt') || isempty(params.residual.dt)
    params.residual.dt = 1;
end
if ~isfield(params.residual, 'first_u_mode') || isempty(params.residual.first_u_mode)
    params.residual.first_u_mode = 'phi1';
end
if ~isfield(params.residual, 'max_modes') || isempty(params.residual.max_modes)
    params.residual.max_modes = 20;
end

if ~isfield(params, 'session_filter_ids')
    params.session_filter_ids = [];
end
if ~isfield(params, 'spike_channels') || isempty(params.spike_channels)
    params.spike_channels = 'selected';
end
if ~isfield(params, 'dictionary_file')
    params.dictionary_file = '';
end
if ~isfield(params, 'output_root') || isempty(params.output_root)
    params.output_root = io_project.get_project_processed_root();
end
if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = io_project.get_pipeline_stage_dir( ...
        params.output_root, cfg, 6, 'spkt_residual_cross_correlation');
end
if ~isfield(params, 'verbose')
    params.verbose = true;
end
if ~isfield(params, 'progress_every') || isempty(params.progress_every)
    params.progress_every = 50;
end
if ~isfield(params, 'top_n_rows') || isempty(params.top_n_rows)
    params.top_n_rows = 200;
end
end


function [source_dir, files, selection_info] = local_resolve_source_dir(cfg, params)
selection_info = struct();

if isfield(params.source_cfg, 'data_dir') && ~isempty(params.source_cfg.data_dir)
    source_dir = params.source_cfg.data_dir;
    if strcmpi(params.source_cfg.mode, 'residual_workspace')
        selection_info.method = 'explicit_residual_workspace';
    else
        selection_info.method = 'explicit';
    end
    selection_info.base_dir = '';
    selection_info.name_contains = '';
    selection_info.prefer_non_smoke = false;
else
    base_dir = params.source.base_dir;
    if exist(base_dir, 'dir') ~= 7
        error('Auto source base directory does not exist: %s', base_dir);
    end

    L = dir(base_dir);
    is_dir = [L.isdir];
    names = {L.name};
    keep = is_dir & ~ismember(names, {'.', '..'});
    L = L(keep);

    if isempty(L)
        error('No candidate EDMD output directories were found in %s.', base_dir);
    end

    if ~isempty(params.source.name_contains)
        mask = contains({L.name}, params.source.name_contains, 'IgnoreCase', true);
        L = L(mask);
    end

    if isempty(L)
        error('No EDMD output directories in %s matched "%s".', ...
            base_dir, params.source.name_contains);
    end

    if params.source.prefer_non_smoke
        nonsmoke_mask = ~contains({L.name}, 'smoke', 'IgnoreCase', true);
        if any(nonsmoke_mask)
            L = L(nonsmoke_mask);
        end
    end

    [~, order] = sort([L.datenum], 'descend');
    L = L(order);
    source_dir = fullfile(L(1).folder, L(1).name);

    selection_info.method = 'auto_latest';
    selection_info.base_dir = base_dir;
    selection_info.name_contains = params.source.name_contains;
    selection_info.prefer_non_smoke = params.source.prefer_non_smoke;
    selection_info.selected_name = L(1).name;
end

files = local_collect_chunk_files(source_dir, params.filename_pattern);
end


function metadata_file = local_resolve_dictionary_file(cfg, params, chunk_files)
if isfield(params, 'dictionary_file') && ~isempty(params.dictionary_file)
    metadata_file = params.dictionary_file;
    if exist(metadata_file, 'file') ~= 2
        error('Dictionary file does not exist: %s', metadata_file);
    end
    return;
end

dict_dir = io_project.get_pipeline_stage_dir(params.output_root, cfg, 1, 'dictionary');
prefix = chunk_files(1).prefix;
tokens = regexp(prefix, '^(.*)_Python_', 'tokens', 'once');

candidate_list = {};
if ~isempty(tokens)
    candidate_list{end+1, 1} = fullfile(dict_dir, [tokens{1} '.mat']); %#ok<AGROW>
end
candidate_list{end+1, 1} = fullfile(dict_dir, [cfg.file_stem '_low50_high250_g2_abs_single.mat']); %#ok<AGROW>

for i = 1:numel(candidate_list)
    if exist(candidate_list{i}, 'file') == 2 && local_dictionary_file_has_required_metadata(candidate_list{i})
        metadata_file = candidate_list{i};
        return;
    end
end

L = dir(fullfile(dict_dir, '*.mat'));
if numel(L) == 1
    candidate_file = fullfile(L(1).folder, L(1).name);
    if local_dictionary_file_has_required_metadata(candidate_file)
        metadata_file = candidate_file;
        return;
    end
end

error(['Could not resolve the dictionary metadata file automatically.\n' ...
    'Checked candidates under:\n  %s'], dict_dir);
end


function tf = local_dictionary_file_has_required_metadata(metadata_file)
required_fields = { ...
    'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', ...
    'selected_channels', 'channel_sites'};

tf = false;
try
    vars = {whos('-file', metadata_file).name};
catch
    return;
end

tf = all(ismember(required_fields, vars));
end


function meta = local_load_dictionary_metadata(metadata_file)
required_fields = { ...
    'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', ...
    'selected_channels', 'channel_sites'};

S = load(metadata_file, required_fields{:});
for i = 1:numel(required_fields)
    name = required_fields{i};
    if ~isfield(S, name)
        error('Field %s is missing from dictionary metadata file %s.', name, metadata_file);
    end
end

meta = struct();
meta.session_ids = double(S.session_ids(:));
meta.session_lengths = double(S.session_lengths(:));
meta.session_dx = double(S.session_dx(:));
meta.session_start_idx = double(S.session_start_idx(:));
meta.session_end_idx = double(S.session_end_idx(:));
meta.selected_channels = double(S.selected_channels(:)).';
meta.channel_sites = cellstr(string(S.channel_sites(:)));
meta.total_length = sum(meta.session_lengths);
end


function Spk = local_load_spike_dataset_for_comparison(cfg, meta, params)
cfg_spk = cfg;
spike_channels = local_resolve_spike_channels(cfg, meta, params.spike_channels);

for i = 1:numel(cfg_spk.sessions)
    if cfg_spk.sessions(i).include
        cfg_spk.sessions(i).selected_channels = spike_channels;
    end
end

Spk = io_raw.load_spike_dataset(cfg_spk);
Spk.channel_selection_label = local_describe_spike_channel_selection(params.spike_channels, spike_channels);
end


function spike_channels = local_resolve_spike_channels(cfg, meta, spike_channel_spec)
if ischar(spike_channel_spec) || isstring(spike_channel_spec)
    spec = lower(strtrim(char(string(spike_channel_spec))));
    switch spec
        case 'selected'
            spike_channels = double(meta.selected_channels(:)).';
        case 'all'
            spike_channels = 1:numel(cfg.channels.sites);
        otherwise
            error('Unknown params.spike_channels = %s. Use ''selected'', ''all'', or a numeric vector.', spec);
    end
elseif isnumeric(spike_channel_spec)
    spike_channels = double(spike_channel_spec(:)).';
else
    error('params.spike_channels must be ''selected'', ''all'', or a numeric vector.');
end

spike_channels = unique(round(spike_channels), 'stable');
if isempty(spike_channels)
    error('params.spike_channels resolved to an empty channel list.');
end

n_total = numel(cfg.channels.sites);
if any(spike_channels < 1) || any(spike_channels > n_total)
    error('params.spike_channels contains indices outside 1:%d.', n_total);
end
end


function label = local_describe_spike_channel_selection(spike_channel_spec, spike_channels)
if ischar(spike_channel_spec) || isstring(spike_channel_spec)
    spec = lower(strtrim(char(string(spike_channel_spec))));
    switch spec
        case 'selected'
            label = 'spkt_selected';
            return;
        case 'all'
            label = 'spkt_all';
            return;
    end
end

if isequal(spike_channels, 1:max(spike_channels))
    label = sprintf('spkt_%dch', numel(spike_channels));
else
    label = sprintf('spkt_custom_%dch', numel(spike_channels));
end
end


function local_validate_metadata_alignment(meta, Spk)
if ~isequal(meta.session_ids(:), double(Spk.session_ids(:)))
    error('Spike session IDs do not match the dictionary metadata session IDs.');
end
end


function ref_info = local_extract_ref_chunk_info(first_outputs)
ref_info = struct();
ref_info.evalues = first_outputs.evalues;
if isfield(first_outputs, 'kpm_modes')
    ref_info.kpm_modes = first_outputs.kpm_modes;
else
    ref_info.kpm_modes = [];
end
if isfield(first_outputs, 'N_dict')
    ref_info.N_dict = first_outputs.N_dict;
else
    ref_info.N_dict = [];
end
if isfield(first_outputs, 'residual_form')
    ref_info.residual_form = first_outputs.residual_form;
else
    ref_info.residual_form = [];
end
end


function files = local_collect_chunk_files(data_dir, filename_pattern)
L = dir(fullfile(data_dir, filename_pattern));
if isempty(L)
    error('No files matching %s were found in %s.', filename_pattern, data_dir);
end

files = repmat(struct('name', '', 'fullpath', '', 'chunk_id', [], 'prefix', ''), numel(L), 1);
keep_count = 0;

for i = 1:numel(L)
    tokens = regexp(L(i).name, '^(.*)_outputs_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tokens)
        continue;
    end

    keep_count = keep_count + 1;
    files(keep_count).name = L(i).name;
    files(keep_count).fullpath = fullfile(L(i).folder, L(i).name);
    files(keep_count).chunk_id = str2double(tokens{2});
    files(keep_count).prefix = tokens{1};
end

files = files(1:keep_count);
if isempty(files)
    error('No files in %s matched the required *_outputs_<chunk>.mat pattern.', data_dir);
end

[~, order] = sort([files.chunk_id]);
files = files(order);

prefixes = {files.prefix};
if numel(unique(prefixes)) > 1
    error('Multiple EDMD output prefixes were found in %s. Narrow the filename pattern first.', data_dir);
end
end
