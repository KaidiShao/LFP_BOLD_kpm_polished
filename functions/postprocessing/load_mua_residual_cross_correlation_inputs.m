function [ctx, params] = load_mua_residual_cross_correlation_inputs(cfg, params)
%LOAD_MUA_RESIDUAL_CROSS_CORRELATION_INPUTS Load MUA xcorr inputs and context.

params = local_apply_defaults(cfg, params);

ctx = struct();
ctx.cfg = cfg;
ctx.params = params;
ctx.source_mode = lower(strtrim(char(string(params.source_cfg.mode))));
local_debug_log(params, 'after_apply_defaults');
ctx.pair_specs = local_resolve_pair_specs(params);
local_debug_log(params, 'after_resolve_pair_specs');

[ctx.source_dir, ctx.chunk_files, ctx.source_select_info] = local_resolve_source_dir(cfg, params);
local_debug_log(params, sprintf('after_resolve_source_dir n_chunks=%d', numel(ctx.chunk_files)));
if params.verbose
    fprintf('EDMD chunk source:\n  %s\n', ctx.source_dir);
    fprintf('Found %d chunk files.\n', numel(ctx.chunk_files));
end

ctx.metadata_file = local_resolve_dictionary_file(cfg, params, ctx.chunk_files);
local_debug_log(params, 'after_resolve_dictionary_file');
if params.verbose
    fprintf('Dictionary metadata file:\n  %s\n', ctx.metadata_file);
end

ctx.meta = local_load_dictionary_metadata(ctx.metadata_file);
local_debug_log(params, 'after_load_dictionary_metadata');
ctx.Blp = local_load_blp_dataset_for_comparison(cfg, ctx.meta, params);
local_debug_log(params, 'after_load_blp_dataset_for_comparison');
local_validate_metadata_alignment(ctx.meta, ctx.Blp);
local_debug_log(params, 'after_validate_metadata_alignment');

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
    local_debug_log(params, 'after_load_first_chunk');
else
    S0 = load(ctx.chunk_files(1).fullpath, params.variable_name);
    if ~isfield(S0, params.variable_name)
        error('File %s does not contain variable %s.', ctx.chunk_files(1).fullpath, params.variable_name);
    end
    ctx.first_outputs = S0.(params.variable_name);
    local_debug_log(params, 'after_load_first_chunk');
    ctx.mode_info = build_blp_residual_mode_info( ...
        ctx.first_outputs.evalues, build_blp_residual_cfg_from_cross_params(params));
    ctx.ref_info = local_extract_ref_chunk_info(ctx.first_outputs);
end

local_debug_log(params, sprintf('after_select_modes n_modes=%d', ctx.mode_info.n_modes));
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
if ~isfield(params, 'blp_channels') || isempty(params.blp_channels)
    params.blp_channels = 'selected';
end
if ~isfield(params, 'blp_band') || isempty(params.blp_band)
    params.blp_band = 'last';
end
if ~isfield(params, 'pairings') || isempty(params.pairings)
    params.pairings = {'abs_abs', 'raw_real', 'raw_imag'};
end
params.pairings = cellstr(string(params.pairings(:)).');
if ~isfield(params, 'dictionary_file')
    params.dictionary_file = '';
end
if ~isfield(params, 'output_root') || isempty(params.output_root)
    params.output_root = io_project.get_project_processed_root();
end
if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = io_project.get_pipeline_stage_dir( ...
        params.output_root, cfg, 6, 'mua_residual_cross_correlation');
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


function pair_specs = local_resolve_pair_specs(params)
pairings = cellstr(string(params.pairings(:)).');
pair_specs = repmat(struct( ...
    'name', '', ...
    'signal_feature', '', ...
    'residual_feature', ''), numel(pairings), 1);

for i = 1:numel(pairings)
    name = lower(strtrim(pairings{i}));
    switch name
        case 'abs_abs'
            signal_feature = 'abs';
            residual_feature = 'abs';
        case 'raw_real'
            signal_feature = 'raw';
            residual_feature = 'real';
        case 'raw_imag'
            signal_feature = 'raw';
            residual_feature = 'imag';
        otherwise
            error('Unsupported params.pairings entry: %s', pairings{i});
    end

    pair_specs(i).name = name;
    pair_specs(i).signal_feature = signal_feature;
    pair_specs(i).residual_feature = residual_feature;
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


function Blp = local_load_blp_dataset_for_comparison(cfg, meta, params)
cfg_blp = cfg;
blp_channels = local_resolve_blp_channels(cfg, meta, params.blp_channels);

for i = 1:numel(cfg_blp.sessions)
    if cfg_blp.sessions(i).include
        cfg_blp.sessions(i).selected_channels = blp_channels;
    end
end

[all_session_ids, all_session_channels] = ...
    io_utils.collect_included_session_entries(cfg_blp, blp_channels);
n_sessions = numel(all_session_ids);

data_dir = fullfile(cfg_blp.raw_data_root, cfg_blp.data_subfolder);
data_cells = cell(n_sessions, 1);
session_lengths = zeros(n_sessions, 1);
session_dx = zeros(n_sessions, 1);
channel_sites = cfg.channels.sites(blp_channels);

band_index = [];
n_bands = [];

for k = 1:n_sessions
    sid = all_session_ids(k);
    this_channels = all_session_channels{k};
    fname = fullfile(data_dir, sprintf('%s_%04d_blp.mat', cfg_blp.file_stem, sid));
    if k <= 2
        local_debug_log(params, sprintf('blp_loader before_load session=%d file=%s', sid, fname));
    end
    S = load(fname, 'blp');
    if k <= 2
        local_debug_log(params, sprintf('blp_loader after_load session=%d', sid));
    end
    if ~isfield(S, 'blp')
        error('File %s does not contain variable blp.', fname);
    end
    if ~isfield(S.blp, 'dat') || ~isfield(S.blp, 'dx')
        error('File %s is missing blp.dat or blp.dx.', fname);
    end

    dat = S.blp.dat;
    if ndims(dat) < 3
        error('Expected blp.dat in %s to have a 3rd dimension.', fname);
    end

    if isempty(n_bands)
        n_bands = size(dat, 3);
        band_index = local_resolve_blp_band_index(params.blp_band, n_bands);
    elseif size(dat, 3) ~= n_bands
        error('BLP band count changed across sessions. Expected %d bands but found %d in %s.', ...
            n_bands, size(dat, 3), fname);
    end

    x = dat(:, this_channels, band_index);
    x = reshape(x, size(dat, 1), numel(this_channels));
    if k <= 2
        local_debug_log(params, sprintf('blp_loader after_slice session=%d size=%d x %d band=%d', ...
            sid, size(x, 1), size(x, 2), band_index));
    end
    data_cells{k} = single(x);
    session_lengths(k) = size(x, 1);
    session_dx(k) = double(S.blp.dx);
end

Blp = struct();
Blp.data_cells = data_cells;
Blp.session_ids = double(all_session_ids(:));
Blp.session_lengths = double(session_lengths(:));
Blp.session_dx = double(session_dx(:));
Blp.selected_channels = double(blp_channels(:)).';
Blp.channel_sites = cellstr(string(channel_sites(:)));
Blp.channel_selection_label = local_describe_blp_channel_selection(params.blp_channels, blp_channels);
Blp.band_index = band_index;
Blp.n_bands = n_bands;
Blp.band_selection_label = local_describe_blp_band_selection(params.blp_band, band_index, n_bands);
end


function blp_channels = local_resolve_blp_channels(cfg, meta, blp_channel_spec)
if ischar(blp_channel_spec) || isstring(blp_channel_spec)
    spec = lower(strtrim(char(string(blp_channel_spec))));
    switch spec
        case 'selected'
            blp_channels = double(meta.selected_channels(:)).';
        case 'all'
            blp_channels = 1:numel(cfg.channels.sites);
        otherwise
            error('Unknown params.blp_channels = %s. Use ''selected'', ''all'', or a numeric vector.', spec);
    end
elseif isnumeric(blp_channel_spec)
    blp_channels = double(blp_channel_spec(:)).';
else
    error('params.blp_channels must be ''selected'', ''all'', or a numeric vector.');
end

blp_channels = unique(round(blp_channels), 'stable');
if isempty(blp_channels)
    error('params.blp_channels resolved to an empty channel list.');
end

n_total = numel(cfg.channels.sites);
if any(blp_channels < 1) || any(blp_channels > n_total)
    error('params.blp_channels contains indices outside 1:%d.', n_total);
end
end


function band_index = local_resolve_blp_band_index(spec_in, n_bands)
if ischar(spec_in) || isstring(spec_in)
    spec = lower(strtrim(char(string(spec_in))));
    switch spec
        case 'last'
            band_index = n_bands;
        otherwise
            error('Unknown params.blp_band = %s. Use ''last'' or a numeric index.', spec);
    end
elseif isnumeric(spec_in) && isscalar(spec_in) && isfinite(spec_in) && spec_in == floor(spec_in)
    band_index = double(spec_in);
else
    error('params.blp_band must be ''last'' or a positive integer.');
end

if band_index < 1 || band_index > n_bands
    error('Resolved BLP band index %d is outside 1:%d.', band_index, n_bands);
end
end


function label = local_describe_blp_channel_selection(blp_channel_spec, blp_channels)
if ischar(blp_channel_spec) || isstring(blp_channel_spec)
    spec = lower(strtrim(char(string(blp_channel_spec))));
    switch spec
        case 'selected'
            label = 'mua_selected';
            return;
        case 'all'
            label = 'mua_all';
            return;
    end
end

if isequal(blp_channels, 1:max(blp_channels))
    label = sprintf('mua_%dch', numel(blp_channels));
else
    label = sprintf('mua_custom_%dch', numel(blp_channels));
end
end


function label = local_describe_blp_band_selection(blp_band_spec, band_index, n_bands)
if ischar(blp_band_spec) || isstring(blp_band_spec)
    spec = lower(strtrim(char(string(blp_band_spec))));
    if strcmp(spec, 'last')
        label = sprintf('band_last_idx%dof%d', band_index, n_bands);
        return;
    end
end
label = sprintf('band_%02dof%d', band_index, n_bands);
end


function local_validate_metadata_alignment(meta, Blp)
if ~isequal(meta.session_ids(:), double(Blp.session_ids(:)))
    error('BLP session IDs do not match the dictionary metadata session IDs.');
end
if ~isequal(meta.session_lengths(:), double(Blp.session_lengths(:)))
    error('BLP session lengths do not match the dictionary metadata session lengths.');
end
if any(abs(meta.session_dx(:) - double(Blp.session_dx(:))) > 1e-12)
    error('BLP session dx does not match the dictionary metadata session dx.');
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


function local_debug_log(params, message)
if ~isfield(params, 'debug_log_file') || isempty(params.debug_log_file)
    return;
end

fid = fopen(params.debug_log_file, 'a');
if fid < 0
    return;
end
cleanup_obj = onCleanup(@() fclose(fid)); %#ok<NASGU>
fprintf(fid, '[%s] %s\n', char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss')), char(message));
end
