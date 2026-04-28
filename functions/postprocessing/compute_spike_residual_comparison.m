function result = compute_spike_residual_comparison(cfg, params)
%COMPUTE_SPIKE_RESIDUAL_COMPARISON Stream Koopman residuals onto spike bins.
%
% This function avoids concatenating the full efun matrix in memory. It
% reads EDMD chunk files one-by-one, computes Koopman residuals for the
% selected modes, aggregates them to the spike time base, and computes
% channel-by-mode zero-lag correlations.

if nargin < 1 || isempty(cfg)
    error('cfg must be provided.');
end

if nargin < 2
    params = struct();
end

params = local_apply_defaults(cfg, params);

[source_dir, chunk_files, source_select_info] = local_resolve_source_dir(cfg, params);
if params.verbose
    fprintf('EDMD chunk source:\n  %s\n', source_dir);
    fprintf('Found %d chunk files.\n', numel(chunk_files));
end

metadata_file = local_resolve_dictionary_file(cfg, params, chunk_files);
if params.verbose
    fprintf('Dictionary metadata file:\n  %s\n', metadata_file);
end

meta = local_load_dictionary_metadata(metadata_file);
Spk = local_load_spike_dataset_for_comparison(cfg, meta, params);
local_validate_metadata_alignment(meta, Spk);

session_mask = true(numel(meta.session_ids), 1);
if ~isempty(params.session_filter_ids)
    session_mask = ismember(meta.session_ids, params.session_filter_ids(:));
    if ~any(session_mask)
        error('None of params.session_filter_ids were found in the metadata session list.');
    end
end

S0 = load(chunk_files(1).fullpath, params.variable_name);
if ~isfield(S0, params.variable_name)
    error('File %s does not contain variable %s.', chunk_files(1).fullpath, params.variable_name);
end
first_outputs = S0.(params.variable_name);

mode_info = local_select_modes(first_outputs.evalues, params);
if params.verbose
    fprintf('Selected %d residual modes after threshold/sort.\n', mode_info.n_modes);
end

session_state = local_initialize_session_state(cfg, meta, Spk, session_mask, mode_info, params);
last_required_global_idx = max([session_state.lfp_end_idx]);

ref_info = local_extract_ref_chunk_info(first_outputs);
cursor = 1;
prev_phi_last = [];
n_chunks = numel(chunk_files);

for k = 1:n_chunks
    if params.verbose && local_should_print_progress(k, n_chunks, params.progress_every)
        fprintf('[stream %d/%d] %s\n', k, n_chunks, chunk_files(k).name);
    end

    if k == 1
        current_outputs = first_outputs;
    else
        S = load(chunk_files(k).fullpath, params.variable_name);
        if ~isfield(S, params.variable_name)
            error('File %s does not contain variable %s.', chunk_files(k).fullpath, params.variable_name);
        end
        current_outputs = S.(params.variable_name);
    end

    local_validate_chunk_against_reference(ref_info, current_outputs, chunk_files(k).name);

    phi_chunk = current_outputs.efuns(:, mode_info.selected_idx);
    if size(phi_chunk, 2) ~= mode_info.n_modes
        error('Chunk %s does not provide the expected selected mode count.', chunk_files(k).name);
    end

    u_chunk = local_compute_chunk_residual( ...
        phi_chunk, mode_info.lambda_d_row, prev_phi_last, params.residual.first_u_mode, cursor == 1);

    session_state = local_accumulate_chunk(session_state, u_chunk, cursor);

    prev_phi_last = phi_chunk(end, :);
    cursor = cursor + size(phi_chunk, 1);

    if cursor > last_required_global_idx
        break;
    end
end

total_streamed_length = cursor - 1;
expected_required_length = max([session_state.lfp_end_idx]);
if total_streamed_length < expected_required_length
    error(['Stopped after %d samples, but the selected sessions require at least %d. ', ...
        'The EDMD source appears incomplete for the chosen sessions.'], ...
        total_streamed_length, expected_required_length);
end

[result_sessions, session_summary] = local_finalize_sessions(cfg, session_state);
mode_table = local_build_mode_table(mode_info);
channel_table = local_build_channel_table(Spk.selected_channels, Spk.channel_sites);
session_corr_table = local_build_session_corr_table(result_sessions, mode_info, channel_table);
pooled_corr_table = local_build_pooled_corr_table(result_sessions, mode_info, channel_table);

if ~isempty(session_corr_table)
    finite_session = isfinite(session_corr_table.abs_corr);
    session_corr_table = [ ...
        sortrows(session_corr_table(finite_session, :), {'abs_corr', 'session_id'}, {'descend', 'ascend'}); ...
        session_corr_table(~finite_session, :)];
end
if ~isempty(pooled_corr_table)
    finite_pooled = isfinite(pooled_corr_table.abs_corr);
    pooled_corr_table = [ ...
        sortrows(pooled_corr_table(finite_pooled, :), 'abs_corr', 'descend'); ...
        pooled_corr_table(~finite_pooled, :)];
end

finite_top = pooled_corr_table(isfinite(pooled_corr_table.abs_corr), :);
pooled_top_corr = finite_top(1:min(height(finite_top), params.top_n_rows), :);

result = struct();
result.cfg = cfg;
result.params = params;
result.source = struct();
result.source.data_dir = source_dir;
result.source.selection = source_select_info;
result.source.spike_channel_selection_label = Spk.channel_selection_label;
result.source.first_chunk = chunk_files(1).fullpath;
result.source.last_chunk = chunk_files(min(n_chunks, k)).fullpath;
result.source.n_chunks_total = n_chunks;
result.source.n_chunks_used = k;
result.source.filename_pattern = params.filename_pattern;
result.source.variable_name = params.variable_name;
result.source.ref_info = ref_info;
result.metadata_file = metadata_file;
result.lfp_selected_channels = meta.selected_channels;
result.lfp_channel_sites = meta.channel_sites;
result.spike_selected_channels = Spk.selected_channels;
result.spike_channel_sites = Spk.channel_sites;
result.total_streamed_length = total_streamed_length;
result.mode_table = mode_table;
result.channel_table = channel_table;
result.session_summary = session_summary;
result.session_corr_table = session_corr_table;
result.pooled_corr_table = pooled_corr_table;
result.pooled_top_corr = pooled_top_corr;
result.sessions = result_sessions;

if params.save_results
    [save_paths, result] = local_save_result(result, params, cfg, source_dir);
else
    save_paths = struct('main_mat', '', 'session_corr_csv', '', 'pooled_corr_csv', '', ...
        'pooled_top_csv', '', 'mode_csv', '', 'session_summary_csv', '');
end

result.save_paths = save_paths;
end


function params = local_apply_defaults(cfg, params)
if ~isfield(params, 'source_cfg') || isempty(params.source_cfg)
    params.source_cfg = struct();
end
if ~isfield(params.source_cfg, 'mode') || isempty(params.source_cfg.mode)
    params.source_cfg.mode = 'chunk_dir';
end
if ~strcmpi(params.source_cfg.mode, 'chunk_dir')
    error('This streamed comparison currently supports only source_cfg.mode = ''chunk_dir''.');
end
if ~isfield(params.source_cfg, 'data_dir')
    params.source_cfg.data_dir = '';
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
    params.save_dir = fullfile(params.output_root, cfg.file_stem, 'spike_residual_comparison');
end
if ~isfield(params, 'save_results')
    params.save_results = true;
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
    selection_info.method = 'explicit';
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

dict_dir = fullfile(params.output_root, cfg.file_stem, 'reskoopnet_dictionary');
prefix = chunk_files(1).prefix;
tokens = regexp(prefix, '^(.*)_Python_', 'tokens', 'once');

candidate_list = {};
if ~isempty(tokens)
    candidate_list{end+1, 1} = fullfile(dict_dir, [tokens{1} '.mat']); %#ok<AGROW>
end
candidate_list{end+1, 1} = fullfile(dict_dir, [cfg.file_stem '_low50_high250_g2_abs_single.mat']); %#ok<AGROW>

for i = 1:numel(candidate_list)
    if exist(candidate_list{i}, 'file') == 2
        metadata_file = candidate_list{i};
        return;
    end
end

L = dir(fullfile(dict_dir, '*.mat'));
if numel(L) == 1
    metadata_file = fullfile(L(1).folder, L(1).name);
    return;
end

error(['Could not resolve the dictionary metadata file automatically.\n' ...
    'Checked candidates under:\n  %s'], dict_dir);
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
            label = 'spk_selected';
            return;
        case 'all'
            label = 'spk_all';
            return;
    end
end

if isequal(spike_channels, 1:max(spike_channels))
    label = sprintf('spk_%dch', numel(spike_channels));
else
    label = sprintf('spk_custom_%dch', numel(spike_channels));
end
end


function local_validate_metadata_alignment(meta, Spk)
if ~isequal(meta.session_ids(:), double(Spk.session_ids(:)))
    error('Spike session IDs do not match the dictionary metadata session IDs.');
end
end


function mode_info = local_select_modes(evalues_in, params)
evalues0 = evalues_in(:);

switch lower(params.post.sort_by)
    case {'modulus', 'abs'}
        key_full = abs(evalues0);
    case {'real', 'realpart'}
        key_full = real(evalues0);
    otherwise
        error('Unknown params.post.sort_by = %s.', params.post.sort_by);
end

[~, ord_full] = sort(key_full, params.post.sort_dir);
evalues_sorted = evalues0(ord_full);
mask_sorted = abs(evalues_sorted) > params.post.abs_thresh;

if ~any(mask_sorted)
    error('Threshold removed all modes. Adjust params.post.abs_thresh.');
end

idx_sorted_in_original = ord_full(mask_sorted);
Kkeep = min(numel(idx_sorted_in_original), params.post.max_basis);
Ksel = min(Kkeep, params.residual.max_modes);

selected_idx = idx_sorted_in_original(1:Ksel);
selected_evalues = evalues0(selected_idx);
lambda_d = local_to_discrete_lambda(selected_evalues, params.residual);

mode_info = struct();
mode_info.selected_idx = selected_idx(:).';
mode_info.selected_evalues = selected_evalues(:);
mode_info.lambda_d = lambda_d(:);
mode_info.lambda_d_row = reshape(lambda_d(:).', 1, []);
mode_info.n_modes = numel(selected_idx);
mode_info.abs_thresh = params.post.abs_thresh;
mode_info.sort_by = params.post.sort_by;
mode_info.sort_dir = params.post.sort_dir;
mode_info.max_basis = params.post.max_basis;
mode_info.max_modes = params.residual.max_modes;
end


function lambda_d = local_to_discrete_lambda(lambda_in, residual_cfg)
switch lower(residual_cfg.lambdaType)
    case 'discrete'
        lambda_d = lambda_in;
    case 'continuous'
        lambda_d = exp(lambda_in .* residual_cfg.dt);
    otherwise
        error('Unknown params.residual.lambdaType = %s.', residual_cfg.lambdaType);
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


function local_validate_chunk_against_reference(ref_info, current_outputs, file_name)
required_fields = {'efuns', 'evalues'};
for i = 1:numel(required_fields)
    if ~isfield(current_outputs, required_fields{i})
        error('Field %s is missing from %s.', required_fields{i}, file_name);
    end
end

if ~isequaln(ref_info.evalues, current_outputs.evalues)
    error('evalues changed in chunk file %s.', file_name);
end

if ~isempty(ref_info.kpm_modes)
    if ~isfield(current_outputs, 'kpm_modes') || ~isequaln(ref_info.kpm_modes, current_outputs.kpm_modes)
        error('kpm_modes changed in chunk file %s.', file_name);
    end
end

if ~isempty(ref_info.N_dict)
    if ~isfield(current_outputs, 'N_dict') || ~isequaln(ref_info.N_dict, current_outputs.N_dict)
        error('N_dict changed in chunk file %s.', file_name);
    end
end

if ~isempty(ref_info.residual_form)
    if ~isfield(current_outputs, 'residual_form') || ~isequaln(ref_info.residual_form, current_outputs.residual_form)
        error('residual_form changed in chunk file %s.', file_name);
    end
end
end


function session_state = local_initialize_session_state(cfg, meta, Spk, session_mask, mode_info, params)
session_ids = meta.session_ids(session_mask);
session_lengths = meta.session_lengths(session_mask);
session_dx = meta.session_dx(session_mask);
session_start_idx = meta.session_start_idx(session_mask);
session_end_idx = meta.session_end_idx(session_mask);

spike_lengths = double(Spk.session_lengths(session_mask));
spike_dx = double(Spk.session_dx(session_mask));
spike_data_cells = Spk.data_cells(session_mask);

n_sessions = numel(session_ids);
session_state = repmat(struct( ...
    'session_id', [], ...
    'session_note', '', ...
    'lfp_start_idx', [], ...
    'lfp_end_idx', [], ...
    'lfp_length', [], ...
    'lfp_dx', [], ...
    'spike_dx', [], ...
    'n_spike_bins_raw', [], ...
    'n_spike_bins_used', [], ...
    'lfp_duration_sec', [], ...
    'spike_duration_sec', [], ...
    'common_duration_sec', [], ...
    'lfp_to_spike_bin', [], ...
    'spike_counts', [], ...
    'spike_rate', [], ...
    'bin_counts', [], ...
    'sum_abs', [], ...
    'sum_abs_sq', [], ...
    'sum_real', [], ...
    'sum_imag', []), n_sessions, 1);

for i = 1:n_sessions
    n_lfp = session_lengths(i);
    dx_lfp = session_dx(i);
    dx_spk = spike_dx(i);
    n_spk_raw = spike_lengths(i);

    n_spk_from_lfp = floor((n_lfp * dx_lfp) / dx_spk);
    n_spk_used = min(n_spk_raw, n_spk_from_lfp);
    if n_spk_used < 1
        error('Session %d has zero overlapping spike bins.', session_ids(i));
    end

    lfp_to_spike_bin = floor(((0:n_lfp-1) .* dx_lfp) ./ dx_spk) + 1;
    lfp_to_spike_bin = int32(lfp_to_spike_bin(:));
    lfp_to_spike_bin(lfp_to_spike_bin > n_spk_used) = 0;

    spike_counts = single(spike_data_cells{i}(1:n_spk_used, :));
    spike_rate = single(double(spike_counts) ./ dx_spk);

    session_state(i).session_id = session_ids(i);
    session_state(i).session_note = local_lookup_session_note(cfg, session_ids(i));
    session_state(i).lfp_start_idx = session_start_idx(i);
    session_state(i).lfp_end_idx = session_end_idx(i);
    session_state(i).lfp_length = n_lfp;
    session_state(i).lfp_dx = dx_lfp;
    session_state(i).spike_dx = dx_spk;
    session_state(i).n_spike_bins_raw = n_spk_raw;
    session_state(i).n_spike_bins_used = n_spk_used;
    session_state(i).lfp_duration_sec = n_lfp * dx_lfp;
    session_state(i).spike_duration_sec = n_spk_raw * dx_spk;
    session_state(i).common_duration_sec = n_spk_used * dx_spk;
    session_state(i).lfp_to_spike_bin = lfp_to_spike_bin;
    session_state(i).spike_counts = spike_counts;
    session_state(i).spike_rate = spike_rate;
    session_state(i).bin_counts = single(zeros(n_spk_used, 1));
    session_state(i).sum_abs = single(zeros(n_spk_used, mode_info.n_modes));
    session_state(i).sum_abs_sq = single(zeros(n_spk_used, mode_info.n_modes));
    session_state(i).sum_real = single(zeros(n_spk_used, mode_info.n_modes));
    session_state(i).sum_imag = single(zeros(n_spk_used, mode_info.n_modes));
end

if params.verbose
    fprintf('Prepared %d session accumulators for spike-bin aggregation.\n', n_sessions);
end
end


function note = local_lookup_session_note(cfg, session_id)
note = '';
for i = 1:numel(cfg.sessions)
    if any(double(cfg.sessions(i).session_id(:)) == double(session_id))
        if isfield(cfg.sessions(i), 'notes') && ~isempty(cfg.sessions(i).notes)
            note = char(string(cfg.sessions(i).notes));
        end
        return;
    end
end
end


function u_chunk = local_compute_chunk_residual(phi_chunk, lambda_d_row, prev_phi_last, first_u_mode, is_first_chunk)
[T, K] = size(phi_chunk);
u_chunk = zeros(T, K, 'like', phi_chunk);

if T < 1
    return;
end

if is_first_chunk
    switch lower(first_u_mode)
        case 'phi1'
            u_chunk(1, :) = phi_chunk(1, :);
        case 'zero'
            u_chunk(1, :) = 0;
        otherwise
            error('Unknown params.residual.first_u_mode = %s.', first_u_mode);
    end
else
    u_chunk(1, :) = phi_chunk(1, :) - prev_phi_last .* lambda_d_row;
end

if T >= 2
    u_chunk(2:end, :) = phi_chunk(2:end, :) - phi_chunk(1:end-1, :) .* lambda_d_row;
end
end


function session_state = local_accumulate_chunk(session_state, u_chunk, global_start_idx)
global_end_idx = global_start_idx + size(u_chunk, 1) - 1;

for i = 1:numel(session_state)
    seg_start = max(global_start_idx, session_state(i).lfp_start_idx);
    seg_end = min(global_end_idx, session_state(i).lfp_end_idx);
    if seg_start > seg_end
        continue;
    end

    chunk_rows = (seg_start:seg_end) - global_start_idx + 1;
    local_lfp_rows = (seg_start:seg_end) - session_state(i).lfp_start_idx + 1;

    spike_bin_idx = double(session_state(i).lfp_to_spike_bin(local_lfp_rows));
    keep = spike_bin_idx > 0;
    if ~any(keep)
        continue;
    end

    spike_bin_idx = spike_bin_idx(keep);
    u_seg = u_chunk(chunk_rows(keep), :);

    B = sparse(spike_bin_idx, 1:numel(spike_bin_idx), 1, ...
        session_state(i).n_spike_bins_used, numel(spike_bin_idx));

    session_state(i).bin_counts = session_state(i).bin_counts + ...
        single(accumarray(spike_bin_idx(:), 1, [session_state(i).n_spike_bins_used, 1], @sum, 0));

    % MATLAB sparse matrices do not support mtimes with single arrays.
    % Accumulate in double for the sparse multiply, then store compactly.
    abs_u = double(abs(u_seg));
    session_state(i).sum_abs = session_state(i).sum_abs + single(B * abs_u);
    session_state(i).sum_abs_sq = session_state(i).sum_abs_sq + single(B * (abs_u .^ 2));
    session_state(i).sum_real = session_state(i).sum_real + single(B * double(real(u_seg)));
    session_state(i).sum_imag = session_state(i).sum_imag + single(B * double(imag(u_seg)));
end
end


function [result_sessions, session_summary] = local_finalize_sessions(cfg, session_state)
n_sessions = numel(session_state);
result_sessions = repmat(struct( ...
    'session_id', [], ...
    'session_note', '', ...
    'lfp_dx', [], ...
    'spike_dx', [], ...
    'n_lfp_samples', [], ...
    'n_spike_bins_raw', [], ...
    'n_spike_bins_used', [], ...
    'lfp_duration_sec', [], ...
    'spike_duration_sec', [], ...
    'common_duration_sec', [], ...
    'spike_counts', [], ...
    'spike_rate', [], ...
    'lfp_samples_per_spike_bin', [], ...
    'residual_abs_mean', [], ...
    'residual_abs_rms', [], ...
    'residual_real_mean', [], ...
    'residual_imag_mean', [], ...
    'corr', []), n_sessions, 1);

session_id_col = zeros(n_sessions, 1);
session_note_col = strings(n_sessions, 1);
lfp_dx_col = zeros(n_sessions, 1);
spike_dx_col = zeros(n_sessions, 1);
n_lfp_col = zeros(n_sessions, 1);
n_spk_raw_col = zeros(n_sessions, 1);
n_spk_used_col = zeros(n_sessions, 1);
lfp_dur_col = zeros(n_sessions, 1);
spk_dur_col = zeros(n_sessions, 1);
common_dur_col = zeros(n_sessions, 1);

for i = 1:n_sessions
    counts = double(session_state(i).bin_counts);
    count_safe = counts;
    count_safe(count_safe == 0) = NaN;

    residual_abs_mean = single(session_state(i).sum_abs ./ count_safe);
    residual_abs_rms = single(sqrt(session_state(i).sum_abs_sq ./ count_safe));
    residual_real_mean = single(session_state(i).sum_real ./ count_safe);
    residual_imag_mean = single(session_state(i).sum_imag ./ count_safe);

    corr_struct = struct();
    corr_struct.abs_mean = local_corr_matrix(session_state(i).spike_rate, residual_abs_mean);
    corr_struct.abs_rms = local_corr_matrix(session_state(i).spike_rate, residual_abs_rms);
    corr_struct.real_mean = local_corr_matrix(session_state(i).spike_rate, residual_real_mean);
    corr_struct.imag_mean = local_corr_matrix(session_state(i).spike_rate, residual_imag_mean);

    result_sessions(i).session_id = session_state(i).session_id;
    result_sessions(i).session_note = session_state(i).session_note;
    result_sessions(i).lfp_dx = session_state(i).lfp_dx;
    result_sessions(i).spike_dx = session_state(i).spike_dx;
    result_sessions(i).n_lfp_samples = session_state(i).lfp_length;
    result_sessions(i).n_spike_bins_raw = session_state(i).n_spike_bins_raw;
    result_sessions(i).n_spike_bins_used = session_state(i).n_spike_bins_used;
    result_sessions(i).lfp_duration_sec = session_state(i).lfp_duration_sec;
    result_sessions(i).spike_duration_sec = session_state(i).spike_duration_sec;
    result_sessions(i).common_duration_sec = session_state(i).common_duration_sec;
    result_sessions(i).spike_counts = session_state(i).spike_counts;
    result_sessions(i).spike_rate = session_state(i).spike_rate;
    result_sessions(i).lfp_samples_per_spike_bin = session_state(i).bin_counts;
    result_sessions(i).residual_abs_mean = residual_abs_mean;
    result_sessions(i).residual_abs_rms = residual_abs_rms;
    result_sessions(i).residual_real_mean = residual_real_mean;
    result_sessions(i).residual_imag_mean = residual_imag_mean;
    result_sessions(i).corr = corr_struct;

    session_id_col(i) = session_state(i).session_id;
    session_note_col(i) = string(session_state(i).session_note);
    lfp_dx_col(i) = session_state(i).lfp_dx;
    spike_dx_col(i) = session_state(i).spike_dx;
    n_lfp_col(i) = session_state(i).lfp_length;
    n_spk_raw_col(i) = session_state(i).n_spike_bins_raw;
    n_spk_used_col(i) = session_state(i).n_spike_bins_used;
    lfp_dur_col(i) = session_state(i).lfp_duration_sec;
    spk_dur_col(i) = session_state(i).spike_duration_sec;
    common_dur_col(i) = session_state(i).common_duration_sec;
end

session_summary = table( ...
    session_id_col, ...
    session_note_col, ...
    lfp_dx_col, ...
    spike_dx_col, ...
    n_lfp_col, ...
    n_spk_raw_col, ...
    n_spk_used_col, ...
    lfp_dur_col, ...
    spk_dur_col, ...
    common_dur_col, ...
    'VariableNames', { ...
        'session_id', 'session_note', 'lfp_dx_sec', 'spike_dx_sec', ...
        'n_lfp_samples', 'n_spike_bins_raw', 'n_spike_bins_used', ...
        'lfp_duration_sec', 'spike_duration_sec', 'common_duration_sec'});
end


function mode_table = local_build_mode_table(mode_info)
mode_rank = (1:mode_info.n_modes).';
mode_original_idx = double(mode_info.selected_idx(:));
evalue_real = real(mode_info.selected_evalues(:));
evalue_imag = imag(mode_info.selected_evalues(:));
evalue_abs = abs(mode_info.selected_evalues(:));
lambda_d_real = real(mode_info.lambda_d(:));
lambda_d_imag = imag(mode_info.lambda_d(:));

mode_table = table( ...
    mode_rank, mode_original_idx, ...
    evalue_real, evalue_imag, evalue_abs, ...
    lambda_d_real, lambda_d_imag, ...
    'VariableNames', { ...
        'mode_rank', 'mode_original_idx', ...
        'evalue_real', 'evalue_imag', 'evalue_abs', ...
        'lambda_d_real', 'lambda_d_imag'});
end


function channel_table = local_build_channel_table(selected_channels, channel_sites)
channel_rank = (1:numel(selected_channels)).';
channel_index = double(selected_channels(:));
channel_site = string(channel_sites(:));
channel_label = strings(numel(selected_channels), 1);
for i = 1:numel(selected_channels)
    channel_label(i) = sprintf('ch%02d_%s', channel_index(i), char(channel_site(i)));
end

channel_table = table(channel_rank, channel_index, channel_site, channel_label, ...
    'VariableNames', {'channel_rank', 'channel_index', 'channel_site', 'channel_label'});
end


function session_corr_table = local_build_session_corr_table(result_sessions, mode_info, channel_table)
tables = cell(numel(result_sessions), 1);
for i = 1:numel(result_sessions)
    tables{i} = local_build_corr_table_from_struct( ...
        result_sessions(i).corr, ...
        result_sessions(i).session_id, ...
        result_sessions(i).session_note, ...
        result_sessions(i).n_spike_bins_used, ...
        mode_info, channel_table);
end

tables = tables(~cellfun(@isempty, tables));
if isempty(tables)
    session_corr_table = table();
else
    session_corr_table = vertcat(tables{:});
end
end


function pooled_corr_table = local_build_pooled_corr_table(result_sessions, mode_info, channel_table)
if isempty(result_sessions)
    pooled_corr_table = table();
    return;
end

pooled_spike = cat(1, result_sessions.spike_rate);
pooled_abs_mean = cat(1, result_sessions.residual_abs_mean);
pooled_abs_rms = cat(1, result_sessions.residual_abs_rms);
pooled_real_mean = cat(1, result_sessions.residual_real_mean);
pooled_imag_mean = cat(1, result_sessions.residual_imag_mean);

corr_struct = struct();
corr_struct.abs_mean = local_corr_matrix(pooled_spike, pooled_abs_mean);
corr_struct.abs_rms = local_corr_matrix(pooled_spike, pooled_abs_rms);
corr_struct.real_mean = local_corr_matrix(pooled_spike, pooled_real_mean);
corr_struct.imag_mean = local_corr_matrix(pooled_spike, pooled_imag_mean);

pooled_corr_table = local_build_corr_table_from_struct( ...
    corr_struct, NaN, 'pooled', size(pooled_spike, 1), mode_info, channel_table);
end


function T = local_build_corr_table_from_struct(corr_struct, session_id, session_note, n_bins, mode_info, channel_table)
feature_names = fieldnames(corr_struct);
tables = cell(numel(feature_names), 1);

for f = 1:numel(feature_names)
    feature_name = feature_names{f};
    C = corr_struct.(feature_name);
    [n_channels, n_modes] = size(C);
    if n_channels == 0 || n_modes == 0
        tables{f} = table();
        continue;
    end

    [channel_grid, mode_grid] = ndgrid(1:n_channels, 1:n_modes);
    n_rows = numel(C);

    session_id_col = repmat(session_id, n_rows, 1);
    session_note_col = repmat(string(session_note), n_rows, 1);
    residual_feature_col = repmat(string(feature_name), n_rows, 1);
    n_bins_col = repmat(n_bins, n_rows, 1);

    mode_rank = mode_grid(:);
    selected_idx_col = mode_info.selected_idx(:);
    selected_evalues_col = mode_info.selected_evalues(:);
    mode_original_idx = double(selected_idx_col(mode_rank));
    evalue_real = real(selected_evalues_col(mode_rank));
    evalue_imag = imag(selected_evalues_col(mode_rank));
    evalue_abs = abs(selected_evalues_col(mode_rank));

    ch_rank = channel_grid(:);
    channel_index = channel_table.channel_index(ch_rank);
    channel_site = channel_table.channel_site(ch_rank);
    channel_label = channel_table.channel_label(ch_rank);

    corr_val = C(:);
    abs_corr = abs(corr_val);

    tables{f} = table( ...
        session_id_col, ...
        session_note_col, ...
        channel_index, ...
        channel_site, ...
        channel_label, ...
        mode_rank, ...
        mode_original_idx, ...
        residual_feature_col, ...
        evalue_real, ...
        evalue_imag, ...
        evalue_abs, ...
        corr_val, ...
        abs_corr, ...
        n_bins_col, ...
        'VariableNames', { ...
            'session_id', 'session_note', ...
            'channel_index', 'channel_site', 'channel_label', ...
            'mode_rank', 'mode_original_idx', 'residual_feature', ...
            'evalue_real', 'evalue_imag', 'evalue_abs', ...
            'corr', 'abs_corr', 'n_bins'});
end

tables = tables(~cellfun(@isempty, tables));
if isempty(tables)
    T = table();
else
    T = vertcat(tables{:});
end
end


function C = local_corr_matrix(X, Y)
X = double(X);
Y = double(Y);

if isempty(X) || isempty(Y)
    C = zeros(size(X, 2), size(Y, 2));
    return;
end

valid = all(isfinite(X), 2) & all(isfinite(Y), 2);
X = X(valid, :);
Y = Y(valid, :);

if size(X, 1) < 2
    C = nan(size(X, 2), size(Y, 2));
    return;
end

X = X - mean(X, 1);
Y = Y - mean(Y, 1);

sx = std(X, 0, 1);
sy = std(Y, 0, 1);
sx(sx == 0) = NaN;
sy(sy == 0) = NaN;

C = (X.' * Y) ./ ((size(X, 1) - 1) .* (sx(:) * sy(:).'));
end


function [save_paths, result] = local_save_result(result, params, cfg, source_dir)
save_dir = params.save_dir;
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

source_label = local_filename_safe(local_basename(source_dir));
spike_label = local_filename_safe(result.source.spike_channel_selection_label);
tag = sprintf('%s_%s_%s_koopman_residual_edmd_top%d', ...
    cfg.file_stem, source_label, spike_label, height(result.mode_table));
tag = local_filename_safe(tag);

save_paths = struct();
save_paths.main_mat = fullfile(save_dir, [tag '.mat']);
save_paths.session_corr_csv = fullfile(save_dir, [tag '_session_corr.csv']);
save_paths.pooled_corr_csv = fullfile(save_dir, [tag '_pooled_corr.csv']);
save_paths.pooled_top_csv = fullfile(save_dir, [tag '_pooled_top.csv']);
save_paths.mode_csv = fullfile(save_dir, [tag '_modes.csv']);
save_paths.session_summary_csv = fullfile(save_dir, [tag '_session_summary.csv']);

result.save_paths = save_paths;
tmp_mat = [save_paths.main_mat, '.tmp.mat'];
if exist(tmp_mat, 'file') == 2
    delete(tmp_mat);
end
save(tmp_mat, 'result', '-v7.3');
if exist(save_paths.main_mat, 'file') == 2
    delete(save_paths.main_mat);
end
movefile(tmp_mat, save_paths.main_mat, 'f');
if ~isempty(result.session_corr_table)
    writetable(result.session_corr_table, save_paths.session_corr_csv);
else
    writetable(table(), save_paths.session_corr_csv);
end
if ~isempty(result.pooled_corr_table)
    writetable(result.pooled_corr_table, save_paths.pooled_corr_csv);
    writetable(result.pooled_top_corr, save_paths.pooled_top_csv);
else
    writetable(table(), save_paths.pooled_corr_csv);
    writetable(table(), save_paths.pooled_top_csv);
end
writetable(result.mode_table, save_paths.mode_csv);
writetable(result.session_summary, save_paths.session_summary_csv);
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


function tf = local_should_print_progress(k, n_total, progress_every)
tf = (k == 1) || (k == n_total) || (mod(k, progress_every) == 0);
end


function out = local_filename_safe(name_in)
out = regexprep(char(name_in), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled';
end
end


function name = local_basename(path_in)
[~, name, ext] = fileparts(path_in);
name = [name, ext];
end
