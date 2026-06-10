function [session_state, total_streamed_length, k_used] = align_mua_residual_to_lfp_resolution(ctx)
%ALIGN_MUA_RESIDUAL_TO_LFP_RESOLUTION Align residuals to native LFP/MUA grid.

session_state = local_initialize_session_state(ctx.cfg, ctx.meta, ctx.Blp, ctx.session_mask, ctx.pair_specs);
local_debug_log(ctx.params, sprintf('after_initialize_session_state n_sessions=%d', numel(session_state)));
last_required_global_idx = max([session_state.lfp_end_idx]);
n_chunks = numel(ctx.chunk_files);
k_used = 0;

if strcmp(ctx.source_mode, 'residual_workspace')
    u_full = ctx.residual_bundle.u_full(:, 1:ctx.mode_info.n_modes);
    local_debug_log(ctx.params, sprintf('after_use_workspace_residual size=%d x %d', size(u_full, 1), size(u_full, 2)));
    session_state = local_accumulate_chunk(session_state, ctx.pair_specs, u_full, 1);
    local_debug_log(ctx.params, 'after_accumulate_workspace_residual');
    total_streamed_length = size(u_full, 1);
    k_used = n_chunks;
    return;
end

local_debug_log(ctx.params, 'after_extract_ref_chunk_info');
cursor = 1;
prev_phi_last = [];
for k = 1:n_chunks
    if ctx.params.verbose && local_should_print_progress(k, n_chunks, ctx.params.progress_every)
        fprintf('[stream %d/%d] %s\n', k, n_chunks, ctx.chunk_files(k).name);
    end

    if k == 1
        current_outputs = ctx.first_outputs;
    else
        S = load(ctx.chunk_files(k).fullpath, ctx.params.variable_name);
        if ~isfield(S, ctx.params.variable_name)
            error('File %s does not contain variable %s.', ctx.chunk_files(k).fullpath, ctx.params.variable_name);
        end
        current_outputs = S.(ctx.params.variable_name);
    end

    local_validate_chunk_against_reference(ctx.ref_info, current_outputs, ctx.chunk_files(k).name);
    if k == 1
        local_debug_log(ctx.params, 'after_validate_first_chunk_against_reference');
    end

    phi_chunk = current_outputs.efuns(:, ctx.mode_info.selected_idx);
    if size(phi_chunk, 2) ~= ctx.mode_info.n_modes
        error('Chunk %s does not provide the expected selected mode count.', ctx.chunk_files(k).name);
    end

    u_chunk = compute_blp_chunk_residual( ...
        phi_chunk, ctx.mode_info.lambda_d_row, prev_phi_last, ctx.params.residual.first_u_mode, cursor == 1);
    if k == 1
        local_debug_log(ctx.params, sprintf('after_compute_first_chunk_residual size=%d x %d', size(u_chunk, 1), size(u_chunk, 2)));
    end

    session_state = local_accumulate_chunk(session_state, ctx.pair_specs, u_chunk, cursor);
    if k == 1
        local_debug_log(ctx.params, 'after_accumulate_first_chunk');
    end

    prev_phi_last = phi_chunk(end, :);
    cursor = cursor + size(phi_chunk, 1);
    if cursor > last_required_global_idx
        break;
    end
end

total_streamed_length = cursor - 1;
k_used = k;
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


function session_state = local_initialize_session_state(cfg, meta, Blp, session_mask, pair_specs)
session_ids = meta.session_ids(session_mask);
session_lengths = meta.session_lengths(session_mask);
session_dx = meta.session_dx(session_mask);
session_start_idx = meta.session_start_idx(session_mask);
session_end_idx = meta.session_end_idx(session_mask);

blp_lengths = double(Blp.session_lengths(session_mask));
blp_dx = double(Blp.session_dx(session_mask));
blp_data_cells = Blp.data_cells(session_mask);

n_pairs = numel(pair_specs);
n_channels = numel(Blp.selected_channels);
n_sessions = numel(session_ids);
pair_template = local_build_pair_stats_template(n_pairs, n_channels, 0);

session_state = repmat(struct( ...
    'session_id', [], ...
    'session_note', '', ...
    'lfp_start_idx', [], ...
    'lfp_end_idx', [], ...
    'lfp_length', [], ...
    'lfp_dx', [], ...
    'n_blp_samples_raw', [], ...
    'n_blp_samples_used', [], ...
    'lfp_duration_sec', [], ...
    'blp_duration_sec', [], ...
    'common_duration_sec', [], ...
    'blp_data', [], ...
    'pair_stats', pair_template), n_sessions, 1);

for i = 1:n_sessions
    n_lfp = session_lengths(i);
    n_blp_raw = blp_lengths(i);
    dx_lfp = session_dx(i);
    dx_blp = blp_dx(i);

    if abs(dx_lfp - dx_blp) > 1e-12
        error('Session %d has mismatched lfp dx (%.12g) and blp dx (%.12g).', ...
            session_ids(i), dx_lfp, dx_blp);
    end

    n_used = min(n_lfp, n_blp_raw);
    if n_used < 1
        error('Session %d has zero overlapping BLP samples.', session_ids(i));
    end

    session_state(i).session_id = session_ids(i);
    session_state(i).session_note = local_lookup_session_note(cfg, session_ids(i));
    session_state(i).lfp_start_idx = session_start_idx(i);
    session_state(i).lfp_end_idx = session_end_idx(i);
    session_state(i).lfp_length = n_lfp;
    session_state(i).lfp_dx = dx_lfp;
    session_state(i).n_blp_samples_raw = n_blp_raw;
    session_state(i).n_blp_samples_used = n_used;
    session_state(i).lfp_duration_sec = n_lfp * dx_lfp;
    session_state(i).blp_duration_sec = n_blp_raw * dx_blp;
    session_state(i).common_duration_sec = n_used * dx_lfp;
    session_state(i).blp_data = single(blp_data_cells{i}(1:n_used, :));
    session_state(i).pair_stats = local_build_pair_stats_template(n_pairs, n_channels, 0);
end
end


function pair_stats = local_build_pair_stats_template(n_pairs, n_channels, n_modes)
pair_stats = repmat(struct( ...
    'n_valid', 0, ...
    'sum_x', zeros(1, n_channels), ...
    'sum_x2', zeros(1, n_channels), ...
    'sum_y', zeros(1, n_modes), ...
    'sum_y2', zeros(1, n_modes), ...
    'sum_xy', zeros(n_channels, n_modes)), n_pairs, 1);
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


function session_state = local_accumulate_chunk(session_state, pair_specs, u_chunk, global_start_idx)
global_end_idx = global_start_idx + size(u_chunk, 1) - 1;
n_modes = size(u_chunk, 2);

for i = 1:numel(session_state)
    seg_start = max(global_start_idx, session_state(i).lfp_start_idx);
    seg_end = min(global_end_idx, session_state(i).lfp_end_idx);
    if seg_start > seg_end
        continue;
    end

    chunk_rows = (seg_start:seg_end) - global_start_idx + 1;
    local_rows = (seg_start:seg_end) - session_state(i).lfp_start_idx + 1;
    keep = local_rows >= 1 & local_rows <= session_state(i).n_blp_samples_used;
    if ~any(keep)
        continue;
    end

    chunk_rows = chunk_rows(keep);
    local_rows = local_rows(keep);

    Xraw = double(session_state(i).blp_data(local_rows, :));
    Useg = u_chunk(chunk_rows, :);
    Yabs = double(abs(Useg));
    Yreal = double(real(Useg));
    Yimag = double(imag(Useg));

    for p = 1:numel(pair_specs)
        stats = session_state(i).pair_stats(p);
        if isempty(stats.sum_y)
            stats.sum_y = zeros(1, n_modes);
            stats.sum_y2 = zeros(1, n_modes);
            stats.sum_xy = zeros(size(Xraw, 2), n_modes);
        elseif numel(stats.sum_y) ~= n_modes
            stats.sum_y = zeros(1, n_modes);
            stats.sum_y2 = zeros(1, n_modes);
            stats.sum_xy = zeros(size(Xraw, 2), n_modes);
        end

        switch pair_specs(p).signal_feature
            case 'raw'
                X = Xraw;
            case 'abs'
                X = abs(Xraw);
            otherwise
                error('Unsupported signal feature %s.', pair_specs(p).signal_feature);
        end

        switch pair_specs(p).residual_feature
            case 'abs'
                Y = Yabs;
            case 'real'
                Y = Yreal;
            case 'imag'
                Y = Yimag;
            otherwise
                error('Unsupported residual feature %s.', pair_specs(p).residual_feature);
        end

        valid = all(isfinite(X), 2) & all(isfinite(Y), 2);
        if ~any(valid)
            session_state(i).pair_stats(p) = stats;
            continue;
        end

        Xv = X(valid, :);
        Yv = Y(valid, :);
        nv = size(Xv, 1);

        stats.n_valid = stats.n_valid + nv;
        stats.sum_x = stats.sum_x + sum(Xv, 1);
        stats.sum_x2 = stats.sum_x2 + sum(Xv .^ 2, 1);
        stats.sum_y = stats.sum_y + sum(Yv, 1);
        stats.sum_y2 = stats.sum_y2 + sum(Yv .^ 2, 1);
        stats.sum_xy = stats.sum_xy + (Xv.' * Yv);

        session_state(i).pair_stats(p) = stats;
    end
end
end


function tf = local_should_print_progress(k, n_total, progress_every)
tf = (k == 1) || (k == n_total) || (mod(k, progress_every) == 0);
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
