function [session_state, total_streamed_length, k_used] = align_spkt_residual_to_spike_resolution(ctx)
%ALIGN_SPKT_RESIDUAL_TO_SPIKE_RESOLUTION Stream residuals onto SPKT bins.

session_state = local_initialize_session_state(ctx.cfg, ctx.meta, ctx.Spk, ctx.session_mask, ctx.mode_info, ctx.params);
last_required_global_idx = max([session_state.lfp_end_idx]);
n_chunks = numel(ctx.chunk_files);
k_used = 0;

if strcmp(ctx.source_mode, 'residual_workspace')
    u_full = ctx.residual_bundle.u_full(:, 1:ctx.mode_info.n_modes);
    session_state = local_accumulate_chunk(session_state, u_full, 1);
    total_streamed_length = size(u_full, 1);
    k_used = n_chunks;
    return;
end

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

    phi_chunk = current_outputs.efuns(:, ctx.mode_info.selected_idx);
    if size(phi_chunk, 2) ~= ctx.mode_info.n_modes
        error('Chunk %s does not provide the expected selected mode count.', ctx.chunk_files(k).name);
    end

    u_chunk = compute_blp_chunk_residual( ...
        phi_chunk, ctx.mode_info.lambda_d_row, prev_phi_last, ctx.params.residual.first_u_mode, cursor == 1);

    session_state = local_accumulate_chunk(session_state, u_chunk, cursor);

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

    abs_u = double(abs(u_seg));
    session_state(i).sum_abs = session_state(i).sum_abs + single(B * abs_u);
    session_state(i).sum_abs_sq = session_state(i).sum_abs_sq + single(B * (abs_u .^ 2));
    session_state(i).sum_real = session_state(i).sum_real + single(B * double(real(u_seg)));
    session_state(i).sum_imag = session_state(i).sum_imag + single(B * double(imag(u_seg)));
end
end


function tf = local_should_print_progress(k, n_total, progress_every)
tf = (k == 1) || (k == n_total) || (mod(k, progress_every) == 0);
end
