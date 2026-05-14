function result = build_spkt_residual_cross_correlation_result(cfg, params, ctx, session_state, total_streamed_length, k_used)
%BUILD_SPKT_RESIDUAL_CROSS_CORRELATION_RESULT Build SPKT xcorr result struct.

[result_sessions, session_summary] = local_finalize_sessions(cfg, session_state);
mode_table = local_build_mode_table(ctx.mode_info);
channel_table = local_build_channel_table(ctx.Spk.selected_channels, ctx.Spk.channel_sites);
session_corr_table = local_build_session_corr_table(result_sessions, ctx.mode_info, channel_table);
pooled_corr_table = local_build_pooled_corr_table(result_sessions, ctx.mode_info, channel_table);

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
result.source.mode = ctx.source_mode;
result.source.data_dir = ctx.source_dir;
result.source.selection = ctx.source_select_info;
result.source.spike_channel_selection_label = ctx.Spk.channel_selection_label;
result.source.first_chunk = ctx.chunk_files(1).fullpath;
result.source.last_chunk = ctx.chunk_files(min(numel(ctx.chunk_files), max(1, k_used))).fullpath;
result.source.n_chunks_total = numel(ctx.chunk_files);
result.source.n_chunks_used = k_used;
result.source.filename_pattern = params.filename_pattern;
result.source.variable_name = params.variable_name;
result.source.ref_info = ctx.ref_info;
result.metadata_file = ctx.metadata_file;
result.lfp_selected_channels = ctx.meta.selected_channels;
result.lfp_channel_sites = ctx.meta.channel_sites;
result.spike_selected_channels = ctx.Spk.selected_channels;
result.spike_channel_sites = ctx.Spk.channel_sites;
result.total_streamed_length = total_streamed_length;
result.mode_table = mode_table;
result.channel_table = channel_table;
result.session_summary = session_summary;
result.session_corr_table = session_corr_table;
result.pooled_corr_table = pooled_corr_table;
result.pooled_top_corr = pooled_top_corr;
result.sessions = result_sessions;
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
