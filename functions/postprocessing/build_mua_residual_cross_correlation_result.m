function result = build_mua_residual_cross_correlation_result(cfg, params, ctx, session_state, total_streamed_length, k_used)
%BUILD_MUA_RESIDUAL_CROSS_CORRELATION_RESULT Build MUA xcorr result struct.

[result_sessions, session_summary, pooled_pair_stats] = local_finalize_sessions(session_state, ctx.pair_specs);
local_debug_log(params, 'after_finalize_sessions');
mode_table = local_build_mode_table(ctx.mode_info);
channel_table = local_build_channel_table(ctx.Blp.selected_channels, ctx.Blp.channel_sites);
session_corr_table = local_build_session_corr_table(result_sessions, ctx.mode_info, channel_table, ctx.pair_specs);
pooled_corr_table = local_build_pooled_corr_table(pooled_pair_stats, ctx.mode_info, channel_table, ctx.pair_specs);
local_debug_log(params, 'after_build_corr_tables');

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

if isempty(pooled_corr_table)
    pooled_top_corr = table();
else
    finite_top = pooled_corr_table(isfinite(pooled_corr_table.abs_corr), :);
    pooled_top_corr = finite_top(1:min(height(finite_top), params.top_n_rows), :);
end

result = struct();
result.cfg = cfg;
result.params = params;
result.pair_specs = ctx.pair_specs;
result.source = struct();
result.source.mode = ctx.source_mode;
result.source.data_dir = ctx.source_dir;
result.source.selection = ctx.source_select_info;
result.source.blp_channel_selection_label = ctx.Blp.channel_selection_label;
result.source.blp_band_selection_label = ctx.Blp.band_selection_label;
result.source.blp_band_index = ctx.Blp.band_index;
result.source.n_blp_bands = ctx.Blp.n_bands;
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
result.blp_selected_channels = ctx.Blp.selected_channels;
result.blp_channel_sites = ctx.Blp.channel_sites;
result.total_streamed_length = total_streamed_length;
result.mode_table = mode_table;
result.channel_table = channel_table;
result.session_summary = session_summary;
result.session_corr_table = session_corr_table;
result.pooled_corr_table = pooled_corr_table;
result.pooled_top_corr = pooled_top_corr;
result.sessions = result_sessions;
end


function [result_sessions, session_summary, pooled_pair_stats] = local_finalize_sessions(session_state, pair_specs)
n_sessions = numel(session_state);
n_pairs = numel(pair_specs);

result_sessions = repmat(struct( ...
    'session_id', [], ...
    'session_note', '', ...
    'lfp_dx', [], ...
    'n_lfp_samples', [], ...
    'n_blp_samples_raw', [], ...
    'n_blp_samples_used', [], ...
    'lfp_duration_sec', [], ...
    'blp_duration_sec', [], ...
    'common_duration_sec', [], ...
    'corr', []), n_sessions, 1);

session_id_col = zeros(n_sessions, 1);
session_note_col = strings(n_sessions, 1);
lfp_dx_col = zeros(n_sessions, 1);
n_lfp_col = zeros(n_sessions, 1);
n_blp_raw_col = zeros(n_sessions, 1);
n_blp_used_col = zeros(n_sessions, 1);
lfp_dur_col = zeros(n_sessions, 1);
blp_dur_col = zeros(n_sessions, 1);
common_dur_col = zeros(n_sessions, 1);

pooled_pair_stats = [];

for i = 1:n_sessions
    corr_struct = struct();
    for p = 1:n_pairs
        pair_name = pair_specs(p).name;
        corr_struct.(pair_name) = local_corr_matrix_from_stats(session_state(i).pair_stats(p));
        if isempty(pooled_pair_stats)
            pooled_pair_stats = session_state(i).pair_stats;
        else
            pooled_pair_stats(p) = local_add_pair_stats(pooled_pair_stats(p), session_state(i).pair_stats(p));
        end
    end

    result_sessions(i).session_id = session_state(i).session_id;
    result_sessions(i).session_note = session_state(i).session_note;
    result_sessions(i).lfp_dx = session_state(i).lfp_dx;
    result_sessions(i).n_lfp_samples = session_state(i).lfp_length;
    result_sessions(i).n_blp_samples_raw = session_state(i).n_blp_samples_raw;
    result_sessions(i).n_blp_samples_used = session_state(i).n_blp_samples_used;
    result_sessions(i).lfp_duration_sec = session_state(i).lfp_duration_sec;
    result_sessions(i).blp_duration_sec = session_state(i).blp_duration_sec;
    result_sessions(i).common_duration_sec = session_state(i).common_duration_sec;
    result_sessions(i).corr = corr_struct;

    session_id_col(i) = session_state(i).session_id;
    session_note_col(i) = string(session_state(i).session_note);
    lfp_dx_col(i) = session_state(i).lfp_dx;
    n_lfp_col(i) = session_state(i).lfp_length;
    n_blp_raw_col(i) = session_state(i).n_blp_samples_raw;
    n_blp_used_col(i) = session_state(i).n_blp_samples_used;
    lfp_dur_col(i) = session_state(i).lfp_duration_sec;
    blp_dur_col(i) = session_state(i).blp_duration_sec;
    common_dur_col(i) = session_state(i).common_duration_sec;
end

if isempty(pooled_pair_stats)
    pooled_pair_stats = local_build_pair_stats_template(n_pairs, 0, 0);
end

session_summary = table( ...
    session_id_col, ...
    session_note_col, ...
    lfp_dx_col, ...
    n_lfp_col, ...
    n_blp_raw_col, ...
    n_blp_used_col, ...
    lfp_dur_col, ...
    blp_dur_col, ...
    common_dur_col, ...
    'VariableNames', { ...
        'session_id', 'session_note', 'lfp_dx_sec', ...
        'n_lfp_samples', 'n_blp_samples_raw', 'n_blp_samples_used', ...
        'lfp_duration_sec', 'blp_duration_sec', 'common_duration_sec'});
end


function out = local_add_pair_stats(a, b)
out = a;
out.n_valid = a.n_valid + b.n_valid;
out.sum_x = a.sum_x + b.sum_x;
out.sum_x2 = a.sum_x2 + b.sum_x2;
out.sum_y = a.sum_y + b.sum_y;
out.sum_y2 = a.sum_y2 + b.sum_y2;
out.sum_xy = a.sum_xy + b.sum_xy;
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


function session_corr_table = local_build_session_corr_table(result_sessions, mode_info, channel_table, pair_specs)
tables = cell(numel(result_sessions), 1);
for i = 1:numel(result_sessions)
    tables{i} = local_build_corr_table_from_struct( ...
        result_sessions(i).corr, ...
        result_sessions(i).session_id, ...
        result_sessions(i).session_note, ...
        result_sessions(i).n_blp_samples_used, ...
        mode_info, channel_table, pair_specs);
end

tables = tables(~cellfun(@isempty, tables));
if isempty(tables)
    session_corr_table = table();
else
    session_corr_table = vertcat(tables{:});
end
end


function pooled_corr_table = local_build_pooled_corr_table(pooled_pair_stats, mode_info, channel_table, pair_specs)
if isempty(pooled_pair_stats)
    pooled_corr_table = table();
    return;
end

corr_struct = struct();
n_valid_vec = zeros(numel(pair_specs), 1);
for p = 1:numel(pair_specs)
    corr_struct.(pair_specs(p).name) = local_corr_matrix_from_stats(pooled_pair_stats(p));
    n_valid_vec(p) = pooled_pair_stats(p).n_valid;
end

pooled_corr_table = local_build_corr_table_from_struct( ...
    corr_struct, NaN, 'pooled', max(n_valid_vec), mode_info, channel_table, pair_specs);
end


function T = local_build_corr_table_from_struct(corr_struct, session_id, session_note, n_samples, mode_info, channel_table, pair_specs)
tables = cell(numel(pair_specs), 1);

for p = 1:numel(pair_specs)
    pair_name = pair_specs(p).name;
    C = corr_struct.(pair_name);
    [n_channels, n_modes] = size(C);
    if n_channels == 0 || n_modes == 0
        tables{p} = table();
        continue;
    end

    [channel_grid, mode_grid] = ndgrid(1:n_channels, 1:n_modes);
    n_rows = numel(C);

    session_id_col = repmat(session_id, n_rows, 1);
    session_note_col = repmat(string(session_note), n_rows, 1);
    pairing_label_col = repmat(string(pair_specs(p).name), n_rows, 1);
    signal_feature_col = repmat(string(pair_specs(p).signal_feature), n_rows, 1);
    residual_feature_col = repmat(string(pair_specs(p).residual_feature), n_rows, 1);
    n_samples_col = repmat(n_samples, n_rows, 1);

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

    tables{p} = table( ...
        session_id_col, ...
        session_note_col, ...
        channel_index, ...
        channel_site, ...
        channel_label, ...
        mode_rank, ...
        mode_original_idx, ...
        pairing_label_col, ...
        signal_feature_col, ...
        residual_feature_col, ...
        evalue_real, ...
        evalue_imag, ...
        evalue_abs, ...
        corr_val, ...
        abs_corr, ...
        n_samples_col, ...
        'VariableNames', { ...
            'session_id', 'session_note', ...
            'channel_index', 'channel_site', 'channel_label', ...
            'mode_rank', 'mode_original_idx', ...
            'pairing_label', 'signal_feature', 'residual_feature', ...
            'evalue_real', 'evalue_imag', 'evalue_abs', ...
            'corr', 'abs_corr', 'n_samples'});
end

tables = tables(~cellfun(@isempty, tables));
if isempty(tables)
    T = table();
else
    T = vertcat(tables{:});
end
end


function C = local_corr_matrix_from_stats(stats)
if stats.n_valid < 2 || isempty(stats.sum_xy)
    C = nan(numel(stats.sum_x), numel(stats.sum_y));
    return;
end

n = double(stats.n_valid);
sum_x = double(stats.sum_x);
sum_x2 = double(stats.sum_x2);
sum_y = double(stats.sum_y);
sum_y2 = double(stats.sum_y2);
sum_xy = double(stats.sum_xy);

var_x_num = sum_x2 - (sum_x .^ 2) ./ n;
var_y_num = sum_y2 - (sum_y .^ 2) ./ n;
cov_num = sum_xy - (sum_x(:) * sum_y(:).') ./ n;

var_x = var_x_num ./ (n - 1);
var_y = var_y_num ./ (n - 1);
cov_xy = cov_num ./ (n - 1);

var_x(var_x <= 0) = NaN;
var_y(var_y <= 0) = NaN;

C = cov_xy ./ sqrt(var_x(:) * var_y(:).');
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
