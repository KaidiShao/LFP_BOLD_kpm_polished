% Compare global-vs-session z-score preprocessing for E10gb1 event detection.
%
% Scope of this test:
%   - same current raw-data root as cfg_E10gb1()
%   - same bandpass / threshold detector settings
%   - compare only the z-score scope:
%       1) global z-score per channel on the concatenated trace
%       2) session-wise z-score per channel before concatenation
%
% Outputs:
%   tests/test24_detection_alignment_checks/outputs_e10gb1_zscore_scope/*.csv
%   tests/test24_detection_alignment_checks/outputs_e10gb1_zscore_scope/*.mat

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(fileparts(this_script_dir));
addpath(genpath(repo_root));

close all force;

cfg = cfg_E10gb1();
D = io_raw.load_blp_dataset(cfg);

output_dir = fullfile(this_script_dir, 'outputs_e10gb1_zscore_scope');
if exist(output_dir, 'dir') ~= 7
    mkdir(output_dir);
end

params = struct();
params.passband = [2, 15; 30, 90; 90, 190];
params.band_labels = {'theta', 'gamma', 'ripple'};
params.L_start_range = [151, 101, 51];
params.L_extract_range = [301, 201, 101];
params.ThresRatio_range = [3.5, 4, 4];
params.max_rel_dx_diff = 1e-3;

[Fs, dx_ref] = local_resolve_fs(D, params);
state_var_raw = transpose(double(D.data));  % channel x time

[state_var_global, global_norm_info] = local_apply_global_zscore(state_var_raw);
[state_var_session, session_norm_info] = local_apply_session_zscore( ...
    state_var_raw, D.session_start_idx, D.session_end_idx, D.session_ids);

R_global = local_detect_from_preprocessed_state(D, cfg, state_var_global, Fs, dx_ref, params, 'global_zscore');
R_session = local_detect_from_preprocessed_state(D, cfg, state_var_session, Fs, dx_ref, params, 'sessionwise_zscore');

norm_summary_table = local_build_normalization_summary_table( ...
    cfg, D, state_var_raw, state_var_global, state_var_session);
count_table = local_build_variant_count_table({R_global, R_session}, params.band_labels);
peak_match_table = local_build_pairwise_peak_match_table(R_global, R_session, params.band_labels);
aggregate_table = local_build_aggregate_table(peak_match_table);
threshold_table = local_build_threshold_table(R_global, R_session, params.band_labels);

save_file = fullfile(output_dir, 'e10gb1_zscore_scope_comparison.mat');
save(save_file, ...
    'cfg', ...
    'params', ...
    'Fs', ...
    'dx_ref', ...
    'global_norm_info', ...
    'session_norm_info', ...
    'R_global', ...
    'R_session', ...
    'norm_summary_table', ...
    'count_table', ...
    'peak_match_table', ...
    'aggregate_table', ...
    'threshold_table', ...
    '-v7.3');

writetable(norm_summary_table, fullfile(output_dir, 'normalization_summary_table.csv'));
writetable(count_table, fullfile(output_dir, 'count_table.csv'));
writetable(peak_match_table, fullfile(output_dir, 'peak_match_table.csv'));
writetable(aggregate_table, fullfile(output_dir, 'aggregate_table.csv'));
writetable(threshold_table, fullfile(output_dir, 'threshold_table.csv'));

disp('Aggregate comparison between global-zscore and sessionwise-zscore:');
disp(aggregate_table);

fprintf('Saved z-score scope comparison outputs to:\n  %s\n', output_dir);


function [Fs, dx_ref] = local_resolve_fs(D, params)
if isfield(D, 'dx') && ~isempty(D.dx)
    dx_ref = double(D.dx);
else
    session_dx = double(D.session_dx(:));
    dx_ref = median(session_dx);
    max_rel_dx_diff = max(abs(session_dx - dx_ref)) / max(abs(dx_ref), eps);
    if max_rel_dx_diff > params.max_rel_dx_diff
        error(['Sampling period varies too much across sessions for a single detection run. ' ...
               'max relative dx difference = %.6g, allowed = %.6g.'], ...
               max_rel_dx_diff, params.max_rel_dx_diff);
    end
end
Fs = 1 / dx_ref;
end


function [state_var_out, info] = local_apply_global_zscore(state_var_in)
channel_mean = mean(state_var_in, 2, 'omitnan');
channel_std = std(state_var_in, 0, 2, 'omitnan');
safe_std = channel_std;
safe_std(~isfinite(safe_std) | safe_std <= 0) = 1;

state_var_out = state_var_in - channel_mean;
state_var_out = state_var_out ./ safe_std;

info = struct();
info.mode = 'global_zscore';
info.channel_mean = channel_mean;
info.channel_std = safe_std;
end


function [state_var_out, info] = local_apply_session_zscore(state_var_in, session_start_idx, session_end_idx, session_ids)
state_var_out = zeros(size(state_var_in), 'like', double(state_var_in));
n_channels = size(state_var_in, 1);
n_sessions = numel(session_start_idx);
session_mean = zeros(n_sessions, n_channels);
session_std = zeros(n_sessions, n_channels);

for k = 1:n_sessions
    idx1 = session_start_idx(k);
    idx2 = session_end_idx(k);
    block = double(state_var_in(:, idx1:idx2));

    mu = mean(block, 2, 'omitnan');
    sigma = std(block, 0, 2, 'omitnan');
    sigma(~isfinite(sigma) | sigma <= 0) = 1;

    block = block - mu;
    block = block ./ sigma;

    state_var_out(:, idx1:idx2) = block;
    session_mean(k, :) = transpose(mu);
    session_std(k, :) = transpose(sigma);
end

info = struct();
info.mode = 'sessionwise_zscore';
info.session_ids = double(session_ids(:));
info.session_mean = session_mean;
info.session_std = session_std;
end


function R = local_detect_from_preprocessed_state(D, cfg, state_var_detector, Fs, dx_ref, params, variant_name)
n_time = size(state_var_detector, 2);
n_channels = size(state_var_detector, 1);
n_bands = size(params.passband, 1);
n_sessions = numel(D.session_ids);

channel_sites = cfg.channels.sites(D.selected_channels);
regions = cfg.channels.selected_labels;
channel_region_idx = zeros(1, n_channels);
for c = 1:n_channels
    channel_region_idx(c) = find(strcmp(regions, channel_sites{c}), 1, 'first');
end

DetectResults = cell(n_bands, n_channels);

for n_band = 1:n_bands
    filtered_band = filterSignal_mirror( ...
        transpose(state_var_detector), Fs, params.passband(n_band, 1), params.passband(n_band, 2), 1);
    filtered_band = transpose(double(filtered_band));

    for n_channel = 1:n_channels
        trace_all = filtered_band(n_channel, :);
        threshold = mean(trace_all, 'omitnan') + ...
            params.ThresRatio_range(n_band) * std(trace_all, 0, 'omitnan');

        loc_pooled = zeros(0, 1);
        event_win = zeros(0, 2);
        event_win_session_local = zeros(0, 2);
        event_session_idx = zeros(0, 1);
        event_session_id = zeros(0, 1);
        loc_peak = zeros(0, 1);

        for k = 1:n_sessions
            idx1 = D.session_start_idx(k);
            idx2 = D.session_end_idx(k);

            trace_session = trace_all(idx1:idx2);
            loc_session = find(trace_session > threshold);
            if ~isempty(loc_session)
                loc_session(1) = [];
            end
            if isempty(loc_session)
                continue;
            end

            loc_pooled = [loc_pooled; idx1 + loc_session(:) - 1]; %#ok<AGROW>
            loc_peak_local = find_peak_loc(trace_session, loc_session, params.L_extract_range(n_band));

            for i_peak = 1:numel(loc_peak_local)
                peak_local = loc_peak_local(i_peak);
                win_start_local = peak_local - params.L_start_range(n_band) + 1;
                win_end_local = win_start_local + params.L_extract_range(n_band) - 1;

                if win_start_local < 1 || win_end_local > numel(trace_session)
                    continue;
                end

                loc_peak(end+1, 1) = idx1 + peak_local - 1; %#ok<AGROW>
                event_win(end+1, :) = [idx1 + win_start_local - 1, idx1 + win_end_local - 1]; %#ok<AGROW>
                event_win_session_local(end+1, :) = [win_start_local, win_end_local]; %#ok<AGROW>
                event_session_idx(end+1, 1) = k; %#ok<AGROW>
                event_session_id(end+1, 1) = D.session_ids(k); %#ok<AGROW>
            end
        end

        out = struct();
        out.loc_pooled = loc_pooled;
        out.loc_peak = loc_peak;
        out.event_win = event_win;
        out.event_win_session_local = event_win_session_local;
        out.event_session_idx = event_session_idx;
        out.event_session_id = event_session_id;
        out.L_total = n_time;
        out.dx = dx_ref;
        out.Fs = Fs;
        out.threshold = threshold;
        out.bandpass_hz = params.passband(n_band, :);
        out.band_index = n_band;
        out.band_label = params.band_labels{n_band};
        out.channel_col = n_channel;
        out.selected_channel = D.selected_channels(n_channel);
        out.channel_site = channel_sites{n_channel};
        out.region_idx = channel_region_idx(n_channel);
        out.region_label = regions{channel_region_idx(n_channel)};
        out.variant_name = char(string(variant_name));
        DetectResults{n_band, n_channel} = out;
    end
end

R = struct();
R.variant_name = char(string(variant_name));
R.DetectResults = DetectResults;
R.params = params;
R.session_ids = D.session_ids;
R.session_lengths = D.session_lengths;
R.session_dx = D.session_dx;
R.session_start_idx = D.session_start_idx;
R.session_end_idx = D.session_end_idx;
R.selected_channels = D.selected_channels;
end


function tbl = local_build_normalization_summary_table(cfg, D, raw_state, global_state, session_state)
n_channels = size(raw_state, 1);
n_sessions = numel(D.session_ids);
rows = cell(n_sessions * n_channels, 1);
idx = 0;

for s = 1:n_sessions
    idx1 = D.session_start_idx(s);
    idx2 = D.session_end_idx(s);
    for c = 1:n_channels
        idx = idx + 1;
        x_raw = double(raw_state(c, idx1:idx2));
        x_global = double(global_state(c, idx1:idx2));
        x_session = double(session_state(c, idx1:idx2));

        rows{idx} = table( ...
            D.session_ids(s), ...
            c, ...
            string(cfg.channels.sites{D.selected_channels(c)}), ...
            mean(x_raw), ...
            std(x_raw, 0), ...
            mean(x_global), ...
            std(x_global, 0), ...
            mean(x_session), ...
            std(x_session, 0), ...
            'VariableNames', { ...
                'session_id', ...
                'channel_col', ...
                'channel_site', ...
                'raw_mean', ...
                'raw_std', ...
                'global_zscore_mean', ...
                'global_zscore_std', ...
                'sessionwise_zscore_mean', ...
                'sessionwise_zscore_std'});
    end
end

tbl = vertcat(rows{:});
end


function tbl = local_build_variant_count_table(results, band_labels)
rows = cell(0, 1);
for r = 1:numel(results)
    DetectResults = results{r}.DetectResults;
    [n_bands, n_channels] = size(DetectResults);
    for b = 1:n_bands
        variant_name = repmat(string(results{r}.variant_name), n_channels, 1);
        band_label = repmat(string(band_labels{b}), n_channels, 1);
        channel_idx = transpose(1:n_channels);
        event_count = zeros(n_channels, 1);
        threshold = zeros(n_channels, 1);
        for c = 1:n_channels
            event_count(c) = size(DetectResults{b, c}.event_win, 1);
            threshold(c) = DetectResults{b, c}.threshold;
        end
        rows{end+1, 1} = table(variant_name, band_label, channel_idx, event_count, threshold); %#ok<AGROW>
    end
end
tbl = vertcat(rows{:});
end


function tbl = local_build_pairwise_peak_match_table(R_global, R_session, band_labels)
[n_bands, n_channels] = size(R_global.DetectResults);
rows = cell(n_bands * n_channels, 1);
idx = 0;

for b = 1:n_bands
    for c = 1:n_channels
        idx = idx + 1;
        global_peaks = double(R_global.DetectResults{b, c}.loc_peak(:));
        session_peaks = double(R_session.DetectResults{b, c}.loc_peak(:));
        shared_peaks = intersect(global_peaks, session_peaks);
        union_peaks = union(global_peaks, session_peaks);

        if isempty(union_peaks)
            peak_jaccard = 1;
        else
            peak_jaccard = numel(shared_peaks) / numel(union_peaks);
        end

        rows{idx} = table( ...
            string(band_labels{b}), ...
            c, ...
            numel(global_peaks), ...
            numel(session_peaks), ...
            numel(shared_peaks), ...
            numel(setdiff(global_peaks, session_peaks)), ...
            numel(setdiff(session_peaks, global_peaks)), ...
            peak_jaccard, ...
            'VariableNames', { ...
                'band_label', ...
                'channel_idx', ...
                'global_peak_count', ...
                'sessionwise_peak_count', ...
                'shared_peak_count', ...
                'global_only_peak_count', ...
                'sessionwise_only_peak_count', ...
                'peak_jaccard'});
    end
end

tbl = vertcat(rows{:});
end


function tbl = local_build_aggregate_table(peak_match_table)
global_peak_count_total = sum(peak_match_table.global_peak_count);
sessionwise_peak_count_total = sum(peak_match_table.sessionwise_peak_count);
shared_peak_count_total = sum(peak_match_table.shared_peak_count);

tbl = table( ...
    global_peak_count_total, ...
    sessionwise_peak_count_total, ...
    shared_peak_count_total, ...
    shared_peak_count_total / max(global_peak_count_total, 1), ...
    shared_peak_count_total / max(sessionwise_peak_count_total, 1), ...
    mean(peak_match_table.peak_jaccard), ...
    min(peak_match_table.peak_jaccard), ...
    max(peak_match_table.peak_jaccard), ...
    'VariableNames', { ...
        'global_peak_count_total', ...
        'sessionwise_peak_count_total', ...
        'shared_peak_count_total', ...
        'recall_vs_global', ...
        'precision_vs_global', ...
        'mean_peak_jaccard', ...
        'min_peak_jaccard', ...
        'max_peak_jaccard'});
end


function tbl = local_build_threshold_table(R_global, R_session, band_labels)
[n_bands, n_channels] = size(R_global.DetectResults);
rows = cell(n_bands * n_channels, 1);
idx = 0;

for b = 1:n_bands
    for c = 1:n_channels
        idx = idx + 1;
        t_global = double(R_global.DetectResults{b, c}.threshold);
        t_session = double(R_session.DetectResults{b, c}.threshold);
        rows{idx} = table( ...
            string(band_labels{b}), ...
            c, ...
            t_global, ...
            t_session, ...
            t_session - t_global, ...
            t_session / max(abs(t_global), eps), ...
            'VariableNames', { ...
                'band_label', ...
                'channel_idx', ...
                'global_threshold', ...
                'sessionwise_threshold', ...
                'threshold_delta', ...
                'threshold_ratio_vs_global'});
    end
end

tbl = vertcat(rows{:});
end
