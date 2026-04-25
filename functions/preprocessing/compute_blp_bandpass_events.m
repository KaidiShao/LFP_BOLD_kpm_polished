function R = compute_blp_bandpass_events(D, cfg, output_root, params)
% Detect bandpass-threshold events from loaded BLP data.
%
% Inputs
%   D           Output struct from load_blp_dataset
%   cfg         Dataset config struct
%   output_root Root folder for processed outputs
%   params      Struct with detection parameters
%
% Output
%   R           Struct containing detection results and metadata
%
% Notes
%   - Detection is performed channel-wise on the selected BLP channels.
%   - Thresholds are computed on the full concatenated trace for each
%     band/channel pair, matching the legacy behavior.
%   - By default the detector first z-scores each selected channel on the
%     full concatenated trace, which better matches the older normalized
%     BLP roots used in test23-style analyses.
%   - Filtering and peak grouping use filterSignal_mirror and
%     find_peak_loc from utils to stay aligned with the legacy pipeline.
%   - Peak picking and event extraction are done session by session so
%     event windows never cross session borders.

if nargin < 3 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 4
    params = struct();
end

params = apply_default_params(params);

if isempty(D.data)
    error('D.data is empty.');
end

if exist('filterSignal_mirror', 'file') ~= 2
    error('filterSignal_mirror.m was not found on the MATLAB path.');
end

if exist('find_peak_loc', 'file') ~= 2
    error('find_peak_loc.m was not found on the MATLAB path.');
end

if isfield(D, 'dx') && ~isempty(D.dx)
    dx_ref = D.dx;
    dx_source = 'D.dx';
else
    if ~isfield(D, 'session_dx') || isempty(D.session_dx)
        error('Neither D.dx nor D.session_dx is available for event detection.');
    end

    session_dx = double(D.session_dx(:));
    dx_ref = median(session_dx);
    max_rel_dx_diff = max(abs(session_dx - dx_ref)) / max(abs(dx_ref), eps);

    if max_rel_dx_diff > params.max_rel_dx_diff
        error(['Sampling period varies too much across sessions for a single detection run. ' ...
               'max relative dx difference = %.6g, allowed = %.6g.'], ...
               max_rel_dx_diff, params.max_rel_dx_diff);
    end

    dx_source = 'median_session_dx';
    warning(['D.dx is empty because session dx values are not exactly equal. ' ...
             'Using median(session_dx) = %.12g for filtering. Max relative dx difference = %.6g.'], ...
             dx_ref, max_rel_dx_diff);
end

Fs = 1 / dx_ref;
n_time = size(D.data, 1);
n_channels = size(D.data, 2);
n_bands = size(params.passband, 1);
n_sessions = numel(D.session_ids);

channel_sites = cfg.channels.sites(D.selected_channels);
if isfield(cfg.channels, 'selected_labels') && ~isempty(cfg.channels.selected_labels)
    regions = cfg.channels.selected_labels;
else
    regions = unique(channel_sites, 'stable');
end

channel_region_idx = zeros(1, n_channels);
for c = 1:n_channels
    ridx = find(strcmp(regions, channel_sites{c}), 1, 'first');
    if isempty(ridx)
        error('Selected channel site "%s" is not present in region labels.', channel_sites{c});
    end
    channel_region_idx(c) = ridx;
end

save_dir = fullfile(output_root, cfg.file_stem, 'event_detection');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_file = fullfile(save_dir, [build_save_tag(cfg.file_stem, params), '.mat']);

if exist(save_file, 'file') == 2 && ~params.force_recompute
    S = load(save_file, 'R');
    if isfield(S, 'R') && is_saved_result_compatible(S.R, params)
        R = S.R;
        return;
    end
    warning(['Existing event-result file does not match the requested detector parameters ' ...
             'or is missing new metadata. Recomputing:\n  %s'], save_file);
end

state_var_raw = transpose(double(D.data));  % channel x time
[state_var_detector, normalization_info] = apply_input_normalization(state_var_raw, params);

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

        event_data = {};
        event_data_detector_input = {};
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

                peak_global = idx1 + peak_local - 1;
                win_start_global = idx1 + win_start_local - 1;
                win_end_global = idx1 + win_end_local - 1;

                loc_peak(end+1, 1) = peak_global;
                event_win(end+1, :) = [win_start_global, win_end_global];
                event_win_session_local(end+1, :) = [win_start_local, win_end_local];
                event_session_idx(end+1, 1) = k;
                event_session_id(end+1, 1) = D.session_ids(k);
                event_data{end+1, 1} = state_var_raw(n_channel, win_start_global:win_end_global); %#ok<AGROW>
                event_data_detector_input{end+1, 1} = ...
                    state_var_detector(n_channel, win_start_global:win_end_global); %#ok<AGROW>
            end
        end

        out = struct();
        out.loc_pooled = loc_pooled;
        out.loc_peak = loc_peak;
        out.event_data = event_data;
        out.event_data_raw = event_data;
        out.event_data_detector_input = event_data_detector_input;
        out.event_win = event_win;
        out.event_win_session_local = event_win_session_local;
        out.event_session_idx = event_session_idx;
        out.event_session_id = event_session_id;
        out.L_total = n_time;
        out.dx = dx_ref;
        out.dx_source = dx_source;
        out.Fs = Fs;
        out.threshold = threshold;
        out.bandpass_hz = params.passband(n_band, :);
        out.input_normalization = normalization_info.mode;
        out.input_channel_mean = normalization_info.channel_mean(n_channel);
        out.input_channel_std = normalization_info.channel_std(n_channel);
        out.band_index = n_band;
        out.band_label = params.band_labels{n_band};
        out.channel_col = n_channel;
        out.selected_channel = D.selected_channels(n_channel);
        out.channel_site = channel_sites{n_channel};
        out.region_idx = channel_region_idx(n_channel);
        out.region_label = regions{channel_region_idx(n_channel)};
        out.session_ids = D.session_ids;
        out.session_dx = D.session_dx;
        out.session_start_idx = D.session_start_idx;
        out.session_end_idx = D.session_end_idx;
        out.params = params;

        DetectResults{n_band, n_channel} = out;
    end
end

R = struct();
R.DetectResults = DetectResults;
R.params = params;
R.save_file = save_file;
R.dataset_id = cfg.dataset_id;
R.file_stem = cfg.file_stem;
R.selected_channels = D.selected_channels;
R.channel_sites = channel_sites;
R.channel_region_idx = channel_region_idx;
R.regions = regions;
R.session_ids = D.session_ids;
R.session_lengths = D.session_lengths;
R.session_dx = D.session_dx;
R.session_start_idx = D.session_start_idx;
R.session_end_idx = D.session_end_idx;
R.border_idx = D.border_idx;
R.dx = dx_ref;
R.dx_source = dx_source;
R.Fs = Fs;
R.input_normalization = normalization_info.mode;
R.input_channel_mean = normalization_info.channel_mean;
R.input_channel_std = normalization_info.channel_std;
R.input_normalization_notes = normalization_info.notes;

save(save_file, 'R', '-v7.3');
end

function params = apply_default_params(params)
if ~isfield(params, 'passband')
    params.passband = [2, 15; 30, 90; 90, 190];
end

if ~isfield(params, 'band_labels') || isempty(params.band_labels)
    params.band_labels = {'band01', 'band02', 'band03'};
end

if ~isfield(params, 'L_start_range')
    params.L_start_range = [151, 101, 51];
end

if ~isfield(params, 'L_extract_range')
    params.L_extract_range = [301, 201, 101];
end

if ~isfield(params, 'ThresRatio_range')
    params.ThresRatio_range = [3.5, 4, 4];
end

if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end

if ~isfield(params, 'max_rel_dx_diff')
    params.max_rel_dx_diff = 1e-3;
end

if ~isfield(params, 'input_normalization') || isempty(params.input_normalization)
    params.input_normalization = 'zscore_per_channel';
end

n_bands = size(params.passband, 1);

if numel(params.band_labels) ~= n_bands
    error('params.band_labels must have one entry per passband.');
end
if numel(params.L_start_range) ~= n_bands
    error('params.L_start_range must match the number of passbands.');
end
if numel(params.L_extract_range) ~= n_bands
    error('params.L_extract_range must match the number of passbands.');
end
if numel(params.ThresRatio_range) ~= n_bands
    error('params.ThresRatio_range must match the number of passbands.');
end
end

function [state_var_out, info] = apply_input_normalization(state_var_in, params)
mode = char(string(params.input_normalization));
state_var_out = double(state_var_in);

info = struct();
info.mode = mode;
info.channel_mean = mean(state_var_out, 2, 'omitnan');
info.channel_std = std(state_var_out, 0, 2, 'omitnan');
info.notes = '';

switch lower(mode)
    case {'none', 'raw'}
        info.notes = 'Detector input uses the raw selected BLP channels without normalization.';

    case {'zscore', 'zscore_per_channel', 'per_channel_zscore'}
        zero_std_mask = ~isfinite(info.channel_std) | info.channel_std <= 0;
        safe_std = info.channel_std;
        safe_std(zero_std_mask) = 1;

        state_var_out = state_var_out - info.channel_mean;
        state_var_out = state_var_out ./ safe_std;
        info.channel_std = safe_std;
        if any(zero_std_mask)
            info.notes = ['Detector input z-scored each channel on the full concatenated trace. ' ...
                          'Channels with non-positive std were left scaled by 1.'];
        else
            info.notes = 'Detector input z-scored each channel on the full concatenated trace.';
        end

    otherwise
        error('Unsupported params.input_normalization: %s', mode);
end
end

function tf = is_saved_result_compatible(R_saved, params)
tf = false;

if ~isstruct(R_saved) || ~isfield(R_saved, 'params') || ~isstruct(R_saved.params)
    return;
end

saved_params = R_saved.params;
fields_to_compare = { ...
    'passband', ...
    'band_labels', ...
    'L_start_range', ...
    'L_extract_range', ...
    'ThresRatio_range', ...
    'max_rel_dx_diff', ...
    'input_normalization'};

for i = 1:numel(fields_to_compare)
    f = fields_to_compare{i};
    if ~isfield(saved_params, f) || ~isfield(params, f)
        return;
    end
    if ~isequaln(saved_params.(f), params.(f))
        return;
    end
end

if ~isfield(R_saved, 'input_normalization')
    return;
end

tf = true;
end

function tag = build_save_tag(file_stem, params)
n_bands = size(params.passband, 1);
tag = sprintf('%s_bandpass_events_%dbands', file_stem, n_bands);
end
