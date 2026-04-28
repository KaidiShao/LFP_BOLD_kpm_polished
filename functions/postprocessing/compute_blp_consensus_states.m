function C = compute_blp_consensus_states(cfg, output_root, R, params, source_event_file)
% Compute consensus state labels from saved bandpass event detections.
%
% Inputs
%   cfg               Dataset config struct
%   output_root       Root folder for processed outputs
%   R                 Loaded event-detection result struct
%   params            Struct with consensus-state parameters
%   source_event_file Optional source event-result file path for metadata
%
% Output
%   C                 Struct containing band-level consensus masks plus
%                     event-window-defined state labels and summaries

%% =========================
%  Input defaults
%  =========================
if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 3
    error('R must be provided. Load it first with load_event_results.');
end

if nargin < 4
    params = struct();
end

if nargin < 5
    source_event_file = '';
end

params = apply_default_params(params);
validate_params(params);

%% =========================
%  Validate event-detection results
%  =========================
if ~isfield(R, 'DetectResults') || isempty(R.DetectResults)
    error('The loaded event result does not contain DetectResults.');
end

if isempty(source_event_file) && isfield(R, 'save_file') && ~isempty(R.save_file)
    source_event_file = R.save_file;
end
meta = resolve_consensus_context(R, params);

%% =========================
%  Prepare save path
%  =========================
save_dir = fullfile(output_root, cfg.file_stem, 'consensus_states');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

save_tag = build_save_tag(cfg.file_stem, meta.min_channel_count, params.require_region_presence, params.required_regions);
save_file = fullfile(save_dir, [save_tag, '.mat']);
source_event_file_signature = io_utils.build_file_signature(source_event_file);

if exist(save_file, 'file') == 2 && ~params.force_recompute
    S = load(save_file);
    if isfield(S, 'C') && can_reuse_saved_result( ...
            S.C, params, meta.min_channel_count, meta.band_role_indices, source_event_file_signature)
        C = S.C;
        return;
    end
end

%% =========================
%  Band-level consensus masks
%  =========================
[consensus_mask_by_band, channel_support_count_by_band, ...
    region_support_count_by_band, band_consensus] = ...
    build_band_consensus_outputs(R.DetectResults, meta, params.require_region_presence);

%% =========================
%  Derived state labels
%  =========================
[state_catalog, state_windows, state_code_by_time, ...
    state_sample_count, state_duration_sec, state_window_count] = ...
    build_state_outputs( ...
        band_consensus, ...
        consensus_mask_by_band, ...
        channel_support_count_by_band, ...
        region_support_count_by_band, ...
        meta, ...
        params.state_merge_gap_samples);

%% =========================
%  Pack output
%  =========================
C = struct();
C.save_file = save_file;
C.source_event_file = source_event_file;
C.source_event_file_signature = source_event_file_signature;
C.dataset_id = cfg.dataset_id;
C.file_stem = cfg.file_stem;

C.band_labels = meta.band_labels;
C.band_role_indices = meta.band_role_indices;
C.selected_channels = get_optional_field(R, 'selected_channels', []);
C.channel_sites = meta.channel_sites;
C.regions = meta.regions;
C.channel_region_idx = meta.channel_region_idx;

C.session_ids = meta.session_ids;
C.session_lengths = meta.session_lengths;
C.session_dx = meta.session_dx;
C.session_start_idx = meta.session_start_idx;
C.session_end_idx = meta.session_end_idx;
C.dx_ref = meta.dx_ref;
C.L_total = meta.L_total;

C.min_channel_count = meta.min_channel_count;
C.majority_rule = sprintf('Band consensus requires at least %d/%d channels active at each sample.', ...
    meta.min_channel_count, meta.n_channels);
C.require_region_presence = logical(params.require_region_presence);
C.required_regions = params.required_regions(:);
C.state_definition_mode = char(string(params.state_definition_mode));
C.state_merge_gap_samples = double(params.state_merge_gap_samples);
C.state_definition_rule = sprintf([ ...
    'State windows are formed by merging overlapping consensus theta/gamma/ripple events within each session ' ...
    '(merge gap = %d samples) and assigning one label per merged window.'], ...
    double(params.state_merge_gap_samples));

C.consensus_mask_by_band = consensus_mask_by_band;
C.channel_support_count_by_band = channel_support_count_by_band;
C.region_support_count_by_band = region_support_count_by_band;
C.band_consensus = band_consensus;

C.state_catalog = state_catalog;
C.state_code_by_time = state_code_by_time;
C.state_windows = state_windows;
C.state_sample_count = state_sample_count;
C.state_duration_sec = state_duration_sec;
C.state_window_count = state_window_count;
C.params = params;

save(save_file, 'C', '-v7.3');
end


function meta = resolve_consensus_context(R, params)
[n_bands, n_channels] = size(R.DetectResults);
band_labels = get_band_labels(R, n_bands);

[session_ids, session_lengths, session_dx, session_start_idx, session_end_idx, ...
    dx_ref, L_total] = resolve_time_metadata(R);

channel_sites = get_channel_sites(R, n_channels);
[regions, channel_region_idx] = resolve_region_metadata(R, channel_sites, n_channels);

[theta_idx, gamma_idx, ripple_idx] = resolve_band_indices(band_labels, params, n_bands);
role_idx = [theta_idx, gamma_idx, ripple_idx];
if numel(unique(role_idx)) ~= 3
    error('Theta, gamma, and ripple must map to three distinct band indices.');
end

if isempty(params.min_channel_count)
    min_channel_count = floor(n_channels / 2) + 1;
else
    min_channel_count = double(params.min_channel_count);
end

if ~isscalar(min_channel_count) || min_channel_count < 1 || min_channel_count > n_channels
    error('params.min_channel_count must be empty or an integer between 1 and %d.', n_channels);
end

required_region_idx = zeros(1, 0);
if params.require_region_presence
    required_region_idx = resolve_required_regions(regions, params.required_regions);
end

meta = struct();
meta.n_bands = n_bands;
meta.n_channels = n_channels;
meta.n_regions = numel(regions);
meta.band_labels = band_labels;
meta.band_role_indices = struct('theta', theta_idx, 'gamma', gamma_idx, 'ripple', ripple_idx);
meta.theta_idx = theta_idx;
meta.gamma_idx = gamma_idx;
meta.ripple_idx = ripple_idx;
meta.session_ids = session_ids;
meta.session_lengths = session_lengths;
meta.session_dx = session_dx;
meta.session_start_idx = session_start_idx;
meta.session_end_idx = session_end_idx;
meta.dx_ref = dx_ref;
meta.L_total = L_total;
meta.channel_sites = channel_sites;
meta.regions = regions;
meta.channel_region_idx = channel_region_idx;
meta.min_channel_count = round(min_channel_count);
meta.required_region_idx = required_region_idx;
end


function [consensus_mask_by_band, channel_support_count_by_band, ...
    region_support_count_by_band, band_consensus] = ...
    build_band_consensus_outputs(DetectResults, meta, require_region_presence)
consensus_mask_by_band = false(meta.L_total, meta.n_bands);
channel_support_count_by_band = zeros(meta.L_total, meta.n_bands, 'uint16');
region_support_count_by_band = zeros(meta.L_total, meta.n_regions, meta.n_bands, 'uint16');
band_consensus = repmat(struct( ...
    'name', '', ...
    'band_index', [], ...
    'event_win', zeros(0, 2), ...
    'event_win_session_local', zeros(0, 2), ...
    'event_session_idx', zeros(0, 1), ...
    'event_session_id', zeros(0, 1), ...
    'duration_samples', zeros(0, 1), ...
    'duration_sec', zeros(0, 1), ...
    'support_count_max', zeros(0, 1), ...
    'support_count_mean', zeros(0, 1), ...
    'region_labels', {cell(0, 1)}, ...
    'region_support_max', zeros(0, 0), ...
    'region_support_mean', zeros(0, 0), ...
    'window_count', 0, ...
    'total_consensus_samples', 0, ...
    'total_consensus_duration_sec', 0), meta.n_bands, 1);

for b = 1:meta.n_bands
    channel_active = build_band_channel_mask(DetectResults(b, :), meta.L_total);
    band_support = sum(channel_active, 2);
    region_support = zeros(meta.L_total, meta.n_regions, 'uint16');

    for r = 1:meta.n_regions
        region_support(:, r) = uint16(sum(channel_active(:, meta.channel_region_idx == r), 2));
    end

    consensus_mask = band_support >= meta.min_channel_count;
    if require_region_presence
        consensus_mask = consensus_mask & all(region_support(:, meta.required_region_idx) > 0, 2);
    end

    consensus_mask_by_band(:, b) = consensus_mask;
    channel_support_count_by_band(:, b) = uint16(band_support);
    region_support_count_by_band(:, :, b) = region_support;

    band_consensus(b) = build_band_consensus_summary( ...
        meta.band_labels{b}, ...
        b, ...
        consensus_mask, ...
        uint16(band_support), ...
        region_support, ...
        meta.regions, ...
        meta.session_ids, ...
        meta.session_dx, ...
        meta.session_start_idx, ...
        meta.session_end_idx);
end
end


function [state_catalog, state_windows, state_code_by_time, ...
    state_sample_count, state_duration_sec, state_window_count] = ...
    build_state_outputs( ...
        band_consensus, consensus_mask_by_band, channel_support_count_by_band, ...
        region_support_count_by_band, meta, state_merge_gap_samples)
state_catalog = build_state_catalog();

[state_windows, state_code_by_time] = build_event_defined_state_windows( ...
    state_catalog, ...
    meta.band_labels, ...
    band_consensus, ...
    consensus_mask_by_band, ...
    channel_support_count_by_band, ...
    region_support_count_by_band, ...
    meta.regions, ...
    meta.theta_idx, ...
    meta.gamma_idx, ...
    meta.ripple_idx, ...
    meta.session_ids, ...
    meta.session_dx, ...
    meta.session_start_idx, ...
    state_merge_gap_samples);

[state_sample_count, state_duration_sec, state_window_count] = ...
    summarize_state_occupancy( ...
        state_code_by_time, ...
        state_catalog, ...
        meta.session_dx, ...
        meta.session_start_idx, ...
        meta.session_end_idx, ...
        state_windows);
end


function params = apply_default_params(params)
if ~isfield(params, 'min_channel_count')
    params.min_channel_count = [];
end

if ~isfield(params, 'theta_band_index')
    params.theta_band_index = [];
end

if ~isfield(params, 'gamma_band_index')
    params.gamma_band_index = [];
end

if ~isfield(params, 'ripple_band_index')
    params.ripple_band_index = [];
end

if ~isfield(params, 'require_region_presence')
    params.require_region_presence = false;
end

if ~isfield(params, 'required_regions') || isempty(params.required_regions)
    params.required_regions = {'hp', 'pl'};
end

if ~isfield(params, 'state_definition_mode') || isempty(params.state_definition_mode)
    params.state_definition_mode = 'event_windows';
end

if ~isfield(params, 'state_merge_gap_samples') || isempty(params.state_merge_gap_samples)
    params.state_merge_gap_samples = 0;
end

if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end
end


function validate_params(params)
mode_name = validatestring(char(string(params.state_definition_mode)), {'event_windows'});
if ~strcmp(mode_name, 'event_windows')
    error('Unsupported params.state_definition_mode: %s', char(string(params.state_definition_mode)));
end

gap_samples = double(params.state_merge_gap_samples);
if ~isscalar(gap_samples) || gap_samples < 0 || gap_samples ~= round(gap_samples)
    error('params.state_merge_gap_samples must be a nonnegative integer.');
end
end


function tf = can_reuse_saved_result(C, params, min_channel_count, band_role_indices, source_event_file_signature)
tf = isstruct(C);
if ~tf
    return;
end

required_fields = {'state_definition_mode', 'state_merge_gap_samples', 'min_channel_count', ...
    'require_region_presence', 'band_role_indices', 'source_event_file_signature'};
for i = 1:numel(required_fields)
    if ~isfield(C, required_fields{i}) || isempty(C.(required_fields{i}))
        tf = false;
        return;
    end
end

if ~isfield(C, 'state_definition_mode') || isempty(C.state_definition_mode)
    tf = false;
    return;
end

tf = strcmpi(char(string(C.state_definition_mode)), char(string(params.state_definition_mode)));
if ~tf
    return;
end

if ~isfield(C, 'state_merge_gap_samples') || isempty(C.state_merge_gap_samples)
    tf = false;
    return;
end

tf = double(C.state_merge_gap_samples) == double(params.state_merge_gap_samples);
if ~tf
    return;
end

tf = double(C.min_channel_count) == double(min_channel_count);
if ~tf
    return;
end

tf = logical(C.require_region_presence) == logical(params.require_region_presence);
if ~tf
    return;
end

if params.require_region_presence
    if ~isfield(C, 'required_regions') || isempty(C.required_regions)
        tf = false;
        return;
    end

    tf = isequal( ...
        normalize_region_list(C.required_regions), ...
        normalize_region_list(params.required_regions));
    if ~tf
        return;
    end
end

tf = isequal(double(C.band_role_indices.theta), double(band_role_indices.theta)) && ...
    isequal(double(C.band_role_indices.gamma), double(band_role_indices.gamma)) && ...
    isequal(double(C.band_role_indices.ripple), double(band_role_indices.ripple));
if ~tf
    return;
end

tf = io_utils.file_signature_matches(C.source_event_file_signature, source_event_file_signature);
end


function band_labels = get_band_labels(R, n_bands)
if isfield(R, 'params') && isfield(R.params, 'band_labels') && numel(R.params.band_labels) >= n_bands
    band_labels = R.params.band_labels(:);
else
    band_labels = cell(n_bands, 1);
    for b = 1:n_bands
        band_labels{b} = sprintf('band%d', b);
    end
end

for b = 1:n_bands
    band_labels{b} = char(string(band_labels{b}));
end
end


function [session_ids, session_lengths, session_dx, session_start_idx, session_end_idx, ...
    dx_ref, L_total] = resolve_time_metadata(R)
if isfield(R, 'session_ids') && isfield(R, 'session_lengths') && ...
        isfield(R, 'session_dx') && isfield(R, 'session_start_idx') && ...
        isfield(R, 'session_end_idx') && ~isempty(R.session_lengths)
    session_ids = double(R.session_ids(:));
    session_lengths = double(R.session_lengths(:));
    session_dx = double(R.session_dx(:));
    session_start_idx = double(R.session_start_idx(:));
    session_end_idx = double(R.session_end_idx(:));
    dx_ref = median(session_dx);
    L_total = sum(session_lengths);
    return;
end

[dx_ref, L_total] = infer_dt_and_length(R.DetectResults);
session_ids = 1;
session_lengths = L_total;
session_dx = dx_ref;
session_start_idx = 1;
session_end_idx = L_total;
end


function channel_sites = get_channel_sites(R, n_channels)
if isfield(R, 'channel_sites') && numel(R.channel_sites) >= n_channels
    channel_sites = R.channel_sites(:);
    return;
end

channel_sites = cell(n_channels, 1);
for c = 1:n_channels
    S = R.DetectResults{1, c};
    if ~isempty(S) && isfield(S, 'channel_site') && ~isempty(S.channel_site)
        channel_sites{c} = char(string(S.channel_site));
    else
        channel_sites{c} = '';
    end
end
end


function [regions, channel_region_idx] = resolve_region_metadata(R, channel_sites, n_channels)
if isfield(R, 'regions') && ~isempty(R.regions) && ...
        isfield(R, 'channel_region_idx') && numel(R.channel_region_idx) >= n_channels
    regions = R.regions(:);
    channel_region_idx = double(R.channel_region_idx(:))';
    return;
end

if ~isempty(channel_sites) && any(~cellfun(@isempty, channel_sites))
    regions = unique(channel_sites, 'stable');
    channel_region_idx = zeros(1, n_channels);
    for c = 1:n_channels
        ridx = find(strcmp(regions, channel_sites{c}), 1, 'first');
        if isempty(ridx)
            error('Could not map channel site "%s" to a region.', channel_sites{c});
        end
        channel_region_idx(c) = ridx;
    end
    return;
end

regions = {'all'};
channel_region_idx = ones(1, n_channels);
end


function [theta_idx, gamma_idx, ripple_idx] = resolve_band_indices(band_labels, params, n_bands)
theta_idx = resolve_single_band_index(params.theta_band_index, band_labels, 'theta', n_bands);
gamma_idx = resolve_single_band_index(params.gamma_band_index, band_labels, 'gamma', n_bands);
ripple_idx = resolve_single_band_index(params.ripple_band_index, band_labels, 'ripple', n_bands);
end


function idx = resolve_single_band_index(explicit_idx, band_labels, target_label, n_bands)
if ~isempty(explicit_idx)
    idx = double(explicit_idx);
    if ~isscalar(idx) || idx < 1 || idx > n_bands || idx ~= round(idx)
        error('Explicit %s band index must be an integer between 1 and %d.', target_label, n_bands);
    end
    return;
end

labels_norm = cellfun(@normalize_token, band_labels, 'UniformOutput', false);
target_norm = normalize_token(target_label);
idx = find(strcmp(labels_norm, target_norm), 1, 'first');

if isempty(idx)
    error('Could not infer the %s band from R.params.band_labels. Please set params.%s_band_index explicitly.', ...
        target_label, target_label);
end
end


function required_region_idx = resolve_required_regions(regions, required_regions)
required_region_idx = zeros(1, numel(required_regions));

for i = 1:numel(required_regions)
    wanted = normalize_token(required_regions{i});
    found = false;

    for r = 1:numel(regions)
        if strcmp(normalize_token(regions{r}), wanted)
            required_region_idx(i) = r;
            found = true;
            break;
        end
    end

    if ~found
        error('Required region "%s" is not present in R.regions.', char(string(required_regions{i})));
    end
end
end


function channel_active = build_band_channel_mask(DetectResultsBand, L_total)
n_channels = numel(DetectResultsBand);
channel_active = false(L_total, n_channels);

for c = 1:n_channels
    S = DetectResultsBand{c};
    if isempty(S) || ~isfield(S, 'event_win') || isempty(S.event_win)
        continue;
    end

    win = round(double(S.event_win));
    win(:, 1) = max(win(:, 1), 1);
    win(:, 2) = min(win(:, 2), L_total);
    win = win(win(:, 2) >= win(:, 1), :);

    if isempty(win)
        continue;
    end

    starts = win(:, 1);
    stops = win(:, 2) + 1;

    delta = zeros(L_total + 1, 1);
    delta = delta + accumarray(starts, 1, [L_total + 1, 1], @sum, 0);

    valid_stop = stops <= L_total;
    if any(valid_stop)
        delta = delta - accumarray(stops(valid_stop), 1, [L_total + 1, 1], @sum, 0);
    end

    channel_active(:, c) = cumsum(delta(1:L_total)) > 0;
end
end


function band_info = build_band_consensus_summary( ...
    band_label, band_index, consensus_mask, band_support, region_support, regions, ...
    session_ids, session_dx, session_start_idx, session_end_idx)
[event_win, event_win_session_local, event_session_idx, event_session_id] = ...
    logical_mask_to_windows(consensus_mask, session_ids, session_start_idx, session_end_idx);

n_windows = size(event_win, 1);
n_regions = size(region_support, 2);

duration_samples = zeros(n_windows, 1);
duration_sec = zeros(n_windows, 1);
support_count_max = zeros(n_windows, 1);
support_count_mean = zeros(n_windows, 1);
region_support_max = zeros(n_windows, n_regions);
region_support_mean = zeros(n_windows, n_regions);

for i = 1:n_windows
    g1 = event_win(i, 1);
    g2 = event_win(i, 2);

    duration_samples(i) = g2 - g1 + 1;
    duration_sec(i) = duration_samples(i) * double(session_dx(event_session_idx(i)));

    support_seg = double(band_support(g1:g2));
    support_count_max(i) = max(support_seg);
    support_count_mean(i) = mean(support_seg);

    region_support_max(i, :) = max(double(region_support(g1:g2, :)), [], 1);
    region_support_mean(i, :) = mean(double(region_support(g1:g2, :)), 1);
end

band_info = struct();
band_info.name = band_label;
band_info.band_index = band_index;
band_info.event_win = event_win;
band_info.event_win_session_local = event_win_session_local;
band_info.event_session_idx = event_session_idx;
band_info.event_session_id = event_session_id;
band_info.duration_samples = duration_samples;
band_info.duration_sec = duration_sec;
band_info.support_count_max = support_count_max;
band_info.support_count_mean = support_count_mean;
band_info.region_labels = regions(:);
band_info.region_support_max = region_support_max;
band_info.region_support_mean = region_support_mean;
band_info.window_count = n_windows;
band_info.total_consensus_samples = nnz(consensus_mask);
band_info.total_consensus_duration_sec = compute_mask_duration_sec( ...
    consensus_mask, session_dx, session_start_idx, session_end_idx);
end


function state_catalog = build_state_catalog()
state_catalog = struct( ...
    'code', num2cell(uint8([1; 2; 3; 4; 5])), ...
    'label', {'theta'; 'gamma'; 'ripple'; 'theta-gamma'; 'sharp-wave-ripple'}, ...
    'description', { ...
        'Merged consensus-event window containing theta only'; ...
        'Merged consensus-event window containing gamma only'; ...
        'Merged consensus-event window containing ripple, including gamma+ripple'; ...
        'Merged consensus-event window containing theta+gamma without ripple'; ...
        'Merged consensus-event window containing theta+ripple, optionally including gamma'});
end


function [state_windows, state_code_by_time] = build_event_defined_state_windows( ...
    state_catalog, band_labels, band_consensus, consensus_mask_by_band, ...
    channel_support_count_by_band, region_support_count_by_band, regions, ...
    theta_idx, gamma_idx, ripple_idx, session_ids, session_dx, session_start_idx, merge_gap_samples)
n_regions = numel(regions);
n_bands = numel(band_labels);
n_sessions = numel(session_ids);
n_time = size(consensus_mask_by_band, 1);

state_windows = struct([]);
state_code_by_time = zeros(n_time, 1, 'uint8');
idx_out = 0;
role_band_indices = [theta_idx, gamma_idx, ripple_idx];

for k = 1:n_sessions
    session_event_rows = collect_session_role_events(band_consensus, role_band_indices, k);
    if isempty(session_event_rows)
        continue;
    end

    merged_windows = merge_session_event_rows(session_event_rows, merge_gap_samples, n_bands);
    for i = 1:numel(merged_windows)
        g1 = merged_windows(i).start_idx;
        g2 = merged_windows(i).end_idx;

        band_any = any(consensus_mask_by_band(g1:g2, :), 1);
        has_theta = band_any(theta_idx);
        has_gamma = band_any(gamma_idx);
        has_ripple = band_any(ripple_idx);
        code = assign_state_code_from_presence(has_theta, has_gamma, has_ripple);
        if code == 0
            continue;
        end

        support_seg = double(channel_support_count_by_band(g1:g2, :));
        region_support_max = zeros(n_regions, n_bands);
        region_support_mean = zeros(n_regions, n_bands);
        for r = 1:n_regions
            for b = 1:n_bands
                trace_rb = double(region_support_count_by_band(g1:g2, r, b));
                region_support_max(r, b) = max(trace_rb);
                region_support_mean(r, b) = mean(trace_rb);
            end
        end

        idx_out = idx_out + 1;
        state_windows(idx_out).state_code = code;
        state_windows(idx_out).state_label = state_catalog(double(code)).label;
        state_windows(idx_out).win_global = [g1, g2];
        state_windows(idx_out).win_session_local = [ ...
            g1 - session_start_idx(k) + 1, ...
            g2 - session_start_idx(k) + 1];
        state_windows(idx_out).session_idx = k;
        state_windows(idx_out).session_id = session_ids(k);
        state_windows(idx_out).duration_samples = g2 - g1 + 1;
        state_windows(idx_out).duration_sec = (g2 - g1 + 1) * double(session_dx(k));
        state_windows(idx_out).consensus_band_any = band_any;
        state_windows(idx_out).consensus_band_all = all(consensus_mask_by_band(g1:g2, :), 1);
        state_windows(idx_out).consensus_band_labels_any = labels_from_mask(band_any, band_labels);
        state_windows(idx_out).consensus_band_labels_all = labels_from_mask(state_windows(idx_out).consensus_band_all, band_labels);
        state_windows(idx_out).channel_support_max_by_band = max(support_seg, [], 1);
        state_windows(idx_out).channel_support_mean_by_band = mean(support_seg, 1);
        state_windows(idx_out).region_labels = regions(:);
        state_windows(idx_out).region_support_max_by_band = region_support_max;
        state_windows(idx_out).region_support_mean_by_band = region_support_mean;
        state_windows(idx_out).region_support_max_any_band = max(region_support_max, [], 2);
        state_windows(idx_out).constituent_event_count_by_band = merged_windows(i).event_count_by_band;
        state_windows(idx_out).constituent_event_total_count = sum(merged_windows(i).event_count_by_band);

        state_code_by_time(g1:g2) = code;
    end
end
end


function session_event_rows = collect_session_role_events(band_consensus, role_band_indices, session_idx)
session_event_rows = zeros(0, 3);

for i = 1:numel(role_band_indices)
    band_idx = double(role_band_indices(i));
    band_info = band_consensus(band_idx);
    if isempty(band_info.event_win)
        continue;
    end

    keep_mask = double(band_info.event_session_idx(:)) == double(session_idx);
    event_win = double(band_info.event_win(keep_mask, :));
    if isempty(event_win)
        continue;
    end

    n_new = size(event_win, 1);
    session_event_rows = [session_event_rows; [event_win, repmat(band_idx, n_new, 1)]]; %#ok<AGROW>
end
end


function merged_windows = merge_session_event_rows(session_event_rows, merge_gap_samples, n_bands)
if isempty(session_event_rows)
    merged_windows = struct('start_idx', {}, 'end_idx', {}, 'event_count_by_band', {});
    return;
end

rows_sorted = sortrows(session_event_rows, [1, 2, 3]);
merged_windows = struct('start_idx', {}, 'end_idx', {}, 'event_count_by_band', {});

current_start = rows_sorted(1, 1);
current_end = rows_sorted(1, 2);
current_counts = zeros(1, n_bands);
current_counts(rows_sorted(1, 3)) = current_counts(rows_sorted(1, 3)) + 1;
idx_out = 0;

for i = 2:size(rows_sorted, 1)
    g1 = rows_sorted(i, 1);
    g2 = rows_sorted(i, 2);
    b = rows_sorted(i, 3);

    if g1 <= current_end + merge_gap_samples
        current_end = max(current_end, g2);
        current_counts(b) = current_counts(b) + 1;
    else
        idx_out = idx_out + 1;
        merged_windows(idx_out).start_idx = current_start;
        merged_windows(idx_out).end_idx = current_end;
        merged_windows(idx_out).event_count_by_band = current_counts;

        current_start = g1;
        current_end = g2;
        current_counts = zeros(1, n_bands);
        current_counts(b) = 1;
    end
end

idx_out = idx_out + 1;
merged_windows(idx_out).start_idx = current_start;
merged_windows(idx_out).end_idx = current_end;
merged_windows(idx_out).event_count_by_band = current_counts;
end


function code = assign_state_code_from_presence(has_theta, has_gamma, has_ripple)
code = uint8(0);

if has_theta && has_ripple
    code = uint8(5);
    return;
end

if has_theta && has_gamma
    code = uint8(4);
    return;
end

if has_theta
    code = uint8(1);
    return;
end

if has_gamma && ~has_ripple
    code = uint8(2);
    return;
end

if has_ripple
    code = uint8(3);
end
end


function [state_sample_count, state_duration_sec, state_window_count] = ...
    summarize_state_occupancy(state_code_by_time, state_catalog, session_dx, session_start_idx, session_end_idx, state_windows)
n_states = numel(state_catalog);
state_sample_count = zeros(n_states, 1);
state_duration_sec = zeros(n_states, 1);
state_window_count = zeros(n_states, 1);

for s = 1:n_states
    code = state_catalog(s).code;
    mask = (state_code_by_time == code);
    state_sample_count(s) = nnz(mask);
    state_duration_sec(s) = compute_mask_duration_sec(mask, session_dx, session_start_idx, session_end_idx);
end

if isempty(state_windows)
    return;
end

window_codes = double([state_windows.state_code]);
for s = 1:n_states
    state_window_count(s) = sum(window_codes == double(state_catalog(s).code));
end
end


function [event_win, event_win_session_local, event_session_idx, event_session_id] = ...
    logical_mask_to_windows(mask, session_ids, session_start_idx, session_end_idx)
n_sessions = numel(session_ids);

event_win = zeros(0, 2);
event_win_session_local = zeros(0, 2);
event_session_idx = zeros(0, 1);
event_session_id = zeros(0, 1);

for k = 1:n_sessions
    idx1 = session_start_idx(k);
    idx2 = session_end_idx(k);
    seg_mask = logical(mask(idx1:idx2));

    if ~any(seg_mask)
        continue;
    end

    d = diff([false; seg_mask(:); false]);
    start_local = find(d == 1);
    end_local = find(d == -1) - 1;
    n_new = numel(start_local);

    event_win = [event_win; [idx1 + start_local - 1, idx1 + end_local - 1]]; %#ok<AGROW>
    event_win_session_local = [event_win_session_local; [start_local, end_local]]; %#ok<AGROW>
    event_session_idx = [event_session_idx; repmat(k, n_new, 1)]; %#ok<AGROW>
    event_session_id = [event_session_id; repmat(session_ids(k), n_new, 1)]; %#ok<AGROW>
end
end


function duration_sec = compute_mask_duration_sec(mask, session_dx, session_start_idx, session_end_idx)
duration_sec = 0;

for k = 1:numel(session_dx)
    idx1 = session_start_idx(k);
    idx2 = session_end_idx(k);
    duration_sec = duration_sec + nnz(mask(idx1:idx2)) * double(session_dx(k));
end
end


function label_list = labels_from_mask(mask, labels)
label_list = labels(logical(mask(:)));
label_list = label_list(:);
end


function value = get_optional_field(S, field_name, default_value)
if isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function [dx_ref, L_total] = infer_dt_and_length(DetectResults)
dx_ref = [];
L_total = [];

for i = 1:numel(DetectResults)
    S = DetectResults{i};
    if isempty(S)
        continue;
    end

    if isfield(S, 'dx') && ~isempty(S.dx)
        dx_val = double(S.dx);
        if dx_val > 1
            dx_ref = 1 / dx_val;
        else
            dx_ref = dx_val;
        end
    end

    if isfield(S, 'L_total') && ~isempty(S.L_total)
        L_total = double(S.L_total);
    end

    if ~isempty(dx_ref) && ~isempty(L_total)
        return;
    end
end

if isempty(dx_ref) || isempty(L_total)
    error('Could not infer dt and total length from DetectResults.');
end
end


function token = normalize_token(label)
token = lower(strtrim(char(string(label))));
token = regexprep(token, '[^a-z0-9]+', '');
end


function labels = normalize_region_list(labels_in)
labels = cellfun(@normalize_token, cellstr(string(labels_in(:))), 'UniformOutput', false);
labels = labels(:);
end


function tag = build_save_tag(file_stem, min_channel_count, require_region_presence, required_regions)
tag = sprintf('%s_consensus_states_min%dch', file_stem, min_channel_count);

if require_region_presence
    region_tokens = cell(numel(required_regions), 1);
    for i = 1:numel(required_regions)
        region_tokens{i} = normalize_token(required_regions{i});
    end
    tag = sprintf('%s_regions_%s', tag, strjoin(region_tokens, '-'));
end
end
