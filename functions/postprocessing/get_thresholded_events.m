function [E, figs] = get_thresholded_events(phi, dt, evalues, params)
%GET_THRESHOLDED_EVENTS Detect thresholded eigenfunction events with NMS.
%
% Event candidates are contiguous above-threshold runs after optional short
% gap merging and minimum-duration filtering. Accepted events are selected
% with temporal NMS so that each mode has at most one event peak within the
% mode-specific NMS window. By default, that window is derived from the
% Koopman eigenvalue period.

if nargin < 2
    dt = [];
end

if nargin < 3
    evalues = [];
end

if nargin < 4 || isempty(params)
    if isstruct(evalues)
        params = evalues;
        evalues = [];
    else
        params = struct();
    end
end

params = local_apply_defaults(params);
local_validate_phi(phi);

if isempty(dt) || ~isnumeric(dt) || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('A positive scalar dt is required for thresholded event detection.');
end
dt = double(dt);

[X, value_transform_used] = local_transform_phi(phi, params.value_transform);
[T, K] = size(X);
evalues = local_evalues(evalues, K);
mode_index = local_mode_index(params, K);

[session_info, session_source] = local_resolve_sessions(params, T);
if params.require_session_metadata && strcmp(session_source, 'single_trace_fallback')
    error(['Thresholded event detection requires session metadata, but none ', ...
        'was found in params or params.observable_file.']);
end

[bin_len, step_len, bin_sec_used, step_sec_used] = ...
    local_resolve_density_window_lengths(params, dt, T);
[win_start, win_end, win_center, win_session_idx, win_session_id] = ...
    local_session_window_indices(T, bin_len, step_len, session_info);
[t_start, t_end, t_center] = local_window_times( ...
    win_start, win_end, win_center, dt, params.time_axis);

threshold_ratio = local_resolve_threshold_ratio(params.threshold_ratio, K);
threshold_mode = lower(char(params.threshold_mode));

threshold_by_mode = nan(1, K);
mode_timescales = local_empty_mode_timescale(K);
store = local_empty_event_store();
event_counter = 0;
event_count_by_mode = zeros(1, K);
candidate_count_by_mode = zeros(1, K);
suppressed_count_by_mode = zeros(1, K);

for k = 1:K
    x = double(X(:, k));
    valid_values = x(isfinite(x));
    threshold_by_mode(k) = local_compute_threshold( ...
        valid_values, threshold_ratio(k), threshold_mode);

    mode_timescales(k) = local_mode_timescale(evalues(k), dt, params.nms);

    if isnan(threshold_by_mode(k))
        continue;
    end

    mode_candidates_total = 0;
    mode_accepted_total = 0;

    for s = 1:numel(session_info.start_idx)
        idx1 = session_info.start_idx(s);
        idx2 = session_info.end_idx(s);
        x_session = x(idx1:idx2);
        high = isfinite(x_session) & x_session > threshold_by_mode(k);

        candidates = local_detect_session_candidates( ...
            x, x_session, high, idx1, threshold_by_mode(k), ...
            mode_timescales(k), dt, params);
        if isempty(candidates.peak_idx)
            continue;
        end

        mode_candidates_total = mode_candidates_total + numel(candidates.peak_idx);
        keep = local_keep_candidates(candidates, mode_timescales(k), params);
        keep_idx = find(keep(:));
        mode_accepted_total = mode_accepted_total + numel(keep_idx);

        for ii = 1:numel(keep_idx)
            j = keep_idx(ii);
            event_counter = event_counter + 1;
            store = local_append_event(store, event_counter, k, mode_index(k), ...
                local_get_selected_mode_idx(params, k), s, ...
                session_info.session_ids(s), candidates, j, ...
                threshold_by_mode(k), mode_timescales(k), params, dt);
        end
    end

    candidate_count_by_mode(k) = mode_candidates_total;
    event_count_by_mode(k) = mode_accepted_total;
    suppressed_count_by_mode(k) = mode_candidates_total - mode_accepted_total;
end

event_table = local_store_to_table(store);

[event_rate, event_count, event_occupancy, event_presence] = ...
    local_compute_event_density(event_table, K, win_start, win_end, dt, params.output_class);

if params.smooth_rate
    g = local_gaussian_kernel(params.smooth_sigma_sec, step_sec_used);
    smoothed_event_rate = zeros(size(event_rate), params.output_class);
    for k = 1:K
        smoothed_event_rate(:, k) = cast( ...
            conv(double(event_rate(:, k)), g, 'same'), params.output_class);
    end
else
    smoothed_event_rate = event_rate;
end

E = struct();
E.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
E.event_table = event_table;
E.event_rate_time_by_mode = event_rate;
E.smoothed_event_rate_time_by_mode = smoothed_event_rate;
E.event_count_time_by_mode = event_count;
E.event_occupancy_time_by_mode = event_occupancy;
E.event_presence_time_by_mode = event_presence;
E.threshold_by_mode = threshold_by_mode(:).';
E.mode_timescales = mode_timescales;
E.window_start_idx = win_start(:);
E.window_end_idx = win_end(:);
E.window_center_idx = win_center(:);
E.window_session_idx = win_session_idx(:);
E.window_session_id = win_session_id(:);
E.t_start = t_start(:);
E.t_end = t_end(:);
E.t_centers = t_center(:);
E.mode_index = mode_index(:);
E.selected_mode_idx_in_original = local_get_field(params, ...
    'selected_mode_idx_in_original', []);

E.input = struct();
E.input.n_samples = T;
E.input.n_modes = K;
E.input.dt = dt;
E.input.dt_source = local_get_field(params, 'dt_source', []);
E.input.value_transform_requested = params.value_transform;
E.input.value_transform_used = value_transform_used;
E.input.axis_order = 'time_by_mode';

E.session = session_info;
E.session.source = session_source;

E.params = params;
E.params.bin_samples = bin_len;
E.params.step_samples = step_len;
E.params.bin_sec_resolved = bin_sec_used;
E.params.step_sec_resolved = step_sec_used;
E.params.threshold_ratio_by_mode = threshold_ratio(:).';

E.summary = struct();
E.summary.n_events = height(event_table);
E.summary.n_candidate_events_by_mode = candidate_count_by_mode;
E.summary.n_events_by_mode = event_count_by_mode;
E.summary.n_suppressed_events_by_mode = suppressed_count_by_mode;
E.summary.mean_event_rate_by_mode = local_nanmean(double(event_rate), 1);
E.summary.mean_event_occupancy_by_mode = local_nanmean(double(event_occupancy), 1);
E.summary.max_event_rate_by_mode = local_nanmax(double(event_rate), 1);
E.summary.threshold_min = local_nanmin_vector(threshold_by_mode);
E.summary.threshold_max = local_nanmax_vector(threshold_by_mode);

E.artifacts = struct();
E.artifacts.mat_file = '';
E.artifacts.figure_file = '';

[E, mat_path, fig_path] = local_prepare_artifact_paths(E, params);
E.artifacts.mat_file = mat_path;
E.artifacts.figure_file = fig_path;

figs = struct();
figs.summary = [];
if params.make_figure || params.save_figure
    figs.summary = local_plot_events(E, params);
    if params.save_figure
        if exist(params.figure_dir, 'dir') ~= 7
            mkdir(params.figure_dir);
        end
        exportgraphics(figs.summary, fig_path, 'Resolution', params.figure_resolution);
    end
end

if params.save_results
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end
    event_table = E.event_table; %#ok<NASGU>
    event_rate = E.event_rate_time_by_mode; %#ok<NASGU>
    smoothed_event_rate = E.smoothed_event_rate_time_by_mode; %#ok<NASGU>
    event_occupancy = E.event_occupancy_time_by_mode; %#ok<NASGU>
    event_presence = E.event_presence_time_by_mode; %#ok<NASGU>
    thresholds = E.threshold_by_mode; %#ok<NASGU>
    mode_timescales = E.mode_timescales; %#ok<NASGU>
    t_centers = E.t_centers; %#ok<NASGU>
    params_saved = E.params; %#ok<NASGU>
    if params.save_v7_3
        save(mat_path, 'E', 'event_table', 'event_rate', ...
            'smoothed_event_rate', 'event_occupancy', 'event_presence', ...
            'thresholds', 'mode_timescales', 't_centers', ...
            'params_saved', '-v7.3');
    else
        save(mat_path, 'E', 'event_table', 'event_rate', ...
            'smoothed_event_rate', 'event_occupancy', 'event_presence', ...
            'thresholds', 'mode_timescales', 't_centers', ...
            'params_saved');
    end
end
end


function params = local_apply_defaults(params)
if ~isfield(params, 'threshold_ratio'), params.threshold_ratio = 0.7; end
if ~isfield(params, 'threshold_mode'), params.threshold_mode = 'meanplusstd'; end
if ~isfield(params, 'value_transform'), params.value_transform = 'auto'; end
if ~isfield(params, 'event_detector'), params.event_detector = 'find_peak_loc'; end
if ~isfield(params, 'event_score'), params.event_score = 'area_above_threshold'; end
if ~isfield(params, 'min_duration_sec'), params.min_duration_sec = 0.03; end
if ~isfield(params, 'merge_gap_sec'), params.merge_gap_sec = 0.02; end
if ~isfield(params, 'min_duration_samples'), params.min_duration_samples = []; end
if ~isfield(params, 'merge_gap_samples'), params.merge_gap_samples = []; end
if ~isfield(params, 'find_peak_loc_window_sec')
    params.find_peak_loc_window_sec = [];
end
if ~isfield(params, 'find_peak_loc_window_samples')
    params.find_peak_loc_window_samples = [];
end
if ~isfield(params, 'find_peak_loc_drop_first')
    params.find_peak_loc_drop_first = false;
end
if ~isfield(params, 'find_peak_loc_post_nms')
    params.find_peak_loc_post_nms = false;
end
if ~isfield(params, 'bin_sec'), params.bin_sec = 2; end
if ~isfield(params, 'step_sec'), params.step_sec = params.bin_sec; end
if ~isfield(params, 'bin_samples'), params.bin_samples = []; end
if ~isfield(params, 'step_samples'), params.step_samples = []; end
if ~isfield(params, 'smooth_rate'), params.smooth_rate = true; end
if ~isfield(params, 'smooth_sigma_sec'), params.smooth_sigma_sec = params.bin_sec; end
if ~isfield(params, 'output_class'), params.output_class = 'single'; end
if ~isfield(params, 'time_axis'), params.time_axis = []; end
if ~isfield(params, 'mode_index'), params.mode_index = []; end
if ~isfield(params, 'selected_mode_idx_in_original')
    params.selected_mode_idx_in_original = [];
end
if ~isfield(params, 'observable_file'), params.observable_file = ''; end
if ~isfield(params, 'require_session_metadata')
    params.require_session_metadata = false;
end

if ~isfield(params, 'session_start_idx'), params.session_start_idx = []; end
if ~isfield(params, 'session_end_idx'), params.session_end_idx = []; end
if ~isfield(params, 'session_lengths'), params.session_lengths = []; end
if ~isfield(params, 'session_ids'), params.session_ids = []; end
if ~isfield(params, 'session_dx'), params.session_dx = []; end

if ~isfield(params, 'nms') || ~isstruct(params.nms)
    params.nms = struct();
end
if ~isfield(params.nms, 'enable'), params.nms.enable = true; end
if ~isfield(params.nms, 'mode'), params.nms.mode = 'eigen_period'; end
if ~isfield(params.nms, 'period_fraction'), params.nms.period_fraction = 1.0; end
if ~isfield(params.nms, 'fixed_sec'), params.nms.fixed_sec = 0.25; end
if ~isfield(params.nms, 'default_sec'), params.nms.default_sec = 0.25; end
if ~isfield(params.nms, 'min_sec'), params.nms.min_sec = 0.05; end
if ~isfield(params.nms, 'max_sec'), params.nms.max_sec = 2.0; end
if ~isfield(params.nms, 'angle_eps'), params.nms.angle_eps = 1e-4; end

if ~isfield(params, 'save_results'), params.save_results = false; end
if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = pwd;
end
if ~isfield(params, 'save_stem') || isempty(params.save_stem)
    params.save_stem = 'thresholded_events';
end
if ~isfield(params, 'save_tag'), params.save_tag = ''; end
if ~isfield(params, 'save_v7_3'), params.save_v7_3 = true; end

if ~isfield(params, 'make_figure'), params.make_figure = false; end
if ~isfield(params, 'save_figure'), params.save_figure = false; end
if ~isfield(params, 'figure_dir') || isempty(params.figure_dir)
    params.figure_dir = params.save_dir;
end
if ~isfield(params, 'figure_visible') || isempty(params.figure_visible)
    if params.save_figure
        params.figure_visible = 'off';
    else
        params.figure_visible = 'on';
    end
end
if ~isfield(params, 'figure_position') || isempty(params.figure_position)
    params.figure_position = [100, 80, 1280, 860];
end
if ~isfield(params, 'figure_resolution'), params.figure_resolution = 200; end
if ~isfield(params, 'max_plot_modes') || isempty(params.max_plot_modes)
    params.max_plot_modes = 80;
end
if ~isfield(params, 'max_plot_events') || isempty(params.max_plot_events)
    params.max_plot_events = 5000;
end
if ~isfield(params, 'title') || isempty(params.title)
    params.title = 'Thresholded Eigenfunction Events';
end
if ~isfield(params, 'colormap') || isempty(params.colormap)
    params.colormap = 'turbo';
end
if ~isfield(params, 'background_color') || isempty(params.background_color)
    params.background_color = [0, 0, 0];
end
if ~isfield(params, 'axes_color') || isempty(params.axes_color)
    params.axes_color = params.background_color;
end
if ~isfield(params, 'text_color') || isempty(params.text_color)
    params.text_color = [1, 1, 1];
end
if ~isfield(params, 'grid_color') || isempty(params.grid_color)
    params.grid_color = [0.75, 0.75, 0.75];
end
end


function local_validate_phi(phi)
if ~isnumeric(phi) || ndims(phi) ~= 2 || isempty(phi)
    error('phi must be a nonempty numeric [T x K] matrix.');
end
end


function [X, value_transform_used] = local_transform_phi(phi, value_transform)
value_transform = lower(char(value_transform));

switch value_transform
    case 'auto'
        if ~isreal(phi)
            X = abs(phi);
            value_transform_used = 'abs';
        else
            X = phi;
            value_transform_used = 'none';
        end
    case {'none', 'asis', 'as_is'}
        if ~isreal(phi)
            error(['value_transform="none" requires real-valued phi. ', ...
                'Use "abs", "real", or "auto" for complex eigenfunctions.']);
        end
        X = phi;
        value_transform_used = 'none';
    case 'abs'
        X = abs(phi);
        value_transform_used = 'abs';
    case 'real'
        X = real(phi);
        value_transform_used = 'real';
    otherwise
        error('Unknown value_transform = %s.', value_transform);
end
end


function evalues = local_evalues(evalues, K)
if isempty(evalues)
    evalues = nan(K, 1);
else
    evalues = evalues(:);
    if numel(evalues) ~= K
        error('evalues must be empty or contain one value per mode.');
    end
end
end


function mode_index = local_mode_index(params, K)
if isfield(params, 'mode_index') && ~isempty(params.mode_index)
    mode_index = params.mode_index(:);
    if numel(mode_index) ~= K
        error('params.mode_index must have one entry per mode.');
    end
else
    mode_index = (1:K).';
end
end


function selected_idx = local_get_selected_mode_idx(params, k)
selected_idx_all = local_get_field(params, 'selected_mode_idx_in_original', []);
if isempty(selected_idx_all)
    selected_idx = NaN;
else
    selected_idx = selected_idx_all(k);
end
end


function ratio = local_resolve_threshold_ratio(threshold_ratio, K)
if ~isnumeric(threshold_ratio) || isempty(threshold_ratio)
    error('params.threshold_ratio must be a numeric scalar or 1xK vector.');
end

if isscalar(threshold_ratio)
    ratio = repmat(double(threshold_ratio), 1, K);
else
    ratio = double(threshold_ratio(:)).';
    if numel(ratio) ~= K
        error('params.threshold_ratio must be scalar or have one value per mode.');
    end
end
end


function thr = local_compute_threshold(vv, ratio, mode_str)
if isempty(vv)
    thr = NaN;
    return;
end

switch mode_str
    case 'quantile'
        if ratio <= 0 || ratio >= 1
            error('For quantile mode, threshold_ratio must be in (0,1).');
        end
        thr = quantile(vv, ratio);

    case 'maxfrac'
        if ratio < 0
            error('For maxfrac mode, threshold_ratio must be >= 0.');
        end
        thr = ratio * max(vv);

    case 'meanplusstd'
        thr = mean(vv) + ratio * std(vv, 0);

    otherwise
        error(['Unknown threshold_mode = %s. Use quantile, maxfrac, ', ...
            'or meanplusstd.'], mode_str);
end
end


function [session_info, source] = local_resolve_sessions(params, T)
session_info = struct();
session_info.start_idx = params.session_start_idx(:);
session_info.end_idx = params.session_end_idx(:);
session_info.session_ids = params.session_ids(:);
session_info.session_dx = params.session_dx(:);
source = 'params';

if (isempty(session_info.start_idx) || isempty(session_info.end_idx)) && ...
        isfield(params, 'session_lengths') && ~isempty(params.session_lengths)
    lengths = round(double(params.session_lengths(:)));
    lengths = lengths(isfinite(lengths) & lengths > 0);
    if ~isempty(lengths)
        session_info.end_idx = cumsum(lengths);
        session_info.start_idx = [1; session_info.end_idx(1:end-1) + 1];
        source = 'params.session_lengths';
    end
end

if isempty(session_info.start_idx) || isempty(session_info.end_idx)
    [session_info, loaded] = local_load_sessions_from_observable(params.observable_file);
    if loaded
        source = 'observable_file';
    end
end

if isempty(session_info.start_idx) || isempty(session_info.end_idx)
    session_info.start_idx = 1;
    session_info.end_idx = T;
    session_info.session_ids = 1;
    session_info.session_dx = [];
    source = 'single_trace_fallback';
end

if isempty(session_info.session_ids)
    session_info.session_ids = (1:numel(session_info.start_idx)).';
end

if numel(session_info.start_idx) ~= numel(session_info.end_idx)
    error('session_start_idx and session_end_idx must have the same length.');
end

session_info.start_idx = round(double(session_info.start_idx(:)));
session_info.end_idx = round(double(session_info.end_idx(:)));
session_info.session_ids = double(session_info.session_ids(:));

if numel(session_info.session_ids) ~= numel(session_info.start_idx)
    error('session_ids must have one value per session.');
end

bad = session_info.start_idx < 1 | session_info.end_idx > T | ...
    session_info.end_idx < session_info.start_idx;
if any(bad)
    error('Session boundaries are invalid for a trace with T=%d samples.', T);
end
end


function [session_info, loaded] = local_load_sessions_from_observable(observable_file)
loaded = false;
session_info = struct('start_idx', [], 'end_idx', [], ...
    'session_ids', [], 'session_dx', []);

if isempty(observable_file) || exist(observable_file, 'file') ~= 2
    return;
end

try
    vars = who('-file', observable_file);
catch
    return;
end

fields = intersect({'session_start_idx', 'session_end_idx', ...
    'session_lengths', 'session_ids', 'session_dx'}, vars, 'stable');
if isempty(fields)
    return;
end

S = load(observable_file, fields{:});
if isfield(S, 'session_start_idx') && isfield(S, 'session_end_idx')
    session_info.start_idx = S.session_start_idx(:);
    session_info.end_idx = S.session_end_idx(:);
elseif isfield(S, 'session_lengths')
    lengths = double(S.session_lengths(:));
    session_info.end_idx = cumsum(lengths);
    session_info.start_idx = [1; session_info.end_idx(1:end-1) + 1];
else
    return;
end

if isfield(S, 'session_ids')
    session_info.session_ids = S.session_ids(:);
end
if isfield(S, 'session_dx')
    session_info.session_dx = S.session_dx(:);
end
loaded = true;
end


function samples = local_sec_to_samples(sec_value, sample_value, dt)
if ~isempty(sample_value)
    samples = max(0, round(double(sample_value)));
else
    samples = max(0, round(double(sec_value) / dt));
end
end


function runs = local_mask_to_runs(mask)
mask = logical(mask(:));
edge = diff([false; mask; false]);
runs = [find(edge == 1), find(edge == -1) - 1];
end


function runs_out = local_merge_runs(runs, merge_gap_samples)
if isempty(runs)
    runs_out = runs;
    return;
end

runs_out = zeros(size(runs));
n = 1;
runs_out(n, :) = runs(1, :);

for i = 2:size(runs, 1)
    gap = runs(i, 1) - runs_out(n, 2) - 1;
    if gap <= merge_gap_samples
        runs_out(n, 2) = runs(i, 2);
    else
        n = n + 1;
        runs_out(n, :) = runs(i, :);
    end
end

runs_out = runs_out(1:n, :);
end


function runs = local_filter_min_duration(runs, min_duration_samples)
if isempty(runs)
    return;
end
duration = runs(:, 2) - runs(:, 1) + 1;
runs = runs(duration >= min_duration_samples, :);
end


function candidates = local_detect_session_candidates( ...
    x, x_session, high, session_start_idx, threshold, timescale, dt, params)
detector = lower(char(params.event_detector));

switch detector
    case {'find_peak_loc', 'find_peaks', 'find_peak'}
        candidates = local_find_peak_loc_candidates( ...
            x, x_session, high, session_start_idx, threshold, ...
            timescale, dt, params);

    case {'runs', 'episodes', 'threshold_runs'}
        runs = local_threshold_runs(high, params, dt);
        candidates = local_runs_to_candidates( ...
            x, runs, session_start_idx, threshold, dt, params.event_score);

    otherwise
        error(['Unknown params.event_detector = %s. Use find_peak_loc ', ...
            'or runs.'], params.event_detector);
end
end


function keep = local_keep_candidates(candidates, timescale, params)
detector = lower(char(params.event_detector));

switch detector
    case {'find_peak_loc', 'find_peaks', 'find_peak'}
        if params.find_peak_loc_post_nms
            keep = local_nms_keep(candidates.peak_idx, candidates.score, ...
                timescale.nms_window_samples, params.nms.enable);
        else
            keep = true(numel(candidates.peak_idx), 1);
        end

    otherwise
        keep = local_nms_keep(candidates.peak_idx, candidates.score, ...
            timescale.nms_window_samples, params.nms.enable);
end
end


function runs = local_threshold_runs(high, params, dt)
runs = local_mask_to_runs(high);
runs = local_merge_runs(runs, ...
    local_sec_to_samples(params.merge_gap_sec, params.merge_gap_samples, dt));
runs = local_filter_min_duration(runs, ...
    local_sec_to_samples(params.min_duration_sec, params.min_duration_samples, dt));
end


function candidates = local_find_peak_loc_candidates( ...
    x, x_session, high, session_start_idx, threshold, timescale, dt, params)
if exist('find_peak_loc', 'file') ~= 2
    error(['params.event_detector = ''find_peak_loc'' requires ', ...
        'utils/find_peak_loc.m on the MATLAB path.']);
end

runs = local_threshold_runs(high, params, dt);
if isempty(runs)
    candidates = local_empty_candidates();
    return;
end

loc_session = find(high);
if params.find_peak_loc_drop_first && ~isempty(loc_session)
    loc_session(1) = [];
end

if isempty(loc_session)
    candidates = local_empty_candidates();
    return;
end

L = local_find_peak_loc_window_samples(timescale, params, dt);
if numel(loc_session) == 1
    peak_loc_local = loc_session(:);
else
    peak_loc_local = find_peak_loc(double(x_session(:)).', loc_session(:).', L);
    peak_loc_local = peak_loc_local(:);
end

if isempty(peak_loc_local) && ~isempty(loc_session)
    [~, best_idx] = max(double(x_session(loc_session)));
    peak_loc_local = loc_session(best_idx);
end

peak_loc_local = unique(round(double(peak_loc_local(:))));
peak_loc_local = peak_loc_local(isfinite(peak_loc_local) & ...
    peak_loc_local >= 1 & peak_loc_local <= numel(x_session));

if isempty(peak_loc_local)
    candidates = local_empty_candidates();
    return;
end

candidates = local_peak_locs_to_candidates( ...
    x, peak_loc_local, runs, session_start_idx, threshold, ...
    dt, params.event_score);
end


function L = local_find_peak_loc_window_samples(timescale, params, dt)
if ~isempty(params.find_peak_loc_window_samples)
    L = max(1, round(double(params.find_peak_loc_window_samples)));
elseif ~isempty(params.find_peak_loc_window_sec)
    L = max(1, round(double(params.find_peak_loc_window_sec) / dt));
elseif isfield(timescale, 'nms_window_samples') && ...
        ~isempty(timescale.nms_window_samples) && ...
        isfinite(timescale.nms_window_samples)
    L = max(1, round(double(timescale.nms_window_samples)));
else
    L = max(1, round(double(params.nms.default_sec) / dt));
end
end


function candidates = local_empty_candidates()
candidates = struct();
candidates.start_idx = zeros(0, 1);
candidates.end_idx = zeros(0, 1);
candidates.peak_idx = zeros(0, 1);
candidates.duration_sec = zeros(0, 1);
candidates.peak_value = zeros(0, 1);
candidates.mean_value = zeros(0, 1);
candidates.area_above_threshold = zeros(0, 1);
candidates.score = zeros(0, 1);
end


function candidates = local_runs_to_candidates(x, runs_local, session_start_idx, threshold, dt, score_mode)
n = size(runs_local, 1);
candidates = local_empty_candidates();

candidates.start_idx = zeros(n, 1);
candidates.end_idx = zeros(n, 1);
candidates.peak_idx = zeros(n, 1);
candidates.duration_sec = zeros(n, 1);
candidates.peak_value = zeros(n, 1);
candidates.mean_value = zeros(n, 1);
candidates.area_above_threshold = zeros(n, 1);
candidates.score = zeros(n, 1);

n_keep = 0;
for i = 1:n
    g1 = session_start_idx + runs_local(i, 1) - 1;
    g2 = session_start_idx + runs_local(i, 2) - 1;
    seg = x(g1:g2);
    finite_mask = isfinite(seg);
    if ~any(finite_mask)
        continue;
    end

    finite_idx = find(finite_mask);
    [peak_value, rel_pos] = max(seg(finite_mask));
    peak_idx = g1 + finite_idx(rel_pos) - 1;
    above = max(seg(finite_mask) - threshold, 0);
    area = sum(above) * dt;
    mean_value = mean(seg(finite_mask));
    duration_sec = (g2 - g1 + 1) * dt;

    switch lower(char(score_mode))
        case 'peak_value'
            score = peak_value;
        case 'duration_sec'
            score = duration_sec;
        case {'area', 'area_above_threshold'}
            score = area;
        otherwise
            error('Unknown params.event_score = %s.', score_mode);
    end

    n_keep = n_keep + 1;
    candidates.start_idx(n_keep) = g1;
    candidates.end_idx(n_keep) = g2;
    candidates.peak_idx(n_keep) = peak_idx;
    candidates.duration_sec(n_keep) = duration_sec;
    candidates.peak_value(n_keep) = peak_value;
    candidates.mean_value(n_keep) = mean_value;
    candidates.area_above_threshold(n_keep) = area;
    candidates.score(n_keep) = score;
end

fields = fieldnames(candidates);
for i = 1:numel(fields)
    candidates.(fields{i}) = candidates.(fields{i})(1:n_keep);
end
end


function candidates = local_peak_locs_to_candidates( ...
    x, peak_loc_local, runs_local, session_start_idx, threshold, dt, score_mode)
n = numel(peak_loc_local);
candidates = local_empty_candidates();

candidates.start_idx = zeros(n, 1);
candidates.end_idx = zeros(n, 1);
candidates.peak_idx = zeros(n, 1);
candidates.duration_sec = zeros(n, 1);
candidates.peak_value = zeros(n, 1);
candidates.mean_value = zeros(n, 1);
candidates.area_above_threshold = zeros(n, 1);
candidates.score = zeros(n, 1);

n_keep = 0;
for i = 1:n
    peak_local = peak_loc_local(i);
    run_idx = find(runs_local(:, 1) <= peak_local & ...
        runs_local(:, 2) >= peak_local, 1, 'first');
    if isempty(run_idx)
        continue;
    end

    g1 = session_start_idx + runs_local(run_idx, 1) - 1;
    g2 = session_start_idx + runs_local(run_idx, 2) - 1;
    peak_idx = session_start_idx + peak_local - 1;
    seg = x(g1:g2);
    finite_mask = isfinite(seg);
    if ~any(finite_mask) || ~isfinite(x(peak_idx))
        continue;
    end

    above = max(seg(finite_mask) - threshold, 0);
    area = sum(above) * dt;
    mean_value = mean(seg(finite_mask));
    duration_sec = (g2 - g1 + 1) * dt;
    peak_value = x(peak_idx);

    switch lower(char(score_mode))
        case 'peak_value'
            score = peak_value;
        case 'duration_sec'
            score = duration_sec;
        case {'area', 'area_above_threshold'}
            score = area;
        otherwise
            error('Unknown params.event_score = %s.', score_mode);
    end

    n_keep = n_keep + 1;
    candidates.start_idx(n_keep) = g1;
    candidates.end_idx(n_keep) = g2;
    candidates.peak_idx(n_keep) = peak_idx;
    candidates.duration_sec(n_keep) = duration_sec;
    candidates.peak_value(n_keep) = peak_value;
    candidates.mean_value(n_keep) = mean_value;
    candidates.area_above_threshold(n_keep) = area;
    candidates.score(n_keep) = score;
end

fields = fieldnames(candidates);
for i = 1:numel(fields)
    candidates.(fields{i}) = candidates.(fields{i})(1:n_keep);
end
end


function keep = local_nms_keep(peak_idx, score, nms_window_samples, enable_nms)
n = numel(peak_idx);
keep = true(n, 1);
if ~enable_nms || n <= 1 || isempty(nms_window_samples) || nms_window_samples <= 1
    return;
end

score = double(score(:));
score(~isfinite(score)) = -Inf;
[~, order] = sortrows([-score, double(peak_idx(:))], [1, 2]);
keep = false(n, 1);
accepted_peaks = zeros(0, 1);

for ii = 1:numel(order)
    idx = order(ii);
    if isempty(accepted_peaks) || all(abs(double(peak_idx(idx)) - accepted_peaks) > nms_window_samples)
        keep(idx) = true;
        accepted_peaks(end+1, 1) = double(peak_idx(idx)); %#ok<AGROW>
    end
end
end


function ts = local_empty_mode_timescale(K)
template = struct( ...
    'evalue_discrete', NaN, ...
    'period_sec', NaN, ...
    'decay_sec', NaN, ...
    'nms_window_sec', NaN, ...
    'nms_window_samples', NaN, ...
    'nms_source', '');
ts = repmat(template, K, 1);
end


function ts = local_mode_timescale(lambda, dt, nms)
ts = local_empty_mode_timescale(1);
ts.evalue_discrete = lambda;

rho = abs(lambda);
theta = abs(angle(lambda));
if isfinite(theta) && theta > nms.angle_eps
    ts.period_sec = 2 * pi * dt / theta;
end
if isfinite(rho) && rho > 0 && abs(log(rho)) > eps
    ts.decay_sec = -dt / log(rho);
end

mode = lower(char(nms.mode));
switch mode
    case 'fixed'
        win_sec = nms.fixed_sec;
        source = 'fixed';

    case {'eigen_period', 'period'}
        if isfinite(ts.period_sec) && ts.period_sec > 0
            win_sec = nms.period_fraction * ts.period_sec;
            source = 'eigen_period';
        else
            win_sec = nms.default_sec;
            source = 'default_nonoscillatory';
        end

    otherwise
        error('Unknown params.nms.mode = %s.', nms.mode);
end

win_sec = min(max(win_sec, nms.min_sec), nms.max_sec);
ts.nms_window_sec = win_sec;
ts.nms_window_samples = max(1, round(win_sec / dt));
ts.nms_source = source;
end


function store = local_empty_event_store()
store = struct();
store.event_id = zeros(0, 1);
store.mode_col = zeros(0, 1);
store.mode_index = zeros(0, 1);
store.selected_mode_idx_in_original = zeros(0, 1);
store.session_idx = zeros(0, 1);
store.session_id = zeros(0, 1);
store.start_idx = zeros(0, 1);
store.end_idx = zeros(0, 1);
store.peak_idx = zeros(0, 1);
store.start_time = zeros(0, 1);
store.end_time = zeros(0, 1);
store.peak_time = zeros(0, 1);
store.duration_sec = zeros(0, 1);
store.threshold = zeros(0, 1);
store.peak_value = zeros(0, 1);
store.mean_value = zeros(0, 1);
store.area_above_threshold = zeros(0, 1);
store.score = zeros(0, 1);
store.nms_window_sec = zeros(0, 1);
store.period_sec = zeros(0, 1);
store.decay_sec = zeros(0, 1);
store.nms_source = cell(0, 1);
store.event_score = cell(0, 1);
end


function store = local_append_event(store, event_id, mode_col, mode_index, ...
    selected_mode_idx, session_idx, session_id, candidates, j, threshold, ...
    timescale, params, dt)
store.event_id(end+1, 1) = event_id;
store.mode_col(end+1, 1) = mode_col;
store.mode_index(end+1, 1) = mode_index;
store.selected_mode_idx_in_original(end+1, 1) = selected_mode_idx;
store.session_idx(end+1, 1) = session_idx;
store.session_id(end+1, 1) = session_id;
store.start_idx(end+1, 1) = candidates.start_idx(j);
store.end_idx(end+1, 1) = candidates.end_idx(j);
store.peak_idx(end+1, 1) = candidates.peak_idx(j);
store.start_time(end+1, 1) = (candidates.start_idx(j) - 1) * dt;
store.end_time(end+1, 1) = (candidates.end_idx(j) - 1) * dt;
store.peak_time(end+1, 1) = (candidates.peak_idx(j) - 1) * dt;
store.duration_sec(end+1, 1) = candidates.duration_sec(j);
store.threshold(end+1, 1) = threshold;
store.peak_value(end+1, 1) = candidates.peak_value(j);
store.mean_value(end+1, 1) = candidates.mean_value(j);
store.area_above_threshold(end+1, 1) = candidates.area_above_threshold(j);
store.score(end+1, 1) = candidates.score(j);
store.nms_window_sec(end+1, 1) = timescale.nms_window_sec;
store.period_sec(end+1, 1) = timescale.period_sec;
store.decay_sec(end+1, 1) = timescale.decay_sec;
store.nms_source{end+1, 1} = timescale.nms_source;
store.event_score{end+1, 1} = char(params.event_score);
end


function event_table = local_store_to_table(store)
event_table = table( ...
    store.event_id, ...
    store.mode_col, ...
    store.mode_index, ...
    store.selected_mode_idx_in_original, ...
    store.session_idx, ...
    store.session_id, ...
    store.start_idx, ...
    store.end_idx, ...
    store.peak_idx, ...
    store.start_time, ...
    store.end_time, ...
    store.peak_time, ...
    store.duration_sec, ...
    store.threshold, ...
    store.peak_value, ...
    store.mean_value, ...
    store.area_above_threshold, ...
    store.score, ...
    store.nms_window_sec, ...
    store.period_sec, ...
    store.decay_sec, ...
    store.nms_source, ...
    store.event_score, ...
    'VariableNames', { ...
    'event_id', ...
    'mode_col', ...
    'mode_index', ...
    'selected_mode_idx_in_original', ...
    'session_idx', ...
    'session_id', ...
    'start_idx', ...
    'end_idx', ...
    'peak_idx', ...
    'start_time', ...
    'end_time', ...
    'peak_time', ...
    'duration_sec', ...
    'threshold', ...
    'peak_value', ...
    'mean_value', ...
    'area_above_threshold', ...
    'score', ...
    'nms_window_sec', ...
    'period_sec', ...
    'decay_sec', ...
    'nms_source', ...
    'event_score'});
end


function [event_rate, event_count, event_occupancy, event_presence] = ...
    local_compute_event_density(event_table, K, win_start, win_end, dt, output_class)
Nwin = numel(win_start);
event_count = zeros(Nwin, K, output_class);
event_rate = zeros(Nwin, K, output_class);
event_occupancy = zeros(Nwin, K, output_class);
event_presence = zeros(Nwin, K, output_class);
win_len = double(win_end - win_start + 1);
win_sec = win_len * dt;

if isempty(event_table)
    return;
end

for k = 1:K
    rows = event_table.mode_col == k;
    if ~any(rows)
        continue;
    end

    peaks = double(event_table.peak_idx(rows));
    starts = double(event_table.start_idx(rows));
    ends = double(event_table.end_idx(rows));

    counts = zeros(Nwin, 1);
    coverage = zeros(Nwin, 1);
    for w = 1:Nwin
        counts(w) = sum(peaks >= win_start(w) & peaks <= win_end(w));
        overlaps = max(0, min(ends, win_end(w)) - max(starts, win_start(w)) + 1);
        coverage(w) = min(sum(overlaps), win_len(w));
    end

    event_count(:, k) = cast(counts, output_class);
    event_rate(:, k) = cast(counts ./ win_sec, output_class);
    event_occupancy(:, k) = cast(coverage ./ win_len, output_class);
    event_presence(:, k) = cast(counts > 0, output_class);
end
end


function [win_len, step_len, bin_sec_used, step_sec_used] = ...
    local_resolve_density_window_lengths(params, dt, T)
if ~isempty(params.bin_samples)
    win_len = max(1, round(double(params.bin_samples)));
else
    win_len = max(1, round(double(params.bin_sec) / dt));
end

if ~isempty(params.step_samples)
    step_len = max(1, round(double(params.step_samples)));
else
    step_len = max(1, round(double(params.step_sec) / dt));
end

win_len = min(win_len, T);
step_len = max(1, step_len);
bin_sec_used = win_len * dt;
step_sec_used = step_len * dt;
end


function [start_idx, end_idx, center_idx] = local_window_indices(T, win_len, step_len)
last_start = T - win_len + 1;
start_idx = (1:step_len:last_start).';
if isempty(start_idx)
    start_idx = 1;
end
end_idx = start_idx + win_len - 1;
center_idx = start_idx + (win_len - 1) / 2;
end


function [start_idx, end_idx, center_idx, session_idx, session_id] = ...
    local_session_window_indices(T, win_len, step_len, session_info)
start_idx = zeros(0, 1);
end_idx = zeros(0, 1);
center_idx = zeros(0, 1);
session_idx = zeros(0, 1);
session_id = zeros(0, 1);

if isempty(session_info.start_idx) || isempty(session_info.end_idx)
    [start_idx, end_idx, center_idx] = local_window_indices(T, win_len, step_len);
    session_idx = ones(numel(start_idx), 1);
    session_id = ones(numel(start_idx), 1);
    return;
end

for s = 1:numel(session_info.start_idx)
    idx1 = session_info.start_idx(s);
    idx2 = session_info.end_idx(s);
    session_len = idx2 - idx1 + 1;
    this_win_len = min(win_len, session_len);
    if this_win_len <= 0
        continue;
    end

    local_last_start = session_len - this_win_len + 1;
    local_start = (1:step_len:local_last_start).';
    if isempty(local_start)
        local_start = 1;
    end

    this_start = idx1 + local_start - 1;
    this_end = this_start + this_win_len - 1;
    this_center = this_start + (this_win_len - 1) / 2;
    n_win = numel(this_start);

    start_idx = [start_idx; this_start(:)]; %#ok<AGROW>
    end_idx = [end_idx; this_end(:)]; %#ok<AGROW>
    center_idx = [center_idx; this_center(:)]; %#ok<AGROW>
    session_idx = [session_idx; repmat(s, n_win, 1)]; %#ok<AGROW>
    session_id = [session_id; repmat(session_info.session_ids(s), n_win, 1)]; %#ok<AGROW>
end

if isempty(start_idx)
    [start_idx, end_idx, center_idx] = local_window_indices(T, win_len, step_len);
    session_idx = ones(numel(start_idx), 1);
    session_id = ones(numel(start_idx), 1);
end
end


function [t_start, t_end, t_center] = local_window_times( ...
    start_idx, end_idx, center_idx, dt, time_axis)
if ~isempty(time_axis)
    time_axis = time_axis(:);
    if numel(time_axis) < max(end_idx)
        error('params.time_axis must have at least T samples.');
    end
    t_start = time_axis(start_idx);
    t_end = time_axis(end_idx);
    t_center = interp1((1:numel(time_axis)).', time_axis, center_idx, ...
        'linear', 'extrap');
else
    t_start = (start_idx - 1) * dt;
    t_end = (end_idx - 1) * dt;
    t_center = (center_idx - 1) * dt;
end
end


function g = local_gaussian_kernel(sigma_sec, step_sec)
if isempty(sigma_sec) || sigma_sec <= 0
    g = 1;
    return;
end

sigma_bins = sigma_sec / step_sec;
half_width = max(1, ceil(4 * sigma_bins));
x = (-half_width:half_width);
g = exp(-(x .^ 2) / (2 * sigma_bins ^ 2));
g = g / sum(g);
end


function [E, mat_path, fig_path] = local_prepare_artifact_paths(E, params)
stem = char(params.save_stem);
tag = char(params.save_tag);
ratio_tag = local_ratio_tag(params.threshold_ratio);
pieces = {stem, char(params.threshold_mode), ratio_tag, 'events'};
if ~isempty(tag)
    pieces{end+1} = tag;
end

base = strjoin(pieces, '__');
base = regexprep(base, '[^\w\-]+', '_');
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
filename = sprintf('%s__%s', base, timestamp);

mat_path = '';
fig_path = '';
if params.save_results
    mat_path = fullfile(params.save_dir, [filename, '.mat']);
end
if params.save_figure
    fig_path = fullfile(params.figure_dir, [filename, '.png']);
end

E.artifacts.file_stem = filename;
end


function tag = local_ratio_tag(ratio)
if isscalar(ratio)
    tag = sprintf('ratio_%03d', round(double(ratio) * 100));
else
    tag = sprintf('ratio_vector%d', numel(ratio));
end
end


function fig = local_plot_events(E, params)
K = size(E.event_rate_time_by_mode, 2);
n_modes = min(K, params.max_plot_modes);
mode_idx = E.mode_index(1:n_modes);

fig = figure( ...
    'Color', params.background_color, ...
    'Position', params.figure_position, ...
    'Visible', params.figure_visible, ...
    'Name', params.title, ...
    'NumberTitle', 'off');

ax1 = subplot(3, 1, 1, 'Parent', fig);
imagesc(ax1, E.t_centers, 1:n_modes, ...
    double(E.smoothed_event_rate_time_by_mode(:, 1:n_modes)).');
local_style_heatmap_axis(ax1, params, mode_idx, n_modes);
ylabel(ax1, 'Mode', 'Color', params.text_color);
title(ax1, local_plot_title(E, params), 'Color', params.text_color, ...
    'Interpreter', 'none');
cb1 = colorbar(ax1);
cb1.Color = params.text_color;
cb1.Label.String = 'Smoothed event rate (events/s)';
cb1.Label.Color = params.text_color;
colormap(ax1, local_colormap(params.colormap));

ax2 = subplot(3, 1, 2, 'Parent', fig);
imagesc(ax2, E.t_centers, 1:n_modes, ...
    double(E.event_occupancy_time_by_mode(:, 1:n_modes)).');
local_style_heatmap_axis(ax2, params, mode_idx, n_modes);
ylabel(ax2, 'Mode', 'Color', params.text_color);
cb2 = colorbar(ax2);
cb2.Color = params.text_color;
cb2.Label.String = 'Event occupancy';
cb2.Label.Color = params.text_color;
colormap(ax2, local_colormap(params.colormap));
caxis(ax2, [0, 1]);

ax3 = subplot(3, 1, 3, 'Parent', fig);
set(ax3, ...
    'Color', params.axes_color, ...
    'XColor', params.text_color, ...
    'YColor', params.text_color, ...
    'GridColor', params.grid_color, ...
    'LineWidth', 0.8, ...
    'FontSize', 10);
hold(ax3, 'on');
grid(ax3, 'on');
box(ax3, 'off');

if ~isempty(E.event_table)
    rows = E.event_table.mode_col <= n_modes;
    Tplot = E.event_table(rows, :);
    if height(Tplot) > params.max_plot_events
        Tplot = Tplot(1:params.max_plot_events, :);
    end
    for i = 1:height(Tplot)
        y = Tplot.mode_col(i);
        plot(ax3, [Tplot.start_time(i), Tplot.end_time(i)], [y, y], ...
            'Color', [0.98, 0.86, 0.32], 'LineWidth', 1.0);
        plot(ax3, Tplot.peak_time(i), y, '.', ...
            'Color', [0.62, 0.85, 1.0], 'MarkerSize', 7);
    end
end

ylim(ax3, [0.5, n_modes + 0.5]);
set(ax3, 'YDir', 'reverse');
if n_modes <= 40
    yticks(ax3, 1:n_modes);
    yticklabels(ax3, compose('%d', mode_idx));
end
xlabel(ax3, 'Time (s)', 'Color', params.text_color);
ylabel(ax3, 'Mode', 'Color', params.text_color);
title(ax3, 'Accepted event intervals and peaks', ...
    'Color', params.text_color, 'Interpreter', 'none');

linkaxes([ax1, ax2, ax3], 'x');
end


function local_style_heatmap_axis(ax, params, mode_idx, n_modes)
set(ax, 'YDir', 'reverse');
set(ax, ...
    'Color', params.axes_color, ...
    'XColor', params.text_color, ...
    'YColor', params.text_color, ...
    'GridColor', params.grid_color, ...
    'LineWidth', 0.8, ...
    'FontSize', 10);
grid(ax, 'on');
box(ax, 'off');

if n_modes <= 40
    yticks(ax, 1:n_modes);
    yticklabels(ax, compose('%d', mode_idx));
end
end


function title_str = local_plot_title(E, params)
title_str = sprintf('%s | bin=%.3g s, nms=%s, frac=%.3g, score=%s', ...
    params.title, E.params.bin_sec_resolved, char(params.nms.mode), ...
    params.nms.period_fraction, char(params.event_score));
end


function cmap = local_colormap(name)
name = lower(char(name));
switch name
    case 'turbo'
        if exist('turbo', 'file') == 2 || exist('turbo', 'builtin') == 5
            cmap = turbo(256);
        else
            cmap = parula(256);
        end
    case 'parula'
        cmap = parula(256);
    case 'hot'
        cmap = hot(256);
    otherwise
        cmap = parula(256);
end
end


function value = local_get_field(S, field_name, default_value)
if isfield(S, field_name)
    value = S.(field_name);
else
    value = default_value;
end
end


function y = local_nanmean(A, dim)
valid = isfinite(A);
A(~valid) = 0;
counts = sum(valid, dim);
y = sum(A, dim) ./ max(counts, 1);
y(counts == 0) = NaN;
end


function y = local_nanmax(A, dim)
valid = isfinite(A);
A(~valid) = -Inf;
y = max(A, [], dim);
y(~any(valid, dim)) = NaN;
end


function y = local_nanmin_vector(v)
v = v(isfinite(v));
if isempty(v)
    y = NaN;
else
    y = min(v);
end
end


function y = local_nanmax_vector(v)
v = v(isfinite(v));
if isempty(v)
    y = NaN;
else
    y = max(v);
end
end
