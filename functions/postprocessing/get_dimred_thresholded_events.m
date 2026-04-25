function [E, figs] = get_dimred_thresholded_events(dr_input, dt, params)
%GET_DIMRED_THRESHOLDED_EVENTS Thresholded events for reduced components.
%
%   [E, figs] = get_dimred_thresholded_events(result, dt, params)
%   [E, figs] = get_dimred_thresholded_events(C_time_by_comp, dt, params)
%
% This wrapper applies get_thresholded_events to eigenfunction
% dimension-reduction temporal components. When a full reduction result is
% supplied, it estimates one effective discrete eigenvalue per component
% from the component's source-mode weights so eigen-period NMS can still be
% used at the component level.

if nargin < 2
    dt = [];
end

if nargin < 3 || isempty(params)
    if isstruct(dt)
        params = dt;
        dt = [];
    else
        params = struct();
    end
end

params = local_apply_defaults(params);
[C_time_by_comp, source_meta, inferred_dt, source_evalues, mode_weights] = ...
    local_extract_components(dr_input);

if isempty(dt)
    dt = inferred_dt;
end

if isempty(dt) || ~isnumeric(dt) || ~isscalar(dt) || ~isfinite(dt) || dt <= 0
    error('A positive scalar dt is required for component event detection.');
end

n_components = size(C_time_by_comp, 2);
component_index = local_component_index(params, n_components);
[component_evalues, component_evalue_info] = local_component_evalues( ...
    params, source_evalues, mode_weights, n_components);

core_params = params;
core_params.mode_index = component_index;
core_params.selected_mode_idx_in_original = [];
core_params.make_figure = params.make_figure || params.save_figure;
core_params.save_results = false;
core_params.save_figure = false;

[E, figs] = get_thresholded_events( ...
    C_time_by_comp, dt, component_evalues, core_params);

E.meta = source_meta;
E.meta.feature_family = 'eigenfunction_dimension_reduction';
E.meta.event_source = 'temporal_components_time_by_comp';

E.component_index = component_index(:);
E.component_evalues_discrete = component_evalues(:);
E.component_evalue_info = component_evalue_info;
E.component_timescales = E.mode_timescales;

E.event_rate_time_by_component = E.event_rate_time_by_mode;
E.smoothed_event_rate_time_by_component = E.smoothed_event_rate_time_by_mode;
E.event_count_time_by_component = E.event_count_time_by_mode;
E.event_occupancy_time_by_component = E.event_occupancy_time_by_mode;
E.event_presence_time_by_component = E.event_presence_time_by_mode;
E.threshold_by_component = E.threshold_by_mode;

E.input.axis_order = 'time_by_component';
E.input.n_components = n_components;
E.input.component_source = 'result.core.temporal_components_time_by_comp';

E.params.save_results = params.save_results;
E.params.save_figure = params.save_figure;
E.params.make_figure = params.make_figure;
E.params.close_figure = params.close_figure;
E.params.component_index = component_index(:);
E.params.component_evalue_weight_transform = ...
    params.component_evalue_weight_transform;
E.params.component_evalue_angle_statistic = ...
    params.component_evalue_angle_statistic;

E.summary.n_components = n_components;
E.summary.event_rate_size_time_by_component = size(E.event_rate_time_by_component);

E.artifacts = local_ensure_artifact_paths(E.artifacts, params);

if params.save_figure
    if ~isfield(figs, 'summary') || isempty(figs.summary) || ...
            ~isvalid(figs.summary)
        error('Figure saving was requested, but no valid figure was created.');
    end
    if exist(params.figure_dir, 'dir') ~= 7
        mkdir(params.figure_dir);
    end
    exportgraphics(figs.summary, E.artifacts.figure_file, ...
        'Resolution', params.figure_resolution);
end

if params.save_results
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end

    event_table = E.event_table; %#ok<NASGU>
    event_rate = E.event_rate_time_by_component; %#ok<NASGU>
    event_rate_time_by_component = E.event_rate_time_by_component; %#ok<NASGU>
    smoothed_event_rate = E.smoothed_event_rate_time_by_component; %#ok<NASGU>
    event_occupancy = E.event_occupancy_time_by_component; %#ok<NASGU>
    event_presence = E.event_presence_time_by_component; %#ok<NASGU>
    thresholds = E.threshold_by_component; %#ok<NASGU>
    threshold_by_component = E.threshold_by_component; %#ok<NASGU>
    component_timescales = E.component_timescales; %#ok<NASGU>
    component_evalues = E.component_evalues_discrete; %#ok<NASGU>
    component_evalue_info = E.component_evalue_info; %#ok<NASGU>
    t_centers = E.t_centers; %#ok<NASGU>
    component_index = E.component_index; %#ok<NASGU>
    source_meta = E.meta; %#ok<NASGU>
    params_saved = E.params; %#ok<NASGU>

    if params.save_v7_3
        save(E.artifacts.mat_file, 'E', 'event_table', 'event_rate', ...
            'event_rate_time_by_component', 'smoothed_event_rate', ...
            'event_occupancy', 'event_presence', 'thresholds', ...
            'threshold_by_component', 'component_timescales', ...
            'component_evalues', 'component_evalue_info', 't_centers', ...
            'component_index', 'source_meta', 'params_saved', '-v7.3');
    else
        save(E.artifacts.mat_file, 'E', 'event_table', 'event_rate', ...
            'event_rate_time_by_component', 'smoothed_event_rate', ...
            'event_occupancy', 'event_presence', 'thresholds', ...
            'threshold_by_component', 'component_timescales', ...
            'component_evalues', 'component_evalue_info', 't_centers', ...
            'component_index', 'source_meta', 'params_saved');
    end
end

if params.close_figure && isfield(figs, 'summary') && ...
        ~isempty(figs.summary) && isvalid(figs.summary)
    close(figs.summary);
end
end


function params = local_apply_defaults(params)
if ~isfield(params, 'component_index')
    params.component_index = [];
end
if ~isfield(params, 'component_evalues_discrete')
    params.component_evalues_discrete = [];
end
if ~isfield(params, 'source_evalues_discrete')
    params.source_evalues_discrete = [];
end
if ~isfield(params, 'mode_weights_mode_by_comp')
    params.mode_weights_mode_by_comp = [];
end
if ~isfield(params, 'component_evalue_weight_transform') || ...
        isempty(params.component_evalue_weight_transform)
    params.component_evalue_weight_transform = 'abs';
end
if ~isfield(params, 'component_evalue_angle_statistic') || ...
        isempty(params.component_evalue_angle_statistic)
    params.component_evalue_angle_statistic = 'weighted_median_abs';
end
if ~isfield(params, 'component_evalue_top_modes') || ...
        isempty(params.component_evalue_top_modes)
    params.component_evalue_top_modes = 5;
end

if ~isfield(params, 'threshold_ratio'), params.threshold_ratio = 0.7; end
if ~isfield(params, 'threshold_mode'), params.threshold_mode = 'meanplusstd'; end
if ~isfield(params, 'value_transform'), params.value_transform = 'none'; end
if ~isfield(params, 'event_score'), params.event_score = 'area_above_threshold'; end
if ~isfield(params, 'min_duration_sec'), params.min_duration_sec = 0.03; end
if ~isfield(params, 'merge_gap_sec'), params.merge_gap_sec = 0.02; end
if ~isfield(params, 'bin_sec'), params.bin_sec = 2; end
if ~isfield(params, 'step_sec'), params.step_sec = params.bin_sec; end
if ~isfield(params, 'smooth_rate'), params.smooth_rate = true; end
if ~isfield(params, 'smooth_sigma_sec'), params.smooth_sigma_sec = params.bin_sec; end
if ~isfield(params, 'output_class'), params.output_class = 'single'; end
if ~isfield(params, 'time_axis'), params.time_axis = []; end
if ~isfield(params, 'dt_source'), params.dt_source = []; end
if ~isfield(params, 'observable_file'), params.observable_file = ''; end
if ~isfield(params, 'require_session_metadata')
    params.require_session_metadata = false;
end

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
    params.save_stem = 'dimred_thresholded_events';
end
if ~isfield(params, 'save_tag'), params.save_tag = ''; end
if ~isfield(params, 'save_v7_3'), params.save_v7_3 = true; end

if ~isfield(params, 'make_figure'), params.make_figure = false; end
if ~isfield(params, 'save_figure'), params.save_figure = false; end
if ~isfield(params, 'close_figure'), params.close_figure = true; end
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
if ~isfield(params, 'figure_resolution'), params.figure_resolution = 200; end
if ~isfield(params, 'max_plot_modes') || isempty(params.max_plot_modes)
    params.max_plot_modes = 40;
end
if ~isfield(params, 'max_plot_events') || isempty(params.max_plot_events)
    params.max_plot_events = 5000;
end
if ~isfield(params, 'title') || isempty(params.title)
    params.title = 'Thresholded Eigenfunction Component Events';
end
end


function [C, meta, dt, evalues, weights] = local_extract_components(dr_input)
dt = [];
evalues = [];
weights = [];
meta = struct();
meta.source_type = '';
meta.path_kind = '';
meta.method = '';
meta.component_domain = '';
meta.reconstruction_domain = '';

if isnumeric(dr_input)
    C = dr_input;
    meta.source_type = 'matrix';
    return;
end

if ~isstruct(dr_input)
    error('dr_input must be a result struct or a numeric [T x C] matrix.');
end

if ~isfield(dr_input, 'core') || ~isstruct(dr_input.core) || ...
        ~isfield(dr_input.core, 'temporal_components_time_by_comp')
    error(['Result input must contain ', ...
        'core.temporal_components_time_by_comp.']);
end

C = dr_input.core.temporal_components_time_by_comp;
meta.source_type = 'eigenfunction_reduction_result';

if isfield(dr_input, 'meta') && isstruct(dr_input.meta)
    meta.path_kind = local_get_field(dr_input.meta, 'path_kind', '');
    meta.method = local_get_field(dr_input.meta, 'method', '');
end

if isfield(dr_input.core, 'component_domain')
    meta.component_domain = dr_input.core.component_domain;
end
if isfield(dr_input.core, 'reconstruction_domain')
    meta.reconstruction_domain = dr_input.core.reconstruction_domain;
end
if isfield(dr_input.core, 'mode_weights_mode_by_comp')
    weights = dr_input.core.mode_weights_mode_by_comp;
end
if isfield(dr_input, 'data') && isstruct(dr_input.data) && ...
        isfield(dr_input.data, 'evalues_discrete')
    evalues = dr_input.data.evalues_discrete;
end
if isfield(dr_input, 'input') && isstruct(dr_input.input) && ...
        isfield(dr_input.input, 'dt') && ~isempty(dr_input.input.dt)
    dt = dr_input.input.dt;
end
end


function component_index = local_component_index(params, C)
if isfield(params, 'component_index') && ~isempty(params.component_index)
    component_index = params.component_index(:);
    if numel(component_index) ~= C
        error('params.component_index must contain one entry per component.');
    end
else
    component_index = (1:C).';
end
end


function [component_evalues, info] = local_component_evalues( ...
        params, source_evalues, mode_weights, C)
component_evalues = nan(C, 1);

info = struct();
info.source = 'missing';
info.method = 'weighted_effective_discrete_evalue';
info.weight_transform = params.component_evalue_weight_transform;
info.angle_statistic = params.component_evalue_angle_statistic;
info.message = '';
info.top_source_mode_idx_by_comp = [];
info.top_source_weight_by_comp = [];
info.top_source_evalue_by_comp = [];

if ~isempty(params.component_evalues_discrete)
    component_evalues = params.component_evalues_discrete(:);
    if numel(component_evalues) ~= C
        error('params.component_evalues_discrete must have one value per component.');
    end
    info.source = 'params.component_evalues_discrete';
    info.message = 'Using caller-supplied component eigenvalues.';
    return;
end

if ~isempty(params.source_evalues_discrete)
    source_evalues = params.source_evalues_discrete;
end
if ~isempty(params.mode_weights_mode_by_comp)
    mode_weights = params.mode_weights_mode_by_comp;
end

if isempty(source_evalues) || isempty(mode_weights)
    info.message = ['No source eigenvalues/mode weights were available; ', ...
        'component event NMS will use get_thresholded_events defaults.'];
    return;
end

source_evalues = source_evalues(:);
if size(mode_weights, 1) ~= numel(source_evalues) || size(mode_weights, 2) ~= C
    error(['mode_weights_mode_by_comp must be [n_source_modes x n_components] ', ...
        'and match source_evalues_discrete.']);
end

W = local_transform_weights(mode_weights, params.component_evalue_weight_transform);
info.source = 'source_evalues_weighted_by_component_loadings';

top_n = min(max(0, round(params.component_evalue_top_modes)), size(W, 1));
if top_n > 0
    info.top_source_mode_idx_by_comp = nan(top_n, C);
    info.top_source_weight_by_comp = nan(top_n, C);
    info.top_source_evalue_by_comp = nan(top_n, C);
end

for c = 1:C
    wc = double(W(:, c));
    valid = isfinite(wc) & wc > 0 & isfinite(source_evalues);
    if ~any(valid)
        continue;
    end

    lambda = source_evalues(valid);
    w = wc(valid);
    w = w ./ sum(w);

    radius = exp(sum(w .* log(max(abs(lambda), realmin))));
    angle_abs = abs(angle(lambda));
    angle_eff = local_effective_angle(angle_abs, w, ...
        params.component_evalue_angle_statistic);

    component_evalues(c) = radius .* exp(1i * angle_eff);

    if top_n > 0
        [sorted_w, idx] = sort(wc, 'descend');
        idx = idx(1:top_n);
        info.top_source_mode_idx_by_comp(:, c) = idx(:);
        info.top_source_weight_by_comp(:, c) = sorted_w(1:top_n);
        info.top_source_evalue_by_comp(:, c) = source_evalues(idx);
    end
end
end


function W = local_transform_weights(mode_weights, transform_name)
switch lower(char(transform_name))
    case 'abs'
        W = abs(mode_weights);
    case {'abs2', 'squared_abs'}
        W = abs(mode_weights) .^ 2;
    case {'positive', 'positive_real'}
        W = max(real(mode_weights), 0);
    otherwise
        error(['Unknown component_evalue_weight_transform = %s. ', ...
            'Use abs, abs2, or positive.'], transform_name);
end
end


function angle_eff = local_effective_angle(angle_abs, w, statistic_name)
angle_abs = double(angle_abs(:));
w = double(w(:));

switch lower(char(statistic_name))
    case {'weighted_median_abs', 'median_abs'}
        angle_eff = local_weighted_quantile(angle_abs, w, 0.5);
    case {'weighted_mean_abs', 'mean_abs'}
        angle_eff = sum(w .* angle_abs) ./ max(sum(w), eps);
    otherwise
        error(['Unknown component_evalue_angle_statistic = %s. ', ...
            'Use weighted_median_abs or weighted_mean_abs.'], statistic_name);
end
end


function qv = local_weighted_quantile(x, w, q)
[x_sorted, ord] = sort(x(:));
w_sorted = w(ord);
w_sorted = w_sorted ./ max(sum(w_sorted), eps);
cw = cumsum(w_sorted);
idx = find(cw >= q, 1, 'first');
if isempty(idx)
    qv = x_sorted(end);
else
    qv = x_sorted(idx);
end
end


function value = local_get_field(S, name, default_value)
if isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function artifacts = local_ensure_artifact_paths(artifacts, params)
if ~isfield(artifacts, 'file_stem') || isempty(artifacts.file_stem)
    timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
    pieces = {char(params.save_stem), char(params.threshold_mode), ...
        local_ratio_tag(params.threshold_ratio), 'events'};
    if isfield(params, 'save_tag') && ~isempty(params.save_tag)
        pieces{end+1} = char(params.save_tag);
    end
    base = strjoin(pieces, '__');
    base = regexprep(base, '[^\w\-]+', '_');
    artifacts.file_stem = sprintf('%s__%s', base, timestamp);
end

if params.save_results && ...
        (~isfield(artifacts, 'mat_file') || isempty(artifacts.mat_file))
    artifacts.mat_file = fullfile(params.save_dir, ...
        [artifacts.file_stem, '.mat']);
end

if params.save_figure && ...
        (~isfield(artifacts, 'figure_file') || isempty(artifacts.figure_file))
    artifacts.figure_file = fullfile(params.figure_dir, ...
        [artifacts.file_stem, '.png']);
end
end


function tag = local_ratio_tag(ratio)
if isscalar(ratio)
    tag = sprintf('ratio_%03d', round(double(ratio) * 100));
else
    tag = sprintf('ratio_vector%d', numel(ratio));
end
end
