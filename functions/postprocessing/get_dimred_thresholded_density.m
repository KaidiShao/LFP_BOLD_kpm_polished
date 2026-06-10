function [D, fig] = get_dimred_thresholded_density(dr_input, dt, params)
%GET_DIMRED_THRESHOLDED_DENSITY Thresholded density for reduced components.
%
%   [D, fig] = get_dimred_thresholded_density(result, dt, params)
%   [D, fig] = get_dimred_thresholded_density(C_time_by_comp, dt, params)
%
% This is a thin, metadata-preserving wrapper around get_thresholded_density
% for eigenfunction dimension-reduction outputs. It computes the same
% sample-wise above-threshold window density, but records the axis as
% time-by-component and saves component aliases in the MAT output.

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
[C_time_by_comp, source_meta, inferred_dt] = local_extract_components(dr_input);
params = local_fill_component_source_params(params, dr_input);

if isempty(dt)
    dt = inferred_dt;
end

component_index = local_component_index(params, size(C_time_by_comp, 2));
component_timescale_metadata = local_build_component_timescale_metadata( ...
    params, component_index, size(C_time_by_comp, 2), dt);
params = local_apply_component_envelope_windows( ...
    params, component_timescale_metadata);
component_timescale_metadata = local_add_component_envelope_columns( ...
    component_timescale_metadata, params);

core_params = params;
core_params.mode_index = component_index;
core_params.selected_mode_idx_in_original = [];
core_params.evalues_discrete = [];
core_params.evalues_continuous = [];
core_params.evalues_bilinear = [];
core_params.make_figure = params.make_figure || params.save_figure;
core_params.save_results = false;
core_params.save_figure = false;

[D, fig] = get_thresholded_density(C_time_by_comp, dt, core_params);

D.meta = source_meta;
D.meta.feature_family = 'eigenfunction_dimension_reduction';
D.meta.density_source = 'temporal_components_time_by_comp';

D.component_index = component_index(:);
D.density_time_by_component = D.density_time_by_mode;
D.valid_count_time_by_component = D.valid_count_time_by_mode;
D.threshold_by_component = D.threshold_by_mode;
D.component_timescale_metadata = local_update_component_envelope_from_activity( ...
    component_timescale_metadata, D.activity_metadata);

D.input.axis_order = 'time_by_component';
D.input.n_components = size(C_time_by_comp, 2);
D.input.component_source = 'result.core.temporal_components_time_by_comp';

D.params.save_results = params.save_results;
D.params.save_figure = params.save_figure;
D.params.make_figure = params.make_figure;
D.params.close_figure = params.close_figure;
D.params.component_index = component_index(:);

D.summary.n_components = size(C_time_by_comp, 2);
D.summary.density_size_time_by_component = size(D.density_time_by_component);

D.artifacts = local_ensure_artifact_paths(D.artifacts, params);

if params.save_figure
    if isempty(fig) || ~isvalid(fig)
        error('Figure saving was requested, but no valid figure was created.');
    end
    if exist(params.figure_dir, 'dir') ~= 7
        mkdir(params.figure_dir);
    end
    exportgraphics(fig, D.artifacts.figure_file, ...
        'Resolution', params.figure_resolution);
end

if params.save_results
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end

    density = D.density_time_by_component; %#ok<NASGU>
    density_time_by_component = D.density_time_by_component; %#ok<NASGU>
    thresholds = D.threshold_by_component; %#ok<NASGU>
    threshold_by_component = D.threshold_by_component; %#ok<NASGU>
    t_centers = D.t_centers; %#ok<NASGU>
    component_index = D.component_index; %#ok<NASGU>
    component_timescale_metadata = D.component_timescale_metadata; %#ok<NASGU>
    source_meta = D.meta; %#ok<NASGU>
    activity_metadata = D.activity_metadata; %#ok<NASGU>
    params_saved = D.params; %#ok<NASGU>

    if params.save_v7_3
        save(D.artifacts.mat_file, 'D', 'density', ...
            'density_time_by_component', 'thresholds', ...
            'threshold_by_component', 't_centers', ...
            'component_index', 'component_timescale_metadata', ...
            'source_meta', 'activity_metadata', 'params_saved', '-v7.3');
    else
        save(D.artifacts.mat_file, 'D', 'density', ...
            'density_time_by_component', 'thresholds', ...
            'threshold_by_component', 't_centers', ...
            'component_index', 'component_timescale_metadata', ...
            'source_meta', 'activity_metadata', 'params_saved');
    end
end

if params.close_figure && ~isempty(fig) && isvalid(fig)
    close(fig);
end
end


function params = local_apply_defaults(params)
if ~isfield(params, 'component_index')
    params.component_index = [];
end

if ~isfield(params, 'window_sec'), params.window_sec = 2; end
if ~isfield(params, 'step_sec'), params.step_sec = params.window_sec; end
if ~isfield(params, 'threshold_ratio'), params.threshold_ratio = 0.7; end
if ~isfield(params, 'threshold_mode'), params.threshold_mode = 'meanplusstd'; end
if ~isfield(params, 'value_transform'), params.value_transform = 'none'; end
if ~isfield(params, 'lfp_activity_transform') || isempty(params.lfp_activity_transform)
    params.lfp_activity_transform = 'abs_magnitude';
end
if ~isfield(params, 'lfp_activity_window_policy') || isempty(params.lfp_activity_window_policy)
    params.lfp_activity_window_policy = 'samplewise_abs_no_envelope';
end
if ~isfield(params, 'envelope_enable'), params.envelope_enable = false; end
if ~isfield(params, 'envelope_policy') || isempty(params.envelope_policy)
    params.envelope_policy = 'none';
end
if ~isfield(params, 'envelope_alpha') || isempty(params.envelope_alpha)
    params.envelope_alpha = 0.35;
end
if ~isfield(params, 'envelope_min_window_sec') || isempty(params.envelope_min_window_sec)
    params.envelope_min_window_sec = 0.03;
end
if ~isfield(params, 'envelope_max_window_sec') || isempty(params.envelope_max_window_sec)
    params.envelope_max_window_sec = 1.0;
end
if ~isfield(params, 'envelope_fallback_window_sec') || isempty(params.envelope_fallback_window_sec)
    params.envelope_fallback_window_sec = 0.10;
end
if ~isfield(params, 'envelope_window_sec'), params.envelope_window_sec = []; end
if ~isfield(params, 'envelope_window_sec_by_mode'), params.envelope_window_sec_by_mode = []; end
if ~isfield(params, 'envelope_window_samples'), params.envelope_window_samples = []; end
if ~isfield(params, 'envelope_window_samples_by_mode'), params.envelope_window_samples_by_mode = []; end
if ~isfield(params, 'activity_output_class') || isempty(params.activity_output_class)
    params.activity_output_class = 'single';
end
if ~isfield(params, 'source_evalues_discrete'), params.source_evalues_discrete = []; end
if ~isfield(params, 'source_evalues_continuous'), params.source_evalues_continuous = []; end
if ~isfield(params, 'source_evalues_bilinear'), params.source_evalues_bilinear = []; end
if ~isfield(params, 'mode_weights_mode_by_comp'), params.mode_weights_mode_by_comp = []; end
if ~isfield(params, 'selected_mode_idx_in_original')
    params.selected_mode_idx_in_original = [];
end
if ~isfield(params, 'component_timescale_weight_transform')
    params.component_timescale_weight_transform = 'abs';
end
if ~isfield(params, 'component_timescale_top_modes')
    params.component_timescale_top_modes = 5;
end
if ~isfield(params, 'density_denominator')
    params.density_denominator = 'window_samples';
end
if ~isfield(params, 'output_class'), params.output_class = 'single'; end
if ~isfield(params, 'smooth_density'), params.smooth_density = false; end
if ~isfield(params, 'smooth_window_sec'), params.smooth_window_sec = 0; end
if ~isfield(params, 'time_axis'), params.time_axis = []; end
if ~isfield(params, 'dt_source'), params.dt_source = []; end

if ~isfield(params, 'save_results'), params.save_results = false; end
if ~isfield(params, 'save_dir') || isempty(params.save_dir)
    params.save_dir = pwd;
end
if ~isfield(params, 'save_stem') || isempty(params.save_stem)
    params.save_stem = 'dimred_thresholded_density';
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
if ~isfield(params, 'title') || isempty(params.title)
    params.title = 'Thresholded Eigenfunction Component Density';
end
end


function [C, meta, dt] = local_extract_components(dr_input)
dt = [];
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

if isfield(dr_input, 'input') && isstruct(dr_input.input) && ...
        isfield(dr_input.input, 'dt') && ~isempty(dr_input.input.dt)
    dt = dr_input.input.dt;
end
end


function params = local_fill_component_source_params(params, dr_input)
if ~isstruct(dr_input)
    return;
end

if isempty(params.source_evalues_discrete) && isfield(dr_input, 'data') && ...
        isstruct(dr_input.data) && isfield(dr_input.data, 'evalues_discrete')
    params.source_evalues_discrete = dr_input.data.evalues_discrete;
end

if isempty(params.source_evalues_continuous) && isfield(dr_input, 'data') && ...
        isstruct(dr_input.data) && isfield(dr_input.data, 'evalues_bilinear')
    params.source_evalues_continuous = dr_input.data.evalues_bilinear;
end

if isempty(params.source_evalues_bilinear) && isfield(dr_input, 'data') && ...
        isstruct(dr_input.data) && isfield(dr_input.data, 'evalues_bilinear')
    params.source_evalues_bilinear = dr_input.data.evalues_bilinear;
end

if isempty(params.mode_weights_mode_by_comp) && isfield(dr_input, 'core') && ...
        isstruct(dr_input.core) && isfield(dr_input.core, 'mode_weights_mode_by_comp')
    params.mode_weights_mode_by_comp = dr_input.core.mode_weights_mode_by_comp;
end

if isempty(params.selected_mode_idx_in_original) && isfield(dr_input, 'input') && ...
        isstruct(dr_input.input) && isfield(dr_input.input, 'selected_mode_idx_in_original')
    params.selected_mode_idx_in_original = dr_input.input.selected_mode_idx_in_original;
end
end


function meta = local_build_component_timescale_metadata(params, component_index, C, dt)
density_index = (1:C).';
component_index = component_index(:);

source_evalues_discrete = local_complex_vector_or_nan( ...
    params.source_evalues_discrete, []);
expected_n = numel(source_evalues_discrete);
if expected_n == 0
    expected_n = [];
end
source_evalues_continuous = local_complex_vector_or_nan( ...
    params.source_evalues_continuous, expected_n);
if isempty(source_evalues_continuous)
    source_evalues_continuous = local_complex_vector_or_nan( ...
        params.source_evalues_bilinear, expected_n);
end

if ~isempty(source_evalues_discrete) && ...
        (isempty(source_evalues_continuous) || ...
        all(isnan(real(source_evalues_continuous)))) && ~isempty(dt)
    source_evalues_continuous = (2 / dt) .* ...
        (source_evalues_discrete - 1) ./ (source_evalues_discrete + 1);
end

mode_weights = params.mode_weights_mode_by_comp;
n_source_modes = size(mode_weights, 1);
if isempty(mode_weights) || size(mode_weights, 2) ~= C || ...
        isempty(source_evalues_continuous) || ...
        numel(source_evalues_continuous) ~= n_source_modes
    meta = local_empty_component_timescale_table( ...
        density_index, component_index, params);
    return;
end

selected_idx = params.selected_mode_idx_in_original(:);
if numel(selected_idx) ~= n_source_modes
    selected_idx = (1:n_source_modes).';
end

[source_timescale_discrete_log, source_frequency_discrete_angle] = ...
    local_timescale_from_discrete_lambda(source_evalues_discrete(:), dt);
[source_timescale_bilinear, source_frequency_bilinear] = ...
    local_timescale_from_continuous_lambda(source_evalues_continuous(:));
[source_timescale_sec, source_timescale_source] = ...
    local_prefer_discrete_timescale(source_timescale_discrete_log, source_timescale_bilinear);
[source_frequency_hz, source_frequency_source] = ...
    local_prefer_discrete_frequency(source_frequency_discrete_angle, source_frequency_bilinear);
W = local_transform_component_weights( ...
    mode_weights, params.component_timescale_weight_transform);
top_n = min(max(0, round(params.component_timescale_top_modes)), n_source_modes);

timescale_mean = nan(C, 1);
timescale_median = nan(C, 1);
frequency_mean = nan(C, 1);
n_valid_timescale_modes = zeros(C, 1);
top_source_mode_index = repmat({''}, C, 1);
top_raw_efun_index = repmat({''}, C, 1);
top_source_weight = repmat({''}, C, 1);

for c = 1:C
    wc = double(W(:, c));
    valid = isfinite(wc) & wc > 0 & isfinite(source_timescale_sec) & ...
        source_timescale_sec > 0;
    n_valid_timescale_modes(c) = nnz(valid);

    if any(valid)
        w = wc(valid);
        w = w ./ max(sum(w), eps);
        ts = source_timescale_sec(valid);
        fq = source_frequency_hz(valid);
        timescale_mean(c) = sum(w .* ts);
        timescale_median(c) = local_weighted_quantile(ts, w, 0.5);
        fq_valid = isfinite(fq);
        if any(fq_valid)
            wf = w(fq_valid);
            wf = wf ./ max(sum(wf), eps);
            frequency_mean(c) = sum(wf .* fq(fq_valid));
        end
    end

    if top_n > 0
        [sorted_w, ord] = sort(wc, 'descend');
        keep = ord(1:top_n);
        top_source_mode_index{c} = local_join_numeric(keep(:).');
        top_raw_efun_index{c} = local_join_numeric(selected_idx(keep).');
        top_source_weight{c} = local_join_numeric(sorted_w(1:top_n).');
    end
end

meta = table( ...
    density_index, ...
    component_index, ...
    repmat(n_source_modes, C, 1), ...
    n_valid_timescale_modes, ...
    timescale_mean, ...
    timescale_median, ...
    frequency_mean, ...
    top_source_mode_index, ...
    top_raw_efun_index, ...
    top_source_weight, ...
    repmat({source_timescale_source}, C, 1), ...
    repmat({source_frequency_source}, C, 1), ...
    repmat({char(params.component_timescale_weight_transform)}, C, 1), ...
    repmat({char(params.lfp_activity_transform)}, C, 1), ...
    repmat({char(params.lfp_activity_window_policy)}, C, 1), ...
    repmat(logical(params.envelope_enable), C, 1), ...
    repmat({char(params.envelope_policy)}, C, 1), ...
    'VariableNames', { ...
    'density_index', ...
    'component_index', ...
    'n_source_modes', ...
    'n_valid_timescale_modes', ...
    'weighted_timescale_sec_mean', ...
    'weighted_timescale_sec_median', ...
    'weighted_frequency_hz_mean', ...
    'top_source_mode_index', ...
    'top_raw_efun_index', ...
    'top_source_weight', ...
    'source_timescale_definition', ...
    'source_frequency_definition', ...
    'component_timescale_weight_transform', ...
    'lfp_activity_transform', ...
    'lfp_activity_window_policy', ...
    'envelope_enable', ...
    'envelope_policy'});
end


function meta = local_empty_component_timescale_table(density_index, component_index, params)
C = numel(component_index);
meta = table( ...
    density_index(:), ...
    component_index(:), ...
    zeros(C, 1), ...
    zeros(C, 1), ...
    nan(C, 1), ...
    nan(C, 1), ...
    nan(C, 1), ...
    repmat({''}, C, 1), ...
    repmat({''}, C, 1), ...
    repmat({''}, C, 1), ...
    repmat({'missing'}, C, 1), ...
    repmat({'missing'}, C, 1), ...
    repmat({char(params.component_timescale_weight_transform)}, C, 1), ...
    repmat({char(params.lfp_activity_transform)}, C, 1), ...
    repmat({char(params.lfp_activity_window_policy)}, C, 1), ...
    repmat(logical(params.envelope_enable), C, 1), ...
    repmat({char(params.envelope_policy)}, C, 1), ...
    'VariableNames', { ...
    'density_index', ...
    'component_index', ...
    'n_source_modes', ...
    'n_valid_timescale_modes', ...
    'weighted_timescale_sec_mean', ...
    'weighted_timescale_sec_median', ...
    'weighted_frequency_hz_mean', ...
    'top_source_mode_index', ...
    'top_raw_efun_index', ...
    'top_source_weight', ...
    'source_timescale_definition', ...
    'source_frequency_definition', ...
    'component_timescale_weight_transform', ...
    'lfp_activity_transform', ...
    'lfp_activity_window_policy', ...
    'envelope_enable', ...
    'envelope_policy'});
end


function params = local_apply_component_envelope_windows(params, meta)
needs_envelope = logical(params.envelope_enable) || ...
    contains(lower(char(string(params.lfp_activity_transform))), 'envelope') || ...
    contains(lower(char(string(params.envelope_policy))), 'adaptive');
if ~needs_envelope
    return;
end
if ~isempty(params.envelope_window_sec_by_mode) || ...
        ~isempty(params.envelope_window_samples_by_mode)
    return;
end
if isempty(meta) || ~istable(meta) || ...
        ~ismember('weighted_timescale_sec_median', meta.Properties.VariableNames)
    return;
end

tau = double(meta.weighted_timescale_sec_median(:));
alpha = double(params.envelope_alpha);
min_sec = double(params.envelope_min_window_sec);
max_sec = double(params.envelope_max_window_sec);
fallback_sec = double(params.envelope_fallback_window_sec);
window_sec = alpha .* tau;
invalid = ~isfinite(window_sec) | window_sec <= 0;
window_sec(invalid) = fallback_sec;
window_sec = max(min_sec, min(max_sec, window_sec));
params.envelope_window_sec_by_mode = window_sec;
end


function meta = local_add_component_envelope_columns(meta, params)
if isempty(meta) || ~istable(meta)
    return;
end
C = height(meta);
if ~isempty(params.envelope_window_sec_by_mode)
    win_sec = double(params.envelope_window_sec_by_mode(:));
    if numel(win_sec) ~= C
        win_sec = nan(C, 1);
    end
else
    win_sec = nan(C, 1);
end

if ~isempty(params.envelope_window_samples_by_mode)
    win_samples = double(params.envelope_window_samples_by_mode(:));
    if numel(win_samples) ~= C
        win_samples = zeros(C, 1);
    end
else
    win_samples = zeros(C, 1);
end

source = repmat({'component_weighted_timescale_adaptive'}, C, 1);
status = repmat({'ok'}, C, 1);
invalid = ~isfinite(win_sec) | win_sec <= 0;
status(invalid) = {'not_resolved'};

meta.envelope_window_sec = win_sec;
meta.envelope_window_samples = win_samples;
meta.envelope_window_source = source;
meta.envelope_window_status = status;
end


function meta = local_update_component_envelope_from_activity(meta, activity_meta)
if isempty(meta) || ~istable(meta) || ~isstruct(activity_meta)
    return;
end
C = height(meta);
if isfield(activity_meta, 'envelope_window_sec_by_mode')
    values = double(activity_meta.envelope_window_sec_by_mode(:));
    if numel(values) == C
        meta.envelope_window_sec = values;
    end
end
if isfield(activity_meta, 'envelope_window_samples_by_mode')
    values = double(activity_meta.envelope_window_samples_by_mode(:));
    if numel(values) == C
        meta.envelope_window_samples = values;
    end
end
if isfield(activity_meta, 'envelope_window_source_by_mode')
    values = cellstr(string(activity_meta.envelope_window_source_by_mode(:)));
    if numel(values) == C
        meta.envelope_window_source = values;
    end
end
if isfield(activity_meta, 'envelope_window_status_by_mode')
    values = cellstr(string(activity_meta.envelope_window_status_by_mode(:)));
    if numel(values) == C
        meta.envelope_window_status = values;
    end
end
end


function values = local_complex_vector_or_nan(values_in, expected_n)
if isempty(values_in)
    if isempty(expected_n)
        values = [];
    else
        values = complex(nan(expected_n, 1), nan(expected_n, 1));
    end
    return;
end

values = values_in(:);
if ~isempty(expected_n) && numel(values) ~= expected_n
    values = complex(nan(expected_n, 1), nan(expected_n, 1));
end
end


function W = local_transform_component_weights(mode_weights, transform_name)
switch lower(char(transform_name))
    case 'abs'
        W = abs(mode_weights);
    case {'abs2', 'squared_abs'}
        W = abs(mode_weights) .^ 2;
    case {'positive', 'positive_real'}
        W = max(real(mode_weights), 0);
    otherwise
        error('Unknown component_timescale_weight_transform = %s.', ...
            transform_name);
end
end


function [timescale_sec, frequency_hz] = local_timescale_from_discrete_lambda(lambda, dt)
lambda = lambda(:);
timescale_sec = nan(size(lambda));
frequency_hz = nan(size(lambda));
if isempty(lambda) || isempty(dt) || ~isfinite(dt) || dt <= 0
    return;
end

rho = abs(lambda);
theta = abs(angle(lambda));

valid_decay = isfinite(rho) & rho > 0 & rho < 1;
timescale_sec(valid_decay) = -double(dt) ./ log(rho(valid_decay));

unstable = isfinite(rho) & rho >= 1;
timescale_sec(unstable) = Inf;

valid_freq = isfinite(theta);
frequency_hz(valid_freq) = theta(valid_freq) ./ (2 * pi * double(dt));
end


function [preferred, source] = local_prefer_discrete_timescale(discrete_value, bilinear_value)
if isempty(discrete_value)
    preferred = nan(size(bilinear_value(:)));
else
    preferred = discrete_value(:);
end
source = repmat({'discrete_log_abs_evalue'}, numel(preferred), 1);

missing_discrete = isnan(preferred);
preferred(missing_discrete) = bilinear_value(missing_discrete);
source(missing_discrete) = {'bilinear_realpart'};

still_missing = isnan(preferred);
source(still_missing) = {'missing'};
source_text = unique(string(source));
source = char(strjoin(source_text(source_text ~= ""), '+'));
if isempty(source)
    source = 'missing';
end
end


function [preferred, source] = local_prefer_discrete_frequency(discrete_value, bilinear_value)
if isempty(discrete_value)
    preferred = nan(size(bilinear_value(:)));
else
    preferred = discrete_value(:);
end
source = repmat({'discrete_angle'}, numel(preferred), 1);

missing_discrete = isnan(preferred);
preferred(missing_discrete) = bilinear_value(missing_discrete);
source(missing_discrete) = {'bilinear_imag'};

still_missing = isnan(preferred);
source(still_missing) = {'missing'};
source_text = unique(string(source));
source = char(strjoin(source_text(source_text ~= ""), '+'));
if isempty(source)
    source = 'missing';
end
end


function [timescale_sec, frequency_hz] = local_timescale_from_continuous_lambda(lambda)
lambda = lambda(:);
lambda_real = real(lambda);
lambda_imag = imag(lambda);

timescale_sec = nan(size(lambda_real));
valid_decay = isfinite(lambda_real) & lambda_real < 0;
timescale_sec(valid_decay) = -1 ./ lambda_real(valid_decay);

frequency_hz = nan(size(lambda_imag));
valid_freq = isfinite(lambda_imag);
frequency_hz(valid_freq) = abs(lambda_imag(valid_freq)) ./ (2 * pi);
end


function qv = local_weighted_quantile(x, w, q)
[x_sorted, ord] = sort(double(x(:)));
w_sorted = double(w(ord));
w_sorted = w_sorted ./ max(sum(w_sorted), eps);
cw = cumsum(w_sorted);
idx = find(cw >= q, 1, 'first');
if isempty(idx)
    qv = x_sorted(end);
else
    qv = x_sorted(idx);
end
end


function text = local_join_numeric(values)
values = double(values(:));
parts = cell(numel(values), 1);
for i = 1:numel(values)
    parts{i} = sprintf('%.10g', values(i));
end
text = strjoin(parts, ';');
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
        local_ratio_tag(params.threshold_ratio)};
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
