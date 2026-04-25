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

if isempty(dt)
    dt = inferred_dt;
end

component_index = local_component_index(params, size(C_time_by_comp, 2));

core_params = params;
core_params.mode_index = component_index;
core_params.selected_mode_idx_in_original = [];
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
    source_meta = D.meta; %#ok<NASGU>
    params_saved = D.params; %#ok<NASGU>

    if params.save_v7_3
        save(D.artifacts.mat_file, 'D', 'density', ...
            'density_time_by_component', 'thresholds', ...
            'threshold_by_component', 't_centers', ...
            'component_index', 'source_meta', 'params_saved', '-v7.3');
    else
        save(D.artifacts.mat_file, 'D', 'density', ...
            'density_time_by_component', 'thresholds', ...
            'threshold_by_component', 't_centers', ...
            'component_index', 'source_meta', 'params_saved');
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
