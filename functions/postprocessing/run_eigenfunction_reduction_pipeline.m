function [result, EDMD_outputs, concat_info, source_info] = run_eigenfunction_reduction_pipeline(cfg)
%RUN_EIGENFUNCTION_REDUCTION_PIPELINE Minimal eigenfunction-only reduction pipeline.

if nargin < 1 || isempty(cfg)
    error('cfg must be provided.');
end

cfg = local_apply_defaults(cfg);
progress_enabled = local_progress_enabled(cfg);
pipeline_tic = tic;

stage_tic = local_stage_start(progress_enabled, 'Loading EDMD source');
[EDMD_outputs, concat_info, source_info] = load_edmd_source(cfg.source);
local_stage_done(progress_enabled, 'Loaded EDMD source', stage_tic);

stage_tic = local_stage_start(progress_enabled, 'Preparing eigenfunction inputs');
prep = local_prepare_eigenfunction_inputs(EDMD_outputs, cfg);
local_stage_done(progress_enabled, 'Prepared eigenfunction inputs', stage_tic);
local_print_prepared_input_summary(progress_enabled, prep, EDMD_outputs);

switch lower(cfg.path.kind)
    case 'time'
        method_name = cfg.path.time.method;
        stage_tic = local_stage_start(progress_enabled, ...
            sprintf('Running dimension reduction: %s/%s', cfg.path.kind, method_name));
        reduction = reduce_eigenfunction_time_path( ...
            prep.efun_feature_time_by_mode, cfg.path.time);
        local_stage_done(progress_enabled, 'Finished dimension reduction', stage_tic);

    case 'spectrum'
        spectrum_cfg = cfg.path.spectrum;
        spectrum_cfg.evalues_discrete = prep.evalues_discrete;
        spectrum_cfg.evalues_bilinear = prep.evalues_bilinear;
        spectrum_cfg.feature_variant = cfg.feature.variant;
        method_name = cfg.path.spectrum.method;
        stage_tic = local_stage_start(progress_enabled, ...
            sprintf('Running dimension reduction: %s/%s', cfg.path.kind, method_name));
        reduction = reduce_eigenfunction_spectrum_path( ...
            prep.efun_feature_time_by_mode, spectrum_cfg);
        local_stage_done(progress_enabled, 'Finished dimension reduction', stage_tic);

    otherwise
        error('Unknown cfg.path.kind = %s. Use ''time'' or ''spectrum''.', cfg.path.kind);
end

result = struct();

result.meta = struct();
result.meta.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
result.meta.path_kind = cfg.path.kind;
result.meta.feature_family = 'eigenfunction';
result.meta.feature_variant = cfg.feature.variant;
result.meta.method = method_name;

result.cfg = local_strip_preloaded_source_payload(cfg);
result.source = source_info;
result.concat = concat_info;

result.input = struct();
result.input.dt = prep.dt;
result.input.dt_source = prep.dt_source;
result.input.time_axis = prep.time_axis;
result.input.mode_index = (1:numel(prep.evalues_discrete)).';
result.input.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;
result.input.selected_mode_mask_in_original = prep.selected_mode_mask_in_original;

result.data = struct();
result.data.evalues_discrete = prep.evalues_discrete;
result.data.evalues_bilinear = prep.evalues_bilinear;
result.data.efun_raw_time_by_mode = prep.efun_raw_time_by_mode;
result.data.efun_feature_time_by_mode = prep.efun_feature_time_by_mode;
result.data.kpm_modes_mode_by_dict = prep.kpm_modes_mode_by_dict;

result.feature = struct();
result.feature.family = 'eigenfunction';
result.feature.variant = cfg.feature.variant;
result.feature.normalization = cfg.feature.normalization;
result.feature.axis_order = 'time_by_mode';

result.core = reduction.core;
result.quality = reduction.quality;
result.aux = reduction.aux;

stage_tic = local_stage_start(progress_enabled, 'Building smoothed component summary');
result.summary = local_build_summary(result.core, cfg.summary);
local_stage_done(progress_enabled, 'Built smoothed component summary', stage_tic);

result.artifacts = struct();
result.artifacts.result_mat_file = '';
result.artifacts.save_payload = cfg.save.payload;
result.artifacts.thresholded_density_mat_file = '';
result.artifacts.thresholded_density_figure_file = '';
result.artifacts.thresholded_events_mat_file = '';
result.artifacts.thresholded_events_figure_file = '';
result.artifacts.dimred_thresholded_density_mat_file = '';
result.artifacts.dimred_thresholded_density_figure_file = '';
result.artifacts.dimred_thresholded_events_mat_file = '';
result.artifacts.dimred_thresholded_events_figure_file = '';

result.thresholded_density = struct();
if cfg.thresholded_density.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded eigenfunction density');
    [result.thresholded_density, td_artifacts] = local_run_thresholded_density( ...
        prep, result, cfg.thresholded_density);
    result.artifacts.thresholded_density_mat_file = td_artifacts.mat_file;
    result.artifacts.thresholded_density_figure_file = td_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded eigenfunction density', stage_tic);
end

result.thresholded_events = struct();
if cfg.thresholded_events.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded eigenfunction events');
    [result.thresholded_events, te_artifacts] = local_run_thresholded_events( ...
        prep, result, cfg.thresholded_events);
    result.artifacts.thresholded_events_mat_file = te_artifacts.mat_file;
    result.artifacts.thresholded_events_figure_file = te_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded eigenfunction events', stage_tic);
end

result.dimred_thresholded_density = struct();
if cfg.dimred_thresholded_density.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded component density');
    [result.dimred_thresholded_density, dtd_artifacts] = ...
        local_run_dimred_thresholded_density( ...
        prep, result, cfg.dimred_thresholded_density);
    result.artifacts.dimred_thresholded_density_mat_file = dtd_artifacts.mat_file;
    result.artifacts.dimred_thresholded_density_figure_file = dtd_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded component density', stage_tic);
end

result.dimred_thresholded_events = struct();
if cfg.dimred_thresholded_events.enable
    stage_tic = local_stage_start(progress_enabled, 'Running thresholded component events');
    [result.dimred_thresholded_events, dte_artifacts] = ...
        local_run_dimred_thresholded_events( ...
        prep, result, cfg.dimred_thresholded_events);
    result.artifacts.dimred_thresholded_events_mat_file = dte_artifacts.mat_file;
    result.artifacts.dimred_thresholded_events_figure_file = dte_artifacts.figure_file;
    local_stage_done(progress_enabled, 'Finished thresholded component events', stage_tic);
end

if cfg.save.enable
    stage_tic = local_stage_start(progress_enabled, 'Saving compact reduction result');
    if exist(cfg.save.dir, 'dir') ~= 7
        mkdir(cfg.save.dir);
    end

    save_path = local_build_save_path(cfg, result);
    result.artifacts.result_mat_file = save_path;
    result_to_save = local_build_save_result(result, cfg.save.payload);
    save_vars = struct('result', result_to_save);

    if cfg.save.v7_3
        save(save_path, '-struct', 'save_vars', 'result', '-v7.3');
    else
        save(save_path, '-struct', 'save_vars', 'result');
    end
    local_stage_done(progress_enabled, 'Saved reduction result', stage_tic);
end

local_stage_done(progress_enabled, 'Pipeline method finished', pipeline_tic);
end


function tf = local_progress_enabled(cfg)
tf = true;
if isfield(cfg, 'progress') && isfield(cfg.progress, 'verbose') && ...
        ~isempty(cfg.progress.verbose)
    tf = logical(cfg.progress.verbose);
elseif isfield(cfg, 'source') && isfield(cfg.source, 'concat') && ...
        isfield(cfg.source.concat, 'verbose') && ~isempty(cfg.source.concat.verbose)
    tf = logical(cfg.source.concat.verbose);
end
end


function stage_tic = local_stage_start(enabled, label)
stage_tic = tic;
if enabled
    fprintf('[stage] %s...\n', label);
end
end


function local_stage_done(enabled, label, stage_tic)
if enabled
    fprintf('[stage] %s in %s.\n', label, local_format_seconds(toc(stage_tic)));
end
end


function local_print_prepared_input_summary(enabled, prep, EDMD_outputs)
if ~enabled
    return;
end

raw_size = size(prep.efun_raw_time_by_mode);
feature_size = size(prep.efun_feature_time_by_mode);
n_original_modes = numel(EDMD_outputs.evalues);
n_selected_modes = numel(prep.evalues_discrete);

if isempty(prep.dt)
    dt_text = 'missing';
    duration_text = 'unknown';
else
    dt_text = sprintf('%.10g sec', prep.dt);
    duration_text = local_format_duration_for_summary( ...
        (raw_size(1) - 1) * double(prep.dt));
end

fprintf(['[summary] Eigenfunction length: T=%d samples | selected modes=%d/%d | ', ...
    'raw=[%s] | feature=[%s]\n'], ...
    raw_size(1), n_selected_modes, n_original_modes, ...
    num2str(raw_size), num2str(feature_size));
fprintf('[summary] Time axis: dt=%s | duration=%s | dt_source=%s.%s\n', ...
    dt_text, duration_text, prep.dt_source.source, prep.dt_source.field);
end


function txt = local_format_duration_for_summary(seconds_in)
if ~isfinite(seconds_in)
    txt = 'unknown';
    return;
end

seconds_in = max(0, seconds_in);
days = floor(seconds_in / 86400);
hours = floor(mod(seconds_in, 86400) / 3600);
minutes = floor(mod(seconds_in, 3600) / 60);
seconds = floor(mod(seconds_in, 60));

if days > 0
    txt = sprintf('%dd %02d:%02d:%02d', days, hours, minutes, seconds);
elseif hours > 0
    txt = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
else
    txt = sprintf('%02d:%02d', minutes, seconds);
end
end


function txt = local_format_seconds(seconds_in)
if ~isfinite(seconds_in)
    txt = '--';
    return;
end

seconds_in = max(0, seconds_in);
if seconds_in < 1
    txt = sprintf('%.3fs', seconds_in);
elseif seconds_in < 60
    txt = sprintf('%.1fs', seconds_in);
else
    hours = floor(seconds_in / 3600);
    minutes = floor(mod(seconds_in, 3600) / 60);
    seconds = floor(mod(seconds_in, 60));
    if hours > 0
        txt = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
    else
        txt = sprintf('%02d:%02d', minutes, seconds);
    end
end
end


function cfg = local_apply_defaults(cfg)
if ~isfield(cfg, 'feature') || ~isfield(cfg.feature, 'variant')
    error('cfg.feature.variant must be provided.');
end

if ~isfield(cfg, 'feature') || ~isfield(cfg.feature, 'normalization') || isempty(cfg.feature.normalization)
    cfg.feature.normalization = 'maxabs_per_mode';
end

if ~isfield(cfg, 'input') || ~isfield(cfg.input, 'dt')
    cfg.input.dt = [];
end

if ~isfield(cfg.input, 'observable_file')
    cfg.input.observable_file = '';
end

if ~isfield(cfg, 'selection') || ~isfield(cfg.selection, 'abs_thresh')
    cfg.selection.abs_thresh = 0.01;
end

if ~isfield(cfg.selection, 'sort_by') || isempty(cfg.selection.sort_by)
    cfg.selection.sort_by = 'modulus';
end

if ~isfield(cfg.selection, 'sort_dir') || isempty(cfg.selection.sort_dir)
    cfg.selection.sort_dir = 'descend';
end

if ~isfield(cfg.selection, 'max_modes') || isempty(cfg.selection.max_modes)
    cfg.selection.max_modes = Inf;
end

if ~isfield(cfg, 'summary') || ~isfield(cfg.summary, 'smooth')
    cfg.summary.smooth = struct();
end

if ~isfield(cfg.summary.smooth, 'enable') || isempty(cfg.summary.smooth.enable)
    cfg.summary.smooth.enable = false;
end

if ~isfield(cfg.summary.smooth, 'method') || isempty(cfg.summary.smooth.method)
    cfg.summary.smooth.method = 'movmean';
end

if ~isfield(cfg.summary.smooth, 'window') || isempty(cfg.summary.smooth.window)
    cfg.summary.smooth.window = 10;
end

if ~isfield(cfg, 'save') || ~isfield(cfg.save, 'enable')
    cfg.save.enable = false;
end

if ~isfield(cfg.save, 'dir') || isempty(cfg.save.dir)
    cfg.save.dir = local_default_save_dir(cfg);
end

if ~isfield(cfg.save, 'file_stem') || isempty(cfg.save.file_stem)
    cfg.save.file_stem = local_default_file_stem(cfg);
end

if ~isfield(cfg.save, 'tag')
    cfg.save.tag = '';
end

if ~isfield(cfg.save, 'payload') || isempty(cfg.save.payload)
    cfg.save.payload = 'compact';
end

if ~isfield(cfg.save, 'v7_3') || isempty(cfg.save.v7_3)
    cfg.save.v7_3 = true;
end

cfg = local_apply_thresholded_density_defaults(cfg);
cfg = local_apply_thresholded_events_defaults(cfg);
cfg = local_apply_dimred_thresholded_density_defaults(cfg);
cfg = local_apply_dimred_thresholded_events_defaults(cfg);
end


function prep = local_prepare_eigenfunction_inputs(EDMD_outputs, cfg)
required_fields = {'evalues', 'efuns'};
for i = 1:numel(required_fields)
    if ~isfield(EDMD_outputs, required_fields{i})
        error('EDMD_outputs.%s is required.', required_fields{i});
    end
end

evalues0 = EDMD_outputs.evalues(:);
efuns0 = EDMD_outputs.efuns;

if size(efuns0, 2) ~= numel(evalues0)
    error('size(EDMD_outputs.efuns, 2) must equal numel(EDMD_outputs.evalues).');
end

ord = local_sort_eigenvalues(evalues0, cfg.selection.sort_by, cfg.selection.sort_dir);
evalues_sorted = evalues0(ord);
mask_sorted = abs(evalues_sorted) > cfg.selection.abs_thresh;
idx_sorted_selected = ord(mask_sorted);

if isempty(idx_sorted_selected)
    error('No modes remain after applying abs threshold %g.', cfg.selection.abs_thresh);
end

if isfinite(cfg.selection.max_modes)
    max_modes = min(numel(idx_sorted_selected), cfg.selection.max_modes);
    idx_sorted_selected = idx_sorted_selected(1:max_modes);
end

selected_mode_mask = false(size(evalues0));
selected_mode_mask(idx_sorted_selected) = true;

efun_raw = efuns0(:, idx_sorted_selected);
evalues_selected = evalues0(idx_sorted_selected);

switch lower(cfg.feature.normalization)
    case 'maxabs_per_mode'
        efun_feature = normalize_efun(efun_raw, cfg.feature.variant);
    otherwise
        error('Unsupported cfg.feature.normalization = %s.', cfg.feature.normalization);
end

[dt, dt_source] = local_resolve_dt(cfg.input, EDMD_outputs);
if isempty(dt)
    time_axis = (1:size(efun_raw, 1)).';
    evalues_bilinear = [];
else
    time_axis = (0:size(efun_raw, 1)-1).' * dt;
    evalues_bilinear = (2 / dt) * (evalues_selected - 1) ./ (evalues_selected + 1);
end

if isfield(EDMD_outputs, 'kpm_modes') && size(EDMD_outputs.kpm_modes, 1) == numel(evalues0)
    kpm_modes = EDMD_outputs.kpm_modes(idx_sorted_selected, :);
else
    kpm_modes = [];
end

prep = struct();
prep.dt = dt;
prep.dt_source = dt_source;
prep.time_axis = time_axis;
prep.selected_mode_idx_in_original = idx_sorted_selected(:);
prep.selected_mode_mask_in_original = selected_mode_mask;
prep.evalues_discrete = evalues_selected(:);
prep.evalues_bilinear = evalues_bilinear;
prep.efun_raw_time_by_mode = efun_raw;
prep.efun_feature_time_by_mode = efun_feature;
prep.kpm_modes_mode_by_dict = kpm_modes;
end


function ord = local_sort_eigenvalues(evalues, sort_by, sort_dir)
switch lower(sort_by)
    case {'modulus', 'abs'}
        key = abs(evalues);
    case {'real', 'realpart'}
        key = real(evalues);
    otherwise
        error('Unknown cfg.selection.sort_by = %s.', sort_by);
end

[~, ord] = sort(key, sort_dir);
end


function [dt, dt_source] = local_resolve_dt(input_cfg, EDMD_outputs)
dt = [];
dt_source = struct();
dt_source.source = 'missing';
dt_source.field = '';
dt_source.observable_file = '';
dt_source.message = '';

if isfield(input_cfg, 'dt') && ~isempty(input_cfg.dt)
    dt = local_positive_scalar(input_cfg.dt);
    if isempty(dt)
        error('cfg.input.dt must be empty or a positive finite scalar.');
    end
    dt_source.source = 'cfg.input.dt';
    dt_source.field = 'dt';
    return;
end

candidate_fields = {'dt', 'dx', 'sampling_period', 'sample_period'};
[dt, field_name] = local_read_first_positive_scalar(EDMD_outputs, candidate_fields);
if ~isempty(dt)
    dt_source.source = 'EDMD_outputs';
    dt_source.field = field_name;
    return;
end

observable_file = local_resolve_observable_file(input_cfg, EDMD_outputs);
dt_source.observable_file = observable_file;

if isempty(observable_file)
    dt_source.message = ['No dt/dx was found in EDMD_outputs, and no ', ...
        'cfg.input.observable_file or EDMD observable_file/data_full_path was available.'];
    return;
end

if exist(observable_file, 'file') ~= 2
    dt_source.message = sprintf(['No dt/dx was found in EDMD_outputs. ', ...
        'The configured observable file does not exist: %s'], observable_file);
    return;
end

[dt, observable_field, detail_msg] = local_read_dt_from_observable_file(observable_file);
if ~isempty(dt)
    dt_source.source = 'observable_file';
    dt_source.field = observable_field;
    dt_source.message = sprintf('Read sampling interval from %s.', observable_file);
else
    dt_source.message = sprintf(['No dt/dx was found in EDMD_outputs, and ', ...
        'observable file %s did not contain usable dx/dt/session_dx/fs metadata. %s'], ...
        observable_file, detail_msg);
end
end


function observable_file = local_resolve_observable_file(input_cfg, EDMD_outputs)
observable_file = '';

input_fields = {'observable_file', 'data_file', 'data_full_path'};
for i = 1:numel(input_fields)
    field_name = input_fields{i};
    if isfield(input_cfg, field_name)
        observable_file = local_text_scalar(input_cfg.(field_name));
        if ~isempty(observable_file)
            return;
        end
    end
end

edmd_fields = {'observable_file', 'data_full_path'};
for i = 1:numel(edmd_fields)
    field_name = edmd_fields{i};
    if isfield(EDMD_outputs, field_name)
        observable_file = local_text_scalar(EDMD_outputs.(field_name));
        if ~isempty(observable_file)
            return;
        end
    end
end
end


function [dt, field_name] = local_read_first_positive_scalar(S, field_names)
dt = [];
field_name = '';

for i = 1:numel(field_names)
    this_field = field_names{i};
    if ~isfield(S, this_field)
        continue;
    end

    dt = local_positive_scalar(S.(this_field));
    if ~isempty(dt)
        field_name = this_field;
        return;
    end
end
end


function [dt, field_name, detail_msg] = local_read_dt_from_observable_file(observable_file)
dt = [];
field_name = '';
detail_msg = '';

try
    vars = who('-file', observable_file);
catch ME
    detail_msg = sprintf('Could not inspect observable file: %s', ME.message);
    return;
end

fields_to_load = intersect( ...
    {'dx', 'dt', 'sampling_period', 'sample_period', ...
    'fs', 'sampling_frequency', 'session_dx', 'session_fs'}, ...
    vars, 'stable');

if isempty(fields_to_load)
    detail_msg = 'No candidate time-metadata variables were present.';
    return;
end

try
    S = load(observable_file, fields_to_load{:});
catch ME
    detail_msg = sprintf('Could not load observable time metadata: %s', ME.message);
    return;
end

[dt, field_name] = local_read_first_positive_scalar( ...
    S, {'dx', 'dt', 'sampling_period', 'sample_period'});
if ~isempty(dt)
    return;
end

if isfield(S, 'session_dx')
    dt = local_uniform_positive_scalar(S.session_dx);
    if ~isempty(dt)
        field_name = 'session_dx';
        return;
    end
    detail_msg = 'session_dx exists but is empty, non-positive, or inconsistent across sessions.';
end

[fs, fs_field] = local_read_first_positive_scalar(S, {'fs', 'sampling_frequency'});
if isempty(fs) && isfield(S, 'session_fs')
    fs = local_uniform_positive_scalar(S.session_fs);
    if ~isempty(fs)
        fs_field = 'session_fs';
    end
end

if ~isempty(fs)
    dt = 1 / fs;
    field_name = fs_field;
    return;
end

if isempty(detail_msg)
    detail_msg = 'Candidate variables were present but did not resolve to a positive scalar dt.';
end
end


function value = local_positive_scalar(value_in)
value = [];
if ~isnumeric(value_in) || isempty(value_in)
    return;
end

if ~isreal(value_in)
    value_in = real(value_in);
end

value_in = double(value_in);
if ~isscalar(value_in) || ~isfinite(value_in) || value_in <= 0
    return;
end

value = value_in;
end


function value = local_uniform_positive_scalar(value_in)
value = [];
if ~isnumeric(value_in) || isempty(value_in)
    return;
end

if ~isreal(value_in)
    value_in = real(value_in);
end

values = double(value_in(:));
values = values(isfinite(values) & values > 0);
if isempty(values)
    return;
end

value0 = median(values);
tol = max(1e-12, 1e-5 * abs(value0));
if all(abs(values - value0) <= tol)
    value = value0;
end
end


function txt = local_text_scalar(value)
txt = '';
if isstring(value)
    value = char(value);
end

if iscell(value) && isscalar(value)
    txt = local_text_scalar(value{1});
    return;
end

if ischar(value)
    txt = strtrim(value(:).');
end
end


function summary = local_build_summary(core, summary_cfg)
summary = struct();
summary.temporal_components_smooth_time_by_comp = [];

if ~summary_cfg.smooth.enable || isempty(core.temporal_components_time_by_comp)
    return;
end

switch lower(summary_cfg.smooth.method)
    case 'movmean'
        summary.temporal_components_smooth_time_by_comp = movmean( ...
            core.temporal_components_time_by_comp, ...
            summary_cfg.smooth.window, 1, ...
            'Endpoints', 'shrink');
    otherwise
        error('Unsupported summary smoothing method %s.', summary_cfg.smooth.method);
end
end


function [td_summary, artifacts] = local_run_thresholded_density(prep, result, td_cfg)
if isempty(prep.dt) && local_thresholded_density_needs_dt(td_cfg)
    error(['Thresholded density requires a sampling interval because its ', ...
        'window/step/smoothing parameters are specified in seconds. ', ...
        '%s'], prep.dt_source.message);
end

td_cfg.time_axis = prep.time_axis;
td_cfg.dt_source = prep.dt_source;
td_cfg.mode_index = result.input.mode_index;
td_cfg.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;

if ~isfield(td_cfg, 'title') || isempty(td_cfg.title)
    td_cfg.title = sprintf('%s Thresholded Eigenfunction Density', ...
        result.meta.method);
end

[D, fig] = get_thresholded_density( ...
    prep.efun_feature_time_by_mode, prep.dt, td_cfg);

if td_cfg.close_figure && ~isempty(fig) && isvalid(fig)
    close(fig);
end

td_summary = local_compact_thresholded_density(D);
artifacts = D.artifacts;
end


function tf = local_thresholded_density_needs_dt(td_cfg)
has_window_samples = isfield(td_cfg, 'window_samples') && ~isempty(td_cfg.window_samples);
has_step_samples = isfield(td_cfg, 'step_samples') && ~isempty(td_cfg.step_samples);

tf = ~(has_window_samples && has_step_samples);

smooth_enabled = isfield(td_cfg, 'smooth_density') && ~isempty(td_cfg.smooth_density) && ...
    logical(td_cfg.smooth_density);
has_smooth_bins = isfield(td_cfg, 'smooth_window_bins') && ~isempty(td_cfg.smooth_window_bins);
smooth_sec_positive = isfield(td_cfg, 'smooth_window_sec') && ...
    ~isempty(td_cfg.smooth_window_sec) && td_cfg.smooth_window_sec > 0;

tf = tf || (smooth_enabled && smooth_sec_positive && ~has_smooth_bins);
end


function td_summary = local_compact_thresholded_density(D)
td_summary = struct();
td_summary.created_at = D.created_at;
td_summary.input = D.input;
td_summary.params = D.params;
td_summary.summary = D.summary;
td_summary.artifacts = D.artifacts;
td_summary.threshold_by_mode = D.threshold_by_mode;
td_summary.density_size = size(D.density_time_by_mode);
td_summary.t_range = [D.t_centers(1), D.t_centers(end)];
td_summary.n_windows = numel(D.t_centers);
td_summary.n_modes = size(D.density_time_by_mode, 2);
end


function [te_summary, artifacts] = local_run_thresholded_events(prep, result, te_cfg)
if isempty(prep.dt)
    error('Thresholded events require a sampling interval. %s', ...
        prep.dt_source.message);
end

te_cfg.time_axis = prep.time_axis;
te_cfg.dt_source = prep.dt_source;
te_cfg.mode_index = result.input.mode_index;
te_cfg.selected_mode_idx_in_original = prep.selected_mode_idx_in_original;

if ~isfield(te_cfg, 'title') || isempty(te_cfg.title)
    te_cfg.title = sprintf('%s Thresholded Eigenfunction Events', ...
        result.meta.method);
end

[E, figs] = get_thresholded_events( ...
    prep.efun_feature_time_by_mode, prep.dt, prep.evalues_discrete, te_cfg);

if te_cfg.close_figure && isfield(figs, 'summary') && ...
        ~isempty(figs.summary) && isvalid(figs.summary)
    close(figs.summary);
end

te_summary = local_compact_thresholded_events(E);
artifacts = E.artifacts;
end


function te_summary = local_compact_thresholded_events(E)
te_summary = struct();
te_summary.created_at = E.created_at;
te_summary.input = E.input;
te_summary.params = E.params;
te_summary.summary = E.summary;
te_summary.artifacts = E.artifacts;
te_summary.threshold_by_mode = E.threshold_by_mode;
te_summary.mode_timescales = E.mode_timescales;
te_summary.event_rate_size = size(E.event_rate_time_by_mode);
te_summary.t_range = [E.t_centers(1), E.t_centers(end)];
te_summary.n_windows = numel(E.t_centers);
te_summary.n_modes = size(E.event_rate_time_by_mode, 2);
te_summary.n_events = height(E.event_table);
end


function [td_summary, artifacts] = local_run_dimred_thresholded_density( ...
        prep, result, td_cfg)
if isempty(result.core.temporal_components_time_by_comp)
    error('Dimension-reduced thresholded density requires temporal components.');
end

if isempty(prep.dt) && local_thresholded_density_needs_dt(td_cfg)
    error(['Dimension-reduced thresholded density requires a sampling ', ...
        'interval because its window/step/smoothing parameters are ', ...
        'specified in seconds. %s'], prep.dt_source.message);
end

td_cfg.time_axis = prep.time_axis;
td_cfg.dt_source = prep.dt_source;
td_cfg.component_index = (1:size(result.core.temporal_components_time_by_comp, 2)).';

if ~isfield(td_cfg, 'title') || isempty(td_cfg.title)
    td_cfg.title = sprintf('%s Thresholded Component Density', ...
        result.meta.method);
end

[D, ~] = get_dimred_thresholded_density(result, prep.dt, td_cfg);

td_summary = local_compact_dimred_thresholded_density(D);
artifacts = D.artifacts;
end


function td_summary = local_compact_dimred_thresholded_density(D)
td_summary = struct();
td_summary.created_at = D.created_at;
td_summary.meta = D.meta;
td_summary.input = D.input;
td_summary.params = D.params;
td_summary.summary = D.summary;
td_summary.artifacts = D.artifacts;
td_summary.threshold_by_component = D.threshold_by_component;
td_summary.component_index = D.component_index;
td_summary.density_size = size(D.density_time_by_component);
td_summary.t_range = [D.t_centers(1), D.t_centers(end)];
td_summary.n_windows = numel(D.t_centers);
td_summary.n_components = size(D.density_time_by_component, 2);
end


function [te_summary, artifacts] = local_run_dimred_thresholded_events( ...
        prep, result, te_cfg)
if isempty(result.core.temporal_components_time_by_comp)
    error('Dimension-reduced thresholded events require temporal components.');
end

if isempty(prep.dt)
    error(['Dimension-reduced thresholded events require a sampling ', ...
        'interval. %s'], prep.dt_source.message);
end

te_cfg.time_axis = prep.time_axis;
te_cfg.dt_source = prep.dt_source;
te_cfg.component_index = (1:size(result.core.temporal_components_time_by_comp, 2)).';

if ~isfield(te_cfg, 'title') || isempty(te_cfg.title)
    te_cfg.title = sprintf('%s Thresholded Component Events', ...
        result.meta.method);
end

[E, ~] = get_dimred_thresholded_events(result, prep.dt, te_cfg);

te_summary = local_compact_dimred_thresholded_events(E);
artifacts = E.artifacts;
end


function te_summary = local_compact_dimred_thresholded_events(E)
te_summary = struct();
te_summary.created_at = E.created_at;
te_summary.meta = E.meta;
te_summary.input = E.input;
te_summary.params = E.params;
te_summary.summary = E.summary;
te_summary.artifacts = E.artifacts;
te_summary.threshold_by_component = E.threshold_by_component;
te_summary.component_index = E.component_index;
te_summary.component_evalues_discrete = E.component_evalues_discrete;
te_summary.component_evalue_info = E.component_evalue_info;
te_summary.component_timescales = E.component_timescales;
te_summary.event_rate_size = size(E.event_rate_time_by_component);
te_summary.t_range = [E.t_centers(1), E.t_centers(end)];
te_summary.n_windows = numel(E.t_centers);
te_summary.n_components = size(E.event_rate_time_by_component, 2);
te_summary.n_events = height(E.event_table);
end


function result_to_save = local_build_save_result(result, payload)
if isstring(payload)
    payload = char(payload);
end
payload = lower(strtrim(payload));

switch payload
    case {'full', 'all'}
        result_to_save = result;

    case {'compact', 'minimal'}
        result_to_save = result;
        result_to_save.meta.save_payload = 'compact';

        if isfield(result_to_save, 'data') && isstruct(result_to_save.data)
            result_to_save.data = local_compact_data(result_to_save.data);
        end

        if isfield(result_to_save, 'core') && isstruct(result_to_save.core)
            result_to_save.core = local_compact_core(result_to_save.core);
        end

        if isfield(result_to_save, 'aux') && isstruct(result_to_save.aux)
            result_to_save.aux = local_compact_aux(result_to_save.aux);
        end

        if isfield(result_to_save, 'concat') && isstruct(result_to_save.concat)
            result_to_save.concat = local_compact_concat(result_to_save.concat);
        end

        if isfield(result_to_save, 'cfg') && isstruct(result_to_save.cfg)
            result_to_save.cfg = local_compact_cfg(result_to_save.cfg);
        end

    otherwise
        error('Unknown cfg.save.payload = %s. Use ''compact'' or ''full''.', payload);
end
end


function data = local_compact_data(data)
remove_fields = {'efun_raw_time_by_mode', 'efun_feature_time_by_mode', ...
    'kpm_modes_mode_by_dict'};
omitted = local_existing_fields(data, remove_fields);
data = local_rmfield_if_present(data, remove_fields);
if ~isempty(omitted)
    data.omitted_compact_fields = omitted;
end
end


function core = local_compact_core(core)
remove_fields = {'reconstruction_time_by_mode'};
omitted = local_existing_fields(core, remove_fields);
core = local_rmfield_if_present(core, remove_fields);
if ~isempty(omitted)
    core.omitted_compact_fields = omitted;
end
end


function aux = local_compact_aux(aux)
if isfield(aux, 'spectrum') && isstruct(aux.spectrum)
    remove_fields = {'distance_mode_by_mode'};
    omitted = local_existing_fields(aux.spectrum, remove_fields);
    aux.spectrum = local_rmfield_if_present(aux.spectrum, remove_fields);

    if isfield(aux.spectrum, 'embedding_info') && ...
            isstruct(aux.spectrum.embedding_info)
        aux.spectrum.embedding_info = local_rmfield_if_present( ...
            aux.spectrum.embedding_info, {'gm'});
    end

    if ~isempty(omitted)
        aux.spectrum.omitted_compact_fields = omitted;
    end
end
end


function concat = local_compact_concat(concat)
if isfield(concat, 'chunk_ids') && ~isempty(concat.chunk_ids)
    concat.chunk_id_range = [concat.chunk_ids(1), concat.chunk_ids(end)];
end

remove_fields = {'files', 'chunk_ids', 'chunk_lengths', ...
    'chunk_start_idx', 'chunk_end_idx'};
omitted = local_existing_fields(concat, remove_fields);
concat = local_rmfield_if_present(concat, remove_fields);
if ~isempty(omitted)
    concat.omitted_compact_fields = omitted;
end
end


function cfg = local_compact_cfg(cfg)
cfg = local_strip_preloaded_source_payload(cfg);

if isfield(cfg, 'viz') && isstruct(cfg.viz)
    if isfield(cfg.viz, 'spectrum') && isstruct(cfg.viz.spectrum)
        cfg.viz.spectrum = local_rmfield_if_present(cfg.viz.spectrum, ...
            {'distance_colormap', 'embedding_index_colormap', ...
            'cluster_colormap'});
    end

    if isfield(cfg.viz, 'state_space') && isstruct(cfg.viz.state_space)
        cfg.viz.state_space = local_rmfield_if_present(cfg.viz.state_space, ...
            {'time_colormap', 'value_colormap'});
    end
end
end


function cfg = local_strip_preloaded_source_payload(cfg)
if isfield(cfg, 'source') && isstruct(cfg.source)
    cfg.source = local_rmfield_if_present(cfg.source, ...
        {'preloaded_EDMD_outputs', 'preloaded_concat_info', ...
        'preloaded_source_info'});
end
end


function save_dir = local_default_save_dir(cfg)
source_run_name = local_source_run_name(cfg);

if isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    save_dir = fullfile(get_project_processed_root(), cfg.dataset.name, ...
        'efun', ...
        source_run_name, 'mat');
else
    save_dir = fullfile(get_project_processed_root(), ...
        'efun', ...
        source_run_name, 'mat');
end
end


function file_stem = local_default_file_stem(cfg)
if isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    file_stem = sprintf('%s_efun', cfg.dataset.name);
else
    file_stem = 'efun';
end
end


function source_run_name = local_source_run_name(cfg)
source_run_name = 'unspecified_source';

if isfield(cfg, 'output') && isfield(cfg.output, 'output_run_name') && ...
        ~isempty(cfg.output.output_run_name)
    source_run_name = cfg.output.output_run_name;
elseif isfield(cfg, 'output') && isfield(cfg.output, 'source_run_name') && ...
        ~isempty(cfg.output.source_run_name)
    source_run_name = cfg.output.source_run_name;
elseif isfield(cfg, 'source') && isfield(cfg.source, 'data_dir') && ...
        ~isempty(cfg.source.data_dir)
    [~, source_run_name] = fileparts(cfg.source.data_dir);
elseif isfield(cfg, 'source') && isfield(cfg.source, 'edmd_file') && ...
        ~isempty(cfg.source.edmd_file)
    [~, source_run_name] = fileparts(cfg.source.edmd_file);
end

source_run_name = regexprep(char(source_run_name), '[^\w\-]+', '_');
end


function fields = local_existing_fields(S, names)
fields = {};
for i = 1:numel(names)
    if isfield(S, names{i})
        fields{end+1, 1} = names{i}; %#ok<AGROW>
    end
end
end


function S = local_rmfield_if_present(S, names)
for i = 1:numel(names)
    if isfield(S, names{i})
        S = rmfield(S, names{i});
    end
end
end


function cfg = local_apply_thresholded_density_defaults(cfg)
if ~isfield(cfg, 'thresholded_density') || ~isstruct(cfg.thresholded_density)
    cfg.thresholded_density = struct();
end

td = cfg.thresholded_density;

if ~isfield(td, 'enable') || isempty(td.enable)
    td.enable = false;
end

if ~isfield(td, 'window_sec') || isempty(td.window_sec)
    td.window_sec = 2;
end

if ~isfield(td, 'step_sec') || isempty(td.step_sec)
    td.step_sec = td.window_sec;
end

if ~isfield(td, 'threshold_ratio') || isempty(td.threshold_ratio)
    td.threshold_ratio = 0.7;
end

if ~isfield(td, 'threshold_mode') || isempty(td.threshold_mode)
    td.threshold_mode = 'meanplusstd';
end

if ~isfield(td, 'value_transform') || isempty(td.value_transform)
    td.value_transform = 'none';
end

if ~isfield(td, 'density_denominator') || isempty(td.density_denominator)
    td.density_denominator = 'window_samples';
end

if ~isfield(td, 'output_class') || isempty(td.output_class)
    td.output_class = 'single';
end

if ~isfield(td, 'smooth_density') || isempty(td.smooth_density)
    td.smooth_density = false;
end

if ~isfield(td, 'smooth_window_sec') || isempty(td.smooth_window_sec)
    td.smooth_window_sec = 0;
end

if ~isfield(td, 'require_session_metadata') || isempty(td.require_session_metadata)
    td.require_session_metadata = false;
end

if ~isfield(td, 'observable_file') || isempty(td.observable_file)
    if isfield(cfg, 'input') && isfield(cfg.input, 'observable_file')
        td.observable_file = cfg.input.observable_file;
    else
        td.observable_file = '';
    end
end

if ~isfield(td, 'save_results') || isempty(td.save_results)
    td.save_results = td.enable;
end

if ~isfield(td, 'make_figure') || isempty(td.make_figure)
    td.make_figure = td.enable;
end

if ~isfield(td, 'save_figure') || isempty(td.save_figure)
    td.save_figure = td.enable;
end

if ~isfield(td, 'close_figure') || isempty(td.close_figure)
    td.close_figure = true;
end

if ~isfield(td, 'save_v7_3') || isempty(td.save_v7_3)
    td.save_v7_3 = true;
end

output_root = local_output_root(cfg);
if ~isfield(td, 'save_dir') || isempty(td.save_dir)
    td.save_dir = fullfile(output_root, 'thresholded_density', 'mat');
end

if ~isfield(td, 'figure_dir') || isempty(td.figure_dir)
    td.figure_dir = fullfile(output_root, 'thresholded_density', 'fig');
end

if ~isfield(td, 'save_stem') || isempty(td.save_stem)
    td.save_stem = sprintf('%s_thresholded_density', local_dataset_name(cfg));
end

if ~isfield(td, 'save_tag')
    td.save_tag = '';
end

if ~isfield(td, 'figure_visible') || isempty(td.figure_visible)
    td.figure_visible = 'off';
end

if ~isfield(td, 'figure_resolution') || isempty(td.figure_resolution)
    td.figure_resolution = 200;
end

if ~isfield(td, 'max_plot_modes') || isempty(td.max_plot_modes)
    td.max_plot_modes = 100;
end

if ~isfield(td, 'title') || isempty(td.title)
    td.title = 'Thresholded Eigenfunction Density';
end

if ~isfield(td, 'colormap') || isempty(td.colormap)
    td.colormap = 'turbo';
end

if ~isfield(td, 'background_color') || isempty(td.background_color)
    td.background_color = [0, 0, 0];
end

if ~isfield(td, 'axes_color') || isempty(td.axes_color)
    td.axes_color = td.background_color;
end

if ~isfield(td, 'text_color') || isempty(td.text_color)
    td.text_color = [1, 1, 1];
end

if ~isfield(td, 'grid_color') || isempty(td.grid_color)
    td.grid_color = [0.75, 0.75, 0.75];
end

cfg.thresholded_density = td;
end


function cfg = local_apply_thresholded_events_defaults(cfg)
if ~isfield(cfg, 'thresholded_events') || ~isstruct(cfg.thresholded_events)
    cfg.thresholded_events = struct();
end

te = cfg.thresholded_events;
td = cfg.thresholded_density;

if ~isfield(te, 'enable') || isempty(te.enable)
    te.enable = false;
end

if ~isfield(te, 'threshold_ratio') || isempty(te.threshold_ratio)
    te.threshold_ratio = td.threshold_ratio;
end

if ~isfield(te, 'threshold_mode') || isempty(te.threshold_mode)
    te.threshold_mode = td.threshold_mode;
end

if ~isfield(te, 'value_transform') || isempty(te.value_transform)
    te.value_transform = td.value_transform;
end

if ~isfield(te, 'event_detector') || isempty(te.event_detector)
    te.event_detector = 'find_peak_loc';
end

if ~isfield(te, 'event_score') || isempty(te.event_score)
    te.event_score = 'area_above_threshold';
end

if ~isfield(te, 'min_duration_sec') || isempty(te.min_duration_sec)
    te.min_duration_sec = 0.03;
end

if ~isfield(te, 'merge_gap_sec') || isempty(te.merge_gap_sec)
    te.merge_gap_sec = 0.02;
end

if ~isfield(te, 'find_peak_loc_window_sec')
    te.find_peak_loc_window_sec = [];
end

if ~isfield(te, 'find_peak_loc_window_samples')
    te.find_peak_loc_window_samples = [];
end

if ~isfield(te, 'find_peak_loc_drop_first') || isempty(te.find_peak_loc_drop_first)
    te.find_peak_loc_drop_first = false;
end

if ~isfield(te, 'find_peak_loc_post_nms') || isempty(te.find_peak_loc_post_nms)
    te.find_peak_loc_post_nms = false;
end

if ~isfield(te, 'bin_sec') || isempty(te.bin_sec)
    te.bin_sec = td.window_sec;
end

if ~isfield(te, 'step_sec') || isempty(te.step_sec)
    te.step_sec = te.bin_sec;
end

if ~isfield(te, 'smooth_rate') || isempty(te.smooth_rate)
    te.smooth_rate = true;
end

if ~isfield(te, 'smooth_sigma_sec') || isempty(te.smooth_sigma_sec)
    te.smooth_sigma_sec = te.bin_sec;
end

if ~isfield(te, 'output_class') || isempty(te.output_class)
    te.output_class = 'single';
end

if ~isfield(te, 'require_session_metadata') || isempty(te.require_session_metadata)
    te.require_session_metadata = true;
end

if ~isfield(te, 'observable_file') || isempty(te.observable_file)
    if isfield(cfg, 'input') && isfield(cfg.input, 'observable_file')
        te.observable_file = cfg.input.observable_file;
    else
        te.observable_file = '';
    end
end

if ~isfield(te, 'nms') || ~isstruct(te.nms)
    te.nms = struct();
end
if ~isfield(te.nms, 'enable') || isempty(te.nms.enable)
    te.nms.enable = true;
end
if ~isfield(te.nms, 'mode') || isempty(te.nms.mode)
    te.nms.mode = 'eigen_period';
end
if ~isfield(te.nms, 'period_fraction') || isempty(te.nms.period_fraction)
    te.nms.period_fraction = 1.0;
end
if ~isfield(te.nms, 'fixed_sec') || isempty(te.nms.fixed_sec)
    te.nms.fixed_sec = 0.25;
end
if ~isfield(te.nms, 'default_sec') || isempty(te.nms.default_sec)
    te.nms.default_sec = 0.25;
end
if ~isfield(te.nms, 'min_sec') || isempty(te.nms.min_sec)
    te.nms.min_sec = 0.05;
end
if ~isfield(te.nms, 'max_sec') || isempty(te.nms.max_sec)
    te.nms.max_sec = 2.0;
end
if ~isfield(te.nms, 'angle_eps') || isempty(te.nms.angle_eps)
    te.nms.angle_eps = 1e-4;
end

if ~isfield(te, 'save_results') || isempty(te.save_results)
    te.save_results = te.enable;
end

if ~isfield(te, 'make_figure') || isempty(te.make_figure)
    te.make_figure = te.enable;
end

if ~isfield(te, 'save_figure') || isempty(te.save_figure)
    te.save_figure = te.enable;
end

if ~isfield(te, 'close_figure') || isempty(te.close_figure)
    te.close_figure = true;
end

if ~isfield(te, 'save_v7_3') || isempty(te.save_v7_3)
    te.save_v7_3 = true;
end

output_root = local_output_root(cfg);
if ~isfield(te, 'save_dir') || isempty(te.save_dir)
    te.save_dir = fullfile(output_root, 'thresholded_events', 'mat');
end

if ~isfield(te, 'figure_dir') || isempty(te.figure_dir)
    te.figure_dir = fullfile(output_root, 'thresholded_events', 'fig');
end

if ~isfield(te, 'save_stem') || isempty(te.save_stem)
    te.save_stem = sprintf('%s_thresholded_events', local_dataset_name(cfg));
end

if ~isfield(te, 'save_tag')
    te.save_tag = '';
end

if ~isfield(te, 'figure_visible') || isempty(te.figure_visible)
    te.figure_visible = 'off';
end

if ~isfield(te, 'figure_resolution') || isempty(te.figure_resolution)
    te.figure_resolution = 200;
end

if ~isfield(te, 'max_plot_modes') || isempty(te.max_plot_modes)
    te.max_plot_modes = 80;
end

if ~isfield(te, 'max_plot_events') || isempty(te.max_plot_events)
    te.max_plot_events = 5000;
end

if ~isfield(te, 'title') || isempty(te.title)
    te.title = 'Thresholded Eigenfunction Events';
end

if ~isfield(te, 'colormap') || isempty(te.colormap)
    te.colormap = 'turbo';
end

if ~isfield(te, 'background_color') || isempty(te.background_color)
    te.background_color = [0, 0, 0];
end

if ~isfield(te, 'axes_color') || isempty(te.axes_color)
    te.axes_color = te.background_color;
end

if ~isfield(te, 'text_color') || isempty(te.text_color)
    te.text_color = [1, 1, 1];
end

if ~isfield(te, 'grid_color') || isempty(te.grid_color)
    te.grid_color = [0.75, 0.75, 0.75];
end

cfg.thresholded_events = te;
end


function cfg = local_apply_dimred_thresholded_density_defaults(cfg)
if ~isfield(cfg, 'dimred_thresholded_density') || ...
        ~isstruct(cfg.dimred_thresholded_density)
    cfg.dimred_thresholded_density = struct();
end

dtd = cfg.dimred_thresholded_density;
td = cfg.thresholded_density;

if ~isfield(dtd, 'enable') || isempty(dtd.enable)
    dtd.enable = false;
end

dtd = local_copy_missing_fields(dtd, td, { ...
    'window_sec', 'step_sec', 'window_samples', 'step_samples', ...
    'threshold_ratio', 'threshold_mode', 'value_transform', ...
    'density_denominator', 'output_class', 'smooth_density', ...
    'smooth_window_sec', 'smooth_window_bins', 'save_v7_3', ...
    'require_session_metadata', 'observable_file', ...
    'session_start_idx', 'session_end_idx', 'session_lengths', ...
    'session_ids', 'session_dx', ...
    'figure_visible', 'figure_resolution', 'colormap', ...
    'background_color', 'axes_color', 'text_color', 'grid_color'});

if ~isfield(dtd, 'save_results') || isempty(dtd.save_results)
    dtd.save_results = dtd.enable;
end
if ~isfield(dtd, 'make_figure') || isempty(dtd.make_figure)
    dtd.make_figure = dtd.enable;
end
if ~isfield(dtd, 'save_figure') || isempty(dtd.save_figure)
    dtd.save_figure = dtd.enable;
end
if ~isfield(dtd, 'close_figure') || isempty(dtd.close_figure)
    dtd.close_figure = true;
end

output_root = local_output_root(cfg);
if ~isfield(dtd, 'save_dir') || isempty(dtd.save_dir)
    dtd.save_dir = fullfile(output_root, 'dimred_thresholded_density', 'mat');
end
if ~isfield(dtd, 'figure_dir') || isempty(dtd.figure_dir)
    dtd.figure_dir = fullfile(output_root, 'dimred_thresholded_density', 'fig');
end
if ~isfield(dtd, 'save_stem') || isempty(dtd.save_stem)
    dtd.save_stem = sprintf('%s_dimred_thresholded_density', ...
        local_dataset_name(cfg));
end
if ~isfield(dtd, 'save_tag')
    dtd.save_tag = '';
end
if ~isfield(dtd, 'max_plot_modes') || isempty(dtd.max_plot_modes)
    dtd.max_plot_modes = 40;
end
if ~isfield(dtd, 'title') || isempty(dtd.title)
    dtd.title = 'Thresholded Eigenfunction Component Density';
end

cfg.dimred_thresholded_density = dtd;
end


function cfg = local_apply_dimred_thresholded_events_defaults(cfg)
if ~isfield(cfg, 'dimred_thresholded_events') || ...
        ~isstruct(cfg.dimred_thresholded_events)
    cfg.dimred_thresholded_events = struct();
end

dte = cfg.dimred_thresholded_events;
te = cfg.thresholded_events;

if ~isfield(dte, 'enable') || isempty(dte.enable)
    dte.enable = false;
end

dte = local_copy_missing_fields(dte, te, { ...
    'threshold_ratio', 'threshold_mode', 'value_transform', ...
    'event_detector', 'event_score', 'min_duration_sec', 'merge_gap_sec', ...
    'min_duration_samples', 'merge_gap_samples', 'bin_sec', ...
    'step_sec', 'bin_samples', 'step_samples', 'smooth_rate', ...
    'smooth_sigma_sec', 'output_class', 'require_session_metadata', ...
    'find_peak_loc_window_sec', 'find_peak_loc_window_samples', ...
    'find_peak_loc_drop_first', 'find_peak_loc_post_nms', ...
    'observable_file', 'save_v7_3', 'figure_visible', ...
    'figure_resolution', 'max_plot_events', 'colormap', ...
    'background_color', 'axes_color', 'text_color', 'grid_color'});

if ~isfield(dte, 'component_evalue_weight_transform') || ...
        isempty(dte.component_evalue_weight_transform)
    dte.component_evalue_weight_transform = 'abs';
end
if ~isfield(dte, 'component_evalue_angle_statistic') || ...
        isempty(dte.component_evalue_angle_statistic)
    dte.component_evalue_angle_statistic = 'weighted_median_abs';
end
if ~isfield(dte, 'component_evalue_top_modes') || ...
        isempty(dte.component_evalue_top_modes)
    dte.component_evalue_top_modes = 5;
end

if ~isfield(dte, 'nms') || ~isstruct(dte.nms)
    dte.nms = struct();
end
if isfield(te, 'nms') && isstruct(te.nms)
    dte.nms = local_copy_missing_fields(dte.nms, te.nms, { ...
        'enable', 'mode', 'period_fraction', 'fixed_sec', ...
        'default_sec', 'min_sec', 'max_sec', 'angle_eps'});
end
if ~isfield(dte.nms, 'enable') || isempty(dte.nms.enable)
    dte.nms.enable = true;
end
if ~isfield(dte.nms, 'mode') || isempty(dte.nms.mode)
    dte.nms.mode = 'eigen_period';
end
if ~isfield(dte.nms, 'period_fraction') || isempty(dte.nms.period_fraction)
    dte.nms.period_fraction = 1.0;
end
if ~isfield(dte.nms, 'fixed_sec') || isempty(dte.nms.fixed_sec)
    dte.nms.fixed_sec = 0.25;
end
if ~isfield(dte.nms, 'default_sec') || isempty(dte.nms.default_sec)
    dte.nms.default_sec = 0.25;
end
if ~isfield(dte.nms, 'min_sec') || isempty(dte.nms.min_sec)
    dte.nms.min_sec = 0.05;
end
if ~isfield(dte.nms, 'max_sec') || isempty(dte.nms.max_sec)
    dte.nms.max_sec = 2.0;
end
if ~isfield(dte.nms, 'angle_eps') || isempty(dte.nms.angle_eps)
    dte.nms.angle_eps = 1e-4;
end

if ~isfield(dte, 'save_results') || isempty(dte.save_results)
    dte.save_results = dte.enable;
end
if ~isfield(dte, 'make_figure') || isempty(dte.make_figure)
    dte.make_figure = dte.enable;
end
if ~isfield(dte, 'save_figure') || isempty(dte.save_figure)
    dte.save_figure = dte.enable;
end
if ~isfield(dte, 'close_figure') || isempty(dte.close_figure)
    dte.close_figure = true;
end

output_root = local_output_root(cfg);
if ~isfield(dte, 'save_dir') || isempty(dte.save_dir)
    dte.save_dir = fullfile(output_root, 'dimred_thresholded_events', 'mat');
end
if ~isfield(dte, 'figure_dir') || isempty(dte.figure_dir)
    dte.figure_dir = fullfile(output_root, 'dimred_thresholded_events', 'fig');
end
if ~isfield(dte, 'save_stem') || isempty(dte.save_stem)
    dte.save_stem = sprintf('%s_dimred_thresholded_events', ...
        local_dataset_name(cfg));
end
if ~isfield(dte, 'save_tag')
    dte.save_tag = '';
end
if ~isfield(dte, 'max_plot_modes') || isempty(dte.max_plot_modes)
    dte.max_plot_modes = 40;
end
if ~isfield(dte, 'title') || isempty(dte.title)
    dte.title = 'Thresholded Eigenfunction Component Events';
end

cfg.dimred_thresholded_events = dte;
end


function dst = local_copy_missing_fields(dst, src, names)
for i = 1:numel(names)
    name = names{i};
    if (~isfield(dst, name) || isempty(dst.(name))) && isfield(src, name)
        dst.(name) = src.(name);
    end
end
end


function output_root = local_output_root(cfg)
if isfield(cfg, 'output') && isfield(cfg.output, 'root') && ...
        ~isempty(cfg.output.root)
    output_root = cfg.output.root;
elseif isfield(cfg, 'save') && isfield(cfg.save, 'dir') && ...
        ~isempty(cfg.save.dir)
    output_root = fileparts(cfg.save.dir);
else
    output_root = pwd;
end
end


function dataset_name = local_dataset_name(cfg)
dataset_name = 'efun';
if isfield(cfg, 'dataset') && isfield(cfg.dataset, 'name') && ...
        ~isempty(cfg.dataset.name)
    dataset_name = char(cfg.dataset.name);
end
dataset_name = regexprep(dataset_name, '[^\w\-]+', '_');
end


function save_path = local_build_save_path(cfg, result)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {cfg.save.file_stem, cfg.path.kind, lower(cfg.feature.variant), lower(result.meta.method)};

if isfield(cfg.save, 'tag') && ~isempty(cfg.save.tag)
    pieces{end+1} = cfg.save.tag;
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(cfg.save.dir, sprintf('%s__%s.mat', filename, timestamp));
end
