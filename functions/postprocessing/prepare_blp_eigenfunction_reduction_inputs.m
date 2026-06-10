function prep = prepare_blp_eigenfunction_reduction_inputs(EDMD_outputs, cfg)
%PREPARE_BLP_EIGENFUNCTION_REDUCTION_INPUTS Build aligned reduction inputs from EDMD outputs.

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
