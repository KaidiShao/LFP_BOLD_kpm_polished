function [prep, B] = prepare_bold_eigenfunction_reduction_inputs(bold_post_input, cfg)
%PREPARE_BOLD_EIGENFUNCTION_REDUCTION_INPUTS Build P9 inputs from BOLD_POST.

if nargin < 1 || isempty(bold_post_input)
    error('bold_post_input is required.');
end
if nargin < 2 || isempty(cfg)
    cfg = struct();
end
cfg = local_apply_defaults(cfg);

B = local_load_bold_post(bold_post_input);
if ~isfield(B, 'EDMD_outputs') || ~isstruct(B.EDMD_outputs)
    error('BOLD post input does not contain EDMD_outputs.');
end

E = B.EDMD_outputs;
[X0, evalues0, kpm0, original_idx0, feature_info] = ...
    local_select_feature_matrix(E, cfg.feature.name);

if size(X0, 2) ~= numel(evalues0)
    error('Feature matrix columns (%d) must match evalues (%d).', ...
        size(X0, 2), numel(evalues0));
end

ord = local_sort_modes(evalues0, cfg.selection.sort_by, cfg.selection.sort_dir);
evalues_sorted = evalues0(ord);
mask_sorted = abs(evalues_sorted) > cfg.selection.abs_thresh;
idx_selected = ord(mask_sorted);
if isempty(idx_selected)
    error('No BOLD modes remain after applying abs threshold %g.', ...
        cfg.selection.abs_thresh);
end

if isfinite(cfg.selection.max_modes)
    idx_selected = idx_selected(1:min(numel(idx_selected), cfg.selection.max_modes));
end

X_raw = X0(:, idx_selected);
evalues_selected = evalues0(idx_selected);
original_idx = original_idx0(idx_selected);

switch lower(cfg.feature.normalization)
    case 'maxabs_per_mode'
        X_feature = normalize_efun(X_raw, feature_info.variant);
    case {'none', 'raw'}
        switch feature_info.variant
            case 'abs'
                X_feature = abs(X_raw);
            case 'real'
                X_feature = real(X_raw);
            otherwise
                X_feature = X_raw;
        end
    otherwise
        error('Unsupported cfg.feature.normalization = %s.', cfg.feature.normalization);
end

[dt, dt_source] = local_resolve_dt(B, E);
if isempty(dt)
    time_axis = (1:size(X_raw, 1)).';
    evalues_bilinear = [];
else
    time_axis = (0:size(X_raw, 1)-1).' * dt;
    evalues_bilinear = (2 / dt) * (evalues_selected - 1) ./ (evalues_selected + 1);
end

if isempty(kpm0)
    kpm = [];
else
    kpm = kpm0(idx_selected, :);
end

finite_original_idx = original_idx(isfinite(original_idx) & original_idx >= 1);
if isempty(finite_original_idx)
    mask_len = numel(evalues0);
else
    mask_len = max(numel(evalues0), max(round(finite_original_idx)));
end
selected_mode_mask = false(mask_len, 1);
valid_original_idx = original_idx(isfinite(original_idx) & original_idx >= 1 & ...
    original_idx <= mask_len);
selected_mode_mask(round(valid_original_idx)) = true;

prep = struct();
prep.dt = dt;
prep.dt_source = dt_source;
prep.time_axis = time_axis;
prep.selected_mode_idx_in_original = original_idx(:);
prep.selected_mode_idx_in_feature = idx_selected(:);
prep.selected_mode_mask_in_original = selected_mode_mask;
prep.evalues_discrete = evalues_selected(:);
prep.evalues_bilinear = evalues_bilinear;
prep.bold_feature_raw_time_by_mode = X_raw;
prep.bold_feature_time_by_mode = X_feature;
prep.efun_raw_time_by_mode = X_raw;
prep.efun_feature_time_by_mode = X_feature;
prep.kpm_modes_mode_by_dict = kpm;
prep.feature = feature_info;
prep.run_info = local_get_field(B, 'run_info', struct());
prep.observable_file = local_get_field(B, 'observable_file', '');
end


function cfg = local_apply_defaults(cfg)
if ~isfield(cfg, 'feature') || ~isstruct(cfg.feature)
    cfg.feature = struct();
end
if ~isfield(cfg.feature, 'name') || isempty(cfg.feature.name)
    cfg.feature.name = 'efun_abs';
end
if ~isfield(cfg.feature, 'normalization') || isempty(cfg.feature.normalization)
    cfg.feature.normalization = 'maxabs_per_mode';
end
if ~isfield(cfg, 'selection') || ~isstruct(cfg.selection)
    cfg.selection = struct();
end
cfg.selection = local_set_default(cfg.selection, 'abs_thresh', 0.01);
cfg.selection = local_set_default(cfg.selection, 'sort_by', 'preserve');
cfg.selection = local_set_default(cfg.selection, 'sort_dir', 'descend');
cfg.selection = local_set_default(cfg.selection, 'max_modes', Inf);
end


function B = local_load_bold_post(input)
if ischar(input) || isstring(input)
    file = char(string(input));
    S = load_mat_file_with_short_path(file, 'BOLD_POST');
    if ~isfield(S, 'BOLD_POST')
        error('BOLD_POST variable missing in %s.', file);
    end
    B = S.BOLD_POST;
    B.source_file = file;
elseif isstruct(input)
    B = input;
    if ~isfield(B, 'source_file')
        B.source_file = '';
    end
else
    error('bold_post_input must be a path or struct.');
end
end


function [X, evalues, kpm, original_idx, info] = local_select_feature_matrix(E, feature_name)
feature_name = lower(char(string(feature_name)));
if contains(feature_name, 'deconv')
    source_kind = 'deconv_efun';
else
    source_kind = 'efun';
end
if contains(feature_name, 'real')
    variant = 'real';
else
    variant = 'abs';
end

switch source_kind
    case 'efun'
        if ~isfield(E, 'efuns') || isempty(E.efuns) || ...
                ~isfield(E, 'evalues') || isempty(E.evalues)
            error('EDMD_outputs.efuns and EDMD_outputs.evalues are required.');
        end
        X = E.efuns;
        evalues = E.evalues(:);
        raw_field = 'EDMD_outputs.efuns';
        kpm = local_kpm_from_selected(E, size(X, 2));
        original_idx = local_original_indices_from_selected(E, size(X, 2));

    case 'deconv_efun'
        if ~isfield(E, 'deconv_efuns') || ~isstruct(E.deconv_efuns)
            error('EDMD_outputs.deconv_efuns is required for feature %s.', feature_name);
        end
        D = E.deconv_efuns;
        if isfield(D, 'u_sel') && ~isempty(D.u_sel)
            X = D.u_sel;
            evalues = local_first_field(D, {'evalues_ref_sel', 'evalues_sel', 'evalues_used_sel'});
            raw_field = 'EDMD_outputs.deconv_efuns.u_sel';
            kpm = local_kpm_from_selected(E, size(X, 2));
            original_idx = local_original_indices_from_selected(E, size(X, 2));
        elseif isfield(D, 'u_all') && ~isempty(D.u_all)
            X = D.u_all;
            evalues = local_first_field(D, {'evalues_ref_all', 'evalues_all', 'evalues_used_all'});
            raw_field = 'EDMD_outputs.deconv_efuns.u_all';
            kpm = local_kpm_from_all(E, size(X, 2));
            original_idx = (1:size(X, 2)).';
        else
            error('EDMD_outputs.deconv_efuns must contain u_sel or u_all.');
        end
        if isempty(evalues)
            evalues = local_fallback_evalues(E, size(X, 2));
        end

    otherwise
        error('Unsupported source kind: %s.', source_kind);
end

evalues = evalues(:);
if numel(evalues) < size(X, 2)
    error('Not enough evalues for %s: %d values for %d columns.', ...
        feature_name, numel(evalues), size(X, 2));
end
if numel(evalues) > size(X, 2)
    evalues = evalues(1:size(X, 2));
end
if numel(original_idx) < size(X, 2)
    original_idx = (1:size(X, 2)).';
end
original_idx = double(original_idx(1:size(X, 2)));

info = struct();
info.name = feature_name;
info.source_kind = source_kind;
info.variant = variant;
info.normalized_variant = variant;
info.raw_field = raw_field;
info.axis_order = 'time_by_mode';
end


function values = local_first_field(S, names)
values = [];
for i = 1:numel(names)
    name = names{i};
    if isfield(S, name) && ~isempty(S.(name))
        values = S.(name)(:);
        return;
    end
end
end


function evalues = local_fallback_evalues(E, n)
if isfield(E, 'evalues') && numel(E.evalues) >= n
    evalues = E.evalues(1:n);
elseif isfield(E, 'original_sorted') && isfield(E.original_sorted, 'evalues') && ...
        numel(E.original_sorted.evalues) >= n
    evalues = E.original_sorted.evalues(1:n);
else
    error('Could not resolve evalues for %d BOLD feature columns.', n);
end
evalues = evalues(:);
end


function kpm = local_kpm_from_selected(E, n)
kpm = [];
if isfield(E, 'kpm_modes') && size(E.kpm_modes, 1) >= n
    kpm = E.kpm_modes(1:n, :);
end
end


function kpm = local_kpm_from_all(E, n)
kpm = [];
if isfield(E, 'original_sorted') && isfield(E.original_sorted, 'kpm_modes') && ...
        size(E.original_sorted.kpm_modes, 1) >= n
    kpm = E.original_sorted.kpm_modes(1:n, :);
elseif isfield(E, 'kpm_modes') && size(E.kpm_modes, 1) >= n
    kpm = E.kpm_modes(1:n, :);
end
end


function idx = local_original_indices_from_selected(E, n)
if isfield(E, 'idx_final_in_original') && numel(E.idx_final_in_original) >= n
    idx = double(E.idx_final_in_original(1:n));
else
    idx = (1:n).';
end
idx = idx(:);
end


function ord = local_sort_modes(evalues, sort_by, sort_dir)
switch lower(char(string(sort_by)))
    case {'preserve', 'none', 'input'}
        ord = (1:numel(evalues)).';
    case {'modulus', 'abs'}
        [~, ord] = sort(abs(evalues), sort_dir);
    case {'real', 'realpart'}
        [~, ord] = sort(real(evalues), sort_dir);
    otherwise
        error('Unknown selection.sort_by = %s.', sort_by);
end
end


function [dt, dt_source] = local_resolve_dt(B, E)
dt = [];
dt_source = struct('source', 'missing', 'field', '', 'message', '');
fields = {'dt', 'dx', 'sampling_period', 'sample_period'};
for i = 1:numel(fields)
    name = fields{i};
    if isfield(B, name) && ~isempty(B.(name))
        dt = local_positive_scalar(B.(name));
        if ~isempty(dt)
            dt_source.source = 'BOLD_POST';
            dt_source.field = name;
            return;
        end
    end
    if isfield(E, name) && ~isempty(E.(name))
        dt = local_positive_scalar(E.(name));
        if ~isempty(dt)
            dt_source.source = 'EDMD_outputs';
            dt_source.field = name;
            return;
        end
    end
end
dt_source.message = 'No positive scalar dt/dx was found in BOLD_POST or EDMD_outputs.';
end


function value = local_positive_scalar(value_in)
value = [];
if ~isnumeric(value_in) || isempty(value_in)
    return;
end
value_in = double(real(value_in));
value_in = value_in(:);
value_in = value_in(isfinite(value_in) & value_in > 0);
if ~isempty(value_in)
    value = value_in(1);
end
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
