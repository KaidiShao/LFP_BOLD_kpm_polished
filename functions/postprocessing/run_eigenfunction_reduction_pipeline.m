function [result, EDMD_outputs, concat_info, source_info] = run_eigenfunction_reduction_pipeline(cfg)
%RUN_EIGENFUNCTION_REDUCTION_PIPELINE Minimal eigenfunction-only reduction pipeline.

if nargin < 1 || isempty(cfg)
    error('cfg must be provided.');
end

cfg = local_apply_defaults(cfg);
[EDMD_outputs, concat_info, source_info] = load_edmd_source(cfg.source);
prep = local_prepare_eigenfunction_inputs(EDMD_outputs, cfg);

switch lower(cfg.path.kind)
    case 'time'
        reduction = reduce_eigenfunction_time_path( ...
            prep.efun_feature_time_by_mode, cfg.path.time);
        method_name = cfg.path.time.method;

    case 'spectrum'
        spectrum_cfg = cfg.path.spectrum;
        spectrum_cfg.evalues_discrete = prep.evalues_discrete;
        spectrum_cfg.evalues_bilinear = prep.evalues_bilinear;
        spectrum_cfg.feature_variant = cfg.feature.variant;
        reduction = reduce_eigenfunction_spectrum_path( ...
            prep.efun_feature_time_by_mode, spectrum_cfg);
        method_name = cfg.path.spectrum.method;

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

result.cfg = cfg;
result.source = source_info;
result.concat = concat_info;

result.input = struct();
result.input.dt = prep.dt;
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
result.summary = local_build_summary(result.core, cfg.summary);

result.artifacts = struct();
result.artifacts.result_mat_file = '';

if cfg.save.enable
    if exist(cfg.save.dir, 'dir') ~= 7
        mkdir(cfg.save.dir);
    end

    save_path = local_build_save_path(cfg, result);
    result.artifacts.result_mat_file = save_path;

    if cfg.save.v7_3
        save(save_path, 'result', '-v7.3');
    else
        save(save_path, 'result');
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
    cfg.save.dir = fullfile(get_project_results_root(cfg.repo_root), ...
        'eigenfunction_reduction');
end

if ~isfield(cfg.save, 'file_stem') || isempty(cfg.save.file_stem)
    cfg.save.file_stem = 'eigenfunction_reduction_result';
end

if ~isfield(cfg.save, 'tag')
    cfg.save.tag = '';
end

if ~isfield(cfg.save, 'v7_3') || isempty(cfg.save.v7_3)
    cfg.save.v7_3 = true;
end
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

dt = local_resolve_dt(cfg.input, EDMD_outputs);
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


function dt = local_resolve_dt(input_cfg, EDMD_outputs)
dt = [];

if isfield(input_cfg, 'dt') && ~isempty(input_cfg.dt)
    dt = input_cfg.dt;
    return;
end

candidate_fields = {'dt', 'dx', 'sampling_period', 'sample_period'};
for i = 1:numel(candidate_fields)
    field_name = candidate_fields{i};
    if isfield(EDMD_outputs, field_name)
        value = EDMD_outputs.(field_name);
        if isnumeric(value) && isscalar(value) && isfinite(value) && value > 0
            dt = value;
            return;
        end
    end
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


function save_path = local_build_save_path(cfg, result)
timestamp = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
pieces = {cfg.save.file_stem, cfg.path.kind, lower(cfg.feature.variant), lower(result.meta.method)};

if isfield(cfg.save, 'tag') && ~isempty(cfg.save.tag)
    pieces{end+1} = cfg.save.tag; %#ok<AGROW>
end

filename = strjoin(pieces, '__');
filename = regexprep(filename, '[^\w\-]+', '_');
save_path = fullfile(cfg.save.dir, sprintf('%s__%s.mat', filename, timestamp));
end
