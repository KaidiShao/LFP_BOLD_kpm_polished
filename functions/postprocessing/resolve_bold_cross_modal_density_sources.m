function density_sources = resolve_bold_cross_modal_density_sources(dataset_ref, params)
%RESOLVE_BOLD_CROSS_MODAL_DENSITY_SOURCES Resolve canonical pipeline 8 density inputs.
%
% Default pipeline 8 behavior is to pair each pipeline 7 BOLD_POST run with
% explicit, non-ambiguous density sources:
%   1) pipeline 2 event density
%   2) BLP abs raw thresholded density at one threshold
%   3) BLP abs dimred thresholded density at the same threshold
%   4) BLP complex-split raw thresholded density at the same threshold
%   5) BLP complex-split dimred thresholded density at the same threshold

if nargin < 1 || isempty(dataset_ref)
    error('dataset_ref is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

[dataset_stem, dataset_id] = local_resolve_dataset_identity(dataset_ref);
requested_kinds = local_normalize_density_kinds(params.density_source_kinds);

density_sources = repmat(local_empty_density_source(), 0, 1);
pipeline5_cache = struct('loaded', false, 'result_file', '', 'result', struct());

for i_kind = 1:numel(requested_kinds)
    kind = requested_kinds{i_kind};
    try
        switch kind
            case 'event_density'
                source = local_resolve_event_density_source( ...
                    params.processed_root, dataset_stem, dataset_id);
            case 'raw_abs_density'
                source = local_resolve_explicit_pipeline5_density_source( ...
                    params.processed_root, dataset_stem, dataset_id, ...
                    'raw_thresholded_density', 'abs_projected_vlambda', '', ...
                    params.blp_density_threshold_ratio, ...
                    local_explicit_source_name('raw_abs', '', params.blp_density_threshold_ratio), ...
                    'thresholded_density');
            case 'dimred_abs_density'
                source = local_resolve_explicit_pipeline5_density_source( ...
                    params.processed_root, dataset_stem, dataset_id, ...
                    'dimred_thresholded_density', 'abs_projected_vlambda', ...
                    params.blp_dimred_method_tag, params.blp_density_threshold_ratio, ...
                    local_explicit_source_name('dim_abs', params.blp_dimred_method_tag, ...
                    params.blp_density_threshold_ratio), 'dimred_thresholded_density');
            case 'raw_complex_split_density'
                source = local_resolve_explicit_pipeline5_density_source( ...
                    params.processed_root, dataset_stem, dataset_id, ...
                    'raw_thresholded_density', 'complex_split_projected_vlambda', '', ...
                    params.blp_density_threshold_ratio, ...
                    local_explicit_source_name('raw_csplit', '', params.blp_density_threshold_ratio), ...
                    'thresholded_density');
            case 'dimred_complex_split_density'
                source = local_resolve_explicit_pipeline5_density_source( ...
                    params.processed_root, dataset_stem, dataset_id, ...
                    'dimred_thresholded_density', 'complex_split_projected_vlambda', ...
                    params.blp_dimred_method_tag, params.blp_density_threshold_ratio, ...
                    local_explicit_source_name('dim_csplit', params.blp_dimred_method_tag, ...
                    params.blp_density_threshold_ratio), 'dimred_thresholded_density');
            case 'raw_eigenfunction_density'
                [source, pipeline5_cache] = local_resolve_pipeline5_density_source( ...
                    params.processed_root, dataset_stem, dataset_id, ...
                    'thresholded_density_mat_file', ...
                    'blp_raw_eigenfunction_density', ...
                    'thresholded_density', pipeline5_cache);
            case 'dimred_eigenfunction_density'
                [source, pipeline5_cache] = local_resolve_pipeline5_density_source( ...
                    params.processed_root, dataset_stem, dataset_id, ...
                    'dimred_thresholded_density_mat_file', ...
                    'blp_dimred_eigenfunction_density', ...
                    'dimred_thresholded_density', pipeline5_cache);
            otherwise
                error('Unsupported density source kind: %s', kind);
        end
        density_sources(end + 1, 1) = source; %#ok<AGROW>
    catch ME
        if params.require_all_density_sources
            error(['Failed to resolve density source "%s" for dataset %s.\n' ...
                '%s'], kind, dataset_stem, ME.message);
        end
        warning('Skipping density source %s for %s: %s', ...
            kind, dataset_stem, ME.message);
    end
end

if isempty(density_sources)
    error('No pipeline 8 density sources were resolved for dataset %s.', dataset_stem);
end
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'processed_root', io_project.get_project_processed_root());
params = local_set_default(params, 'density_source_kinds', { ...
    'event_density', ...
    'raw_abs_density', ...
    'dimred_abs_density', ...
    'raw_complex_split_density', ...
    'dimred_complex_split_density'});
params = local_set_default(params, 'require_all_density_sources', false);
params = local_set_default(params, 'blp_density_threshold_ratio', 0.70);
params = local_set_default(params, 'blp_dimred_method_tag', 'umap_k08');
end


function [dataset_stem, dataset_id] = local_resolve_dataset_identity(dataset_ref)
dataset_stem = '';
dataset_id = '';

if ischar(dataset_ref) || isstring(dataset_ref)
    dataset_stem = char(string(dataset_ref));
elseif isstruct(dataset_ref)
    if isfield(dataset_ref, 'dataset_stem') && ~isempty(dataset_ref.dataset_stem)
        dataset_stem = char(string(dataset_ref.dataset_stem));
    elseif isfield(dataset_ref, 'file_stem') && ~isempty(dataset_ref.file_stem)
        dataset_stem = char(string(dataset_ref.file_stem));
    end

    if isfield(dataset_ref, 'dataset_id') && ~isempty(dataset_ref.dataset_id)
        dataset_id = char(string(dataset_ref.dataset_id));
    end
else
    error('dataset_ref must be a dataset stem string or a struct-like ref.');
end

if isempty(dataset_stem)
    error('Could not resolve dataset_stem from dataset_ref.');
end
if isempty(dataset_id)
    dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
end
end


function kinds = local_normalize_density_kinds(kinds)
if isempty(kinds)
    kinds = {};
elseif ischar(kinds) || isstring(kinds)
    kinds = cellstr(string(kinds(:)).');
elseif iscell(kinds)
    kinds = cellstr(string(kinds(:)).');
else
    error('density_source_kinds must be a string-like list.');
end

out = {};
for i = 1:numel(kinds)
    key = lower(char(string(kinds{i})));
    switch key
        case {'event', 'event_density', 'blp_event_density', 'evt'}
            out = [out, {'event_density'}]; %#ok<AGROW>
        case {'raw_abs', 'raw_abs_density', 'raw_abs_q070', ...
                'abs_raw', 'abs_raw_density'}
            out = [out, {'raw_abs_density'}]; %#ok<AGROW>
        case {'dimred_abs', 'dimred_abs_density', 'dim_abs', ...
                'dim_abs_umap8_q070', 'abs_dimred', 'abs_dimred_density'}
            out = [out, {'dimred_abs_density'}]; %#ok<AGROW>
        case {'raw_complex_split', 'raw_complex_split_density', ...
                'raw_csplit', 'raw_csplit_q070', 'complex_split_raw'}
            out = [out, {'raw_complex_split_density'}]; %#ok<AGROW>
        case {'dimred_complex_split', 'dimred_complex_split_density', ...
                'dim_csplit', 'dim_csplit_umap8_q070', 'complex_split_dimred'}
            out = [out, {'dimred_complex_split_density'}]; %#ok<AGROW>
        case {'raw', 'raw_density', 'raw_thresholded_density', ...
                'raw_eigenfunction_density'}
            out = [out, {'raw_abs_density', 'raw_complex_split_density'}]; %#ok<AGROW>
        case {'dimred', 'dimred_density', 'dimred_thresholded_density', ...
                'dimred_eigenfunction_density'}
            out = [out, {'dimred_abs_density', 'dimred_complex_split_density'}]; %#ok<AGROW>
        case {'legacy_raw_eigenfunction_density'}
            out = [out, {'raw_eigenfunction_density'}]; %#ok<AGROW>
        case {'legacy_dimred_eigenfunction_density'}
            out = [out, {'dimred_eigenfunction_density'}]; %#ok<AGROW>
        otherwise
            error('Unsupported density source kind: %s', kinds{i});
    end
end
kinds = unique(out, 'stable');
end


function source = local_resolve_event_density_source(processed_root, dataset_stem, dataset_id)
event_dir = io_project.get_pipeline_stage_dir( ...
    processed_root, dataset_stem, 2, 'event_density');
event_file = local_find_event_density_file(event_dir, dataset_stem, dataset_id);

source = local_empty_density_source();
source.name = 'blp_evt';
source.short_name = 'blp_evt';
source.type = 'event_density';
source.file = event_file;
source.dataset_stem = dataset_stem;
source.dataset_id = dataset_id;
source.stage_name = io_project.get_pipeline_stage_name(2, 'event_density');
end


function event_file = local_find_event_density_file(event_dir, dataset_stem, dataset_id)
if exist(event_dir, 'dir') ~= 7
    error('Pipeline 2 event density directory does not exist: %s', event_dir);
end

candidate_names = { ...
    sprintf('%s_event_density_2s.mat', dataset_stem), ...
    sprintf('%s_event_density_2s.mat', dataset_id)};

for i = 1:numel(candidate_names)
    candidate = fullfile(event_dir, candidate_names{i});
    if exist(candidate, 'file') == 2
        event_file = candidate;
        return;
    end
end

L = dir(fullfile(event_dir, '*_event_density_2s.mat'));
L = L(~[L.isdir]);
if isempty(L)
    error(['No pipeline 2 event density MAT matching *_event_density_2s.mat ' ...
        'was found under: %s'], event_dir);
end
if numel(L) > 1
    names = strjoin({L.name}, ', ');
    error(['Multiple event density MAT files were found under %s. ' ...
        'Expected one canonical file. Candidates: %s'], event_dir, names);
end

event_file = fullfile(L(1).folder, L(1).name);
end


function source = local_resolve_explicit_pipeline5_density_source( ...
        processed_root, dataset_stem, dataset_id, stage_key, condition_tag, ...
        method_tag, threshold_ratio, source_name, source_type)
stage_root = io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 5, stage_key);
threshold_tag = local_threshold_tag(threshold_ratio);
artifact_file = local_find_explicit_pipeline5_density_file( ...
    stage_root, dataset_stem, condition_tag, method_tag, threshold_tag);
if isempty(artifact_file)
    error(['No explicit pipeline 5 density artifact was found for %s.\n' ...
        'Stage root: %s\nCondition: %s\nMethod: %s\nThreshold: %s'], ...
        source_name, stage_root, condition_tag, method_tag, threshold_tag);
end

source = local_empty_density_source();
source.name = source_name;
source.type = source_type;
source.file = artifact_file;
source.dataset_stem = dataset_stem;
source.dataset_id = dataset_id;
source.stage_name = io_project.get_pipeline_stage_name(5, stage_key);
source.condition_tag = condition_tag;
source.feature_variant = local_condition_feature_variant(condition_tag);
source.threshold_ratio = threshold_ratio;
source.threshold_tag = threshold_tag;
source.dimred_method_tag = method_tag;
source.short_name = source_name;
end


function artifact_file = local_find_explicit_pipeline5_density_file( ...
        stage_root, dataset_stem, condition_tag, method_tag, threshold_tag)
artifact_file = '';
if exist(stage_root, 'dir') ~= 7
    return;
end

if isempty(method_tag)
    search_root = fullfile(stage_root, condition_tag, 'mat', threshold_tag);
    pattern = sprintf('%s_*ratio_%s*%s*.mat', ...
        dataset_stem, threshold_tag(2:end), condition_tag);
else
    search_root = fullfile(stage_root, condition_tag, method_tag, 'mat', threshold_tag);
    pattern = sprintf('%s_*ratio_%s*%s_%s*.mat', ...
        dataset_stem, threshold_tag(2:end), condition_tag, method_tag);
end

L = dir(fullfile(search_root, pattern));
if isempty(L)
    if isempty(method_tag)
        fallback = fullfile(stage_root, condition_tag, 'mat');
        L = dir(fullfile(fallback, '**', pattern));
    else
        fallback = fullfile(stage_root, condition_tag, method_tag, 'mat');
        L = dir(fullfile(fallback, '**', pattern));
    end
end
L = L(~[L.isdir]);
if isempty(L)
    return;
end
[~, idx] = max([L.datenum]);
artifact_file = fullfile(L(idx).folder, L(idx).name);
end


function tag = local_threshold_tag(threshold_ratio)
value = round(double(threshold_ratio) * 100);
tag = sprintf('q%03d', value);
end


function variant = local_condition_feature_variant(condition_tag)
condition = lower(char(string(condition_tag)));
if contains(condition, 'complex_split')
    variant = 'complex_split';
elseif contains(condition, 'abs')
    variant = 'abs';
else
    variant = condition;
end
end


function name = local_explicit_source_name(kind, method_tag, threshold_ratio)
q = local_threshold_tag(threshold_ratio);
switch lower(char(string(kind)))
    case {'raw_abs', 'raw_csplit'}
        name = sprintf('%s_%s', char(string(kind)), q);
    otherwise
        method_token = local_method_token(method_tag);
        if isempty(method_token)
            name = sprintf('%s_%s', char(string(kind)), q);
        else
            name = sprintf('%s_%s_%s', char(string(kind)), method_token, q);
        end
end
end


function token = local_method_token(method_tag)
token = lower(char(string(method_tag)));
token = regexprep(token, '_k0*', '');
token = regexprep(token, '[^a-z0-9]+', '_');
token = regexprep(token, '^_+|_+$', '');
end


function [source, cache] = local_resolve_pipeline5_density_source( ...
        processed_root, dataset_stem, dataset_id, artifact_field, ...
        source_name, source_type, cache)
cache = local_ensure_pipeline5_result_cache(cache, processed_root, dataset_stem);
stage_key = local_pipeline5_stage_key(source_type);
artifact_file = local_resolve_pipeline5_density_artifact( ...
    cache, processed_root, dataset_stem, artifact_field, stage_key);

source = local_empty_density_source();
source.name = source_name;
source.type = source_type;
source.file = artifact_file;
source.dataset_stem = dataset_stem;
source.dataset_id = dataset_id;
source.stage_name = io_project.get_pipeline_stage_name(5, 'eigenfunction_reduction');
source.pipeline5_result_file = cache.result_file;
end


function stage_key = local_pipeline5_stage_key(source_type)
switch lower(char(string(source_type)))
    case 'thresholded_density'
        stage_key = 'raw_thresholded_density';
    case 'dimred_thresholded_density'
        stage_key = 'dimred_thresholded_density';
    otherwise
        error('Unsupported pipeline 5 density source_type: %s', source_type);
end
end


function artifact_file = local_resolve_pipeline5_density_artifact( ...
        cache, processed_root, dataset_stem, artifact_field, stage_key)
artifact_file = '';
legacy_artifact_file = '';
if isfield(cache.result, 'artifacts') && isstruct(cache.result.artifacts) && ...
        isfield(cache.result.artifacts, artifact_field) && ...
        ~isempty(cache.result.artifacts.(artifact_field))
    legacy_artifact_file = char(string(cache.result.artifacts.(artifact_field)));
    if exist(legacy_artifact_file, 'file') == 2
        artifact_file = legacy_artifact_file;
        return;
    end
end

[condition_tag, method_tag] = local_parse_pipeline5_result_location(cache.result_file);
stage_root = io_project.get_pipeline_stage_dir(processed_root, dataset_stem, 5, stage_key);

artifact_file = local_find_existing_pipeline5_density_artifact( ...
    stage_root, legacy_artifact_file, dataset_stem, condition_tag);
if ~isempty(artifact_file)
    return;
end
error('%s', local_missing_pipeline5_density_message( ...
    stage_root, stage_key, condition_tag, method_tag, legacy_artifact_file));
end


function [condition_tag, method_tag] = local_parse_pipeline5_result_location(result_file)
condition_tag = '';
method_tag = '';
if isempty(result_file)
    return;
end

result_mat_dir = fileparts(result_file);
method_root = fileparts(result_mat_dir);
condition_root = fileparts(method_root);
if ~isempty(method_root)
    [~, method_tag] = fileparts(method_root);
end
if ~isempty(condition_root)
    [~, condition_tag] = fileparts(condition_root);
end
end


function artifact_file = local_find_existing_pipeline5_density_artifact( ...
        stage_root, legacy_artifact_file, dataset_stem, condition_tag)
artifact_file = '';
if exist(stage_root, 'dir') ~= 7
    return;
end

candidate_files = {};
if ~isempty(legacy_artifact_file)
    [~, base_name, ext] = fileparts(legacy_artifact_file);
    exact_hits = dir(fullfile(stage_root, '**', [base_name ext]));
    candidate_files = [candidate_files; local_fullpaths(exact_hits)]; %#ok<AGROW>
end

if isempty(candidate_files) && ~isempty(condition_tag)
    pattern = sprintf('%s_*%s*.mat', dataset_stem, condition_tag);
    condition_hits = dir(fullfile(stage_root, '**', pattern));
    candidate_files = [candidate_files; local_fullpaths(condition_hits)]; %#ok<AGROW>
end

if isempty(candidate_files)
    return;
end

candidate_files = unique(candidate_files, 'stable');
if numel(candidate_files) == 1
    artifact_file = candidate_files{1};
    return;
end

times = nan(numel(candidate_files), 1);
for i = 1:numel(candidate_files)
    info = dir(candidate_files{i});
    if ~isempty(info)
        times(i) = info(1).datenum;
    end
end
[~, best_idx] = max(times);
artifact_file = candidate_files{best_idx};
end


function paths = local_fullpaths(L)
paths = cell(0, 1);
if isempty(L)
    return;
end
L = L(~[L.isdir]);
paths = arrayfun(@(d) fullfile(d.folder, d.name), L, 'UniformOutput', false);
paths = paths(:);
end


function msg = local_missing_pipeline5_density_message( ...
        stage_root, stage_key, condition_tag, method_tag, legacy_artifact_file)
parts = {sprintf('No existing pipeline 5 %s artifact was found.', stage_key)};
parts{end + 1} = sprintf('Searched stage root: %s', stage_root); %#ok<AGROW>
if ~isempty(condition_tag)
    parts{end + 1} = sprintf('Condition tag: %s', condition_tag); %#ok<AGROW>
end
if ~isempty(method_tag)
    parts{end + 1} = sprintf('Method tag: %s', method_tag); %#ok<AGROW>
end
if ~isempty(legacy_artifact_file)
    parts{end + 1} = sprintf('Legacy artifact path was missing: %s', legacy_artifact_file); %#ok<AGROW>
end
parts{end + 1} = 'Pipeline 8 will not backfill pipeline 5 outputs; run or restore the corresponding pipeline 5 density stage first if you need this source.'; %#ok<AGROW>
msg = strjoin(parts, newline);
end


function cache = local_ensure_pipeline5_result_cache(cache, processed_root, dataset_stem)
if isfield(cache, 'loaded') && cache.loaded
    return;
end

result_file = find_latest_blp_eigenfunction_reduction_result(processed_root, dataset_stem);
S = load(result_file, 'result');
if ~isfield(S, 'result')
    error('Pipeline 5 result variable "result" is missing in %s.', result_file);
end

cache.loaded = true;
cache.result_file = result_file;
cache.result = S.result;
end


function source = local_empty_density_source()
source = struct( ...
    'name', '', ...
    'short_name', '', ...
    'type', '', ...
    'file', '', ...
    'dataset_stem', '', ...
    'dataset_id', '', ...
    'stage_name', '', ...
    'pipeline5_result_file', '', ...
    'condition_tag', '', ...
    'feature_variant', '', ...
    'threshold_ratio', NaN, ...
    'threshold_tag', '', ...
    'dimred_method_tag', '');
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end
