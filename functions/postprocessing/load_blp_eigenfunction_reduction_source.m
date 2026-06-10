function [EDMD_outputs, concat_info, source_info] = load_blp_eigenfunction_reduction_source(run_info, params)
%LOAD_BLP_EIGENFUNCTION_REDUCTION_SOURCE Load one completed run once for reuse across methods.

source_cfg = struct();
source_cfg.mode = 'chunk_dir';
source_cfg.data_dir = run_info.output_dir;
source_cfg.concat = struct();
source_cfg.concat.filename_pattern = params.filename_pattern;
source_cfg.concat.variable_name = params.variable_name;
source_cfg.concat.concat_fields = {'efuns'};
source_cfg.concat.concat_dim = 1;
source_cfg.concat.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form', ...
    'dt', 'dx', 'sampling_period', 'sample_period', 'fs', 'sampling_frequency', ...
    'session_dx', 'session_fs', 'session_ids', 'session_lengths', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
source_cfg.concat.allow_missing_chunks = false;
source_cfg.concat.verbose = false;

mode_indices = local_select_source_mode_indices(run_info.output_dir, params);
if ~isempty(mode_indices)
    source_cfg.concat.mode_indices = mode_indices;
    source_cfg.concat.concat_column_indices = struct('efuns', mode_indices);
end

allow_summary_only = isfield(params, 'allow_summary_only_outputs') && ...
    logical(params.allow_summary_only_outputs);
if isempty(dir(fullfile(run_info.output_dir, params.filename_pattern))) && ...
        allow_summary_only
    [EDMD_outputs, concat_info, source_info] = ...
        local_load_summary_only_edmd_outputs(run_info.output_dir, params);
else
    [EDMD_outputs, concat_info, source_info] = io_edmd.load_edmd_source(source_cfg);
end
if ~isempty(mode_indices)
    source_info.mode_indices_in_original = mode_indices(:);
end
end


function [EDMD_outputs, concat_info, source_info] = local_load_summary_only_edmd_outputs(output_dir, params)
summary_files = dir(fullfile(output_dir, '*_summary.mat'));
summary_files = summary_files(~[summary_files.isdir]);
if isempty(summary_files)
    error('No summary MAT files found for summary-only BLP output:\n  %s', output_dir);
end

[~, order] = sort([summary_files.datenum], 'descend');
summary_file = fullfile(summary_files(order(1)).folder, summary_files(order(1)).name);
vars = who('-file', summary_file);
if ~any(strcmp(vars, params.variable_name))
    error('Summary MAT lacks variable %s:\n  %s', params.variable_name, summary_file);
end

S = load(summary_file, params.variable_name);
EDMD_outputs = S.(params.variable_name);
concat_info = struct();
concat_info.data_dir = output_dir;
concat_info.file_pattern = '*_summary.mat';
concat_info.variable_name = params.variable_name;
concat_info.files = {summary_file};
if isfield(EDMD_outputs, 'efuns')
    concat_info.total_length = size(EDMD_outputs.efuns, 1);
else
    concat_info.total_length = NaN;
end
concat_info.n_chunks = 0;
concat_info.summary_only = true;

source_info = struct();
source_info.mode = 'summary_file';
source_info.path = summary_file;
source_info.output_dir = output_dir;
source_info.summary_only = true;
end


function mode_indices = local_select_source_mode_indices(output_dir, params)
mode_indices = [];

files = dir(fullfile(output_dir, params.filename_pattern));
files = files(~[files.isdir]);
if isempty(files)
    return;
end

chunk_ids = nan(numel(files), 1);
for i = 1:numel(files)
    tokens = regexp(files(i).name, '^(.*)_outputs_(\d+)\.mat$', 'tokens', 'once');
    if ~isempty(tokens)
        chunk_ids(i) = str2double(tokens{2});
    end
end
files = files(isfinite(chunk_ids));
chunk_ids = chunk_ids(isfinite(chunk_ids));
if isempty(files)
    return;
end

[~, order] = sort(chunk_ids);
first_file = fullfile(files(order(1)).folder, files(order(1)).name);
S = load(first_file, params.variable_name);
if ~isfield(S, params.variable_name)
    return;
end

EDMD_outputs = S.(params.variable_name);
if ~isfield(EDMD_outputs, 'evalues') || isempty(EDMD_outputs.evalues)
    return;
end

defaults = build_blp_eigenfunction_reduction_defaults();
selection = defaults.selection;
evalues = EDMD_outputs.evalues(:);

switch lower(selection.sort_by)
    case {'modulus', 'abs'}
        key = abs(evalues);
    case {'real', 'realpart'}
        key = real(evalues);
    otherwise
        key = abs(evalues);
end

[~, order] = sort(key, selection.sort_dir);
evalues_sorted = evalues(order);
keep_sorted = abs(evalues_sorted) > selection.abs_thresh;
mode_indices = order(keep_sorted);

if isfinite(selection.max_modes)
    mode_indices = mode_indices(1:min(numel(mode_indices), round(selection.max_modes)));
end

mode_indices = mode_indices(:).';
end
