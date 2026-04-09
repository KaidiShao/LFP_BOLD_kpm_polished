function [EDMD_outputs, concat_info] = load_and_concat_edmd_output_chunks(data_dir, params)
% Load EDMD output chunk files and concatenate selected fields in time.
%
% Example
%   params = struct();
%   params.concat_fields = {'efuns'};
%   [EDMD_outputs, concat_info] = load_and_concat_edmd_output_chunks( ...
%       'E:\autodl_results\initial_point_test1_KV', params);
%
% Inputs
%   data_dir Folder that contains files named like *_outputs_<chunk>.mat
%   params   Optional struct with fields:
%            .filename_pattern       (default '*_outputs_*.mat')
%            .variable_name          (default 'EDMD_outputs')
%            .concat_fields          (default {'efuns'})
%            .concat_dim             (default 1)
%            .required_equal_fields  (default {'evalues','kpm_modes',...
%                                        'N_dict','residual_form'})
%            .allow_missing_chunks   (default false)
%            .verbose                (default true)
%            .progress_every         (default 50)
%
% Outputs
%   EDMD_outputs Concatenated output struct
%   concat_info  Metadata for the concatenation

if nargin < 1 || isempty(data_dir)
    error('data_dir must be provided.');
end

if nargin < 2
    params = struct();
end

params = apply_default_params(params);
files = collect_chunk_files(data_dir, params.filename_pattern);

if isempty(files)
    error('No files matching %s were found in %s.', params.filename_pattern, data_dir);
end

chunk_ids = [files.chunk_id]';
check_chunk_sequence(chunk_ids, params.allow_missing_chunks);

if params.verbose
    fprintf('Found %d EDMD chunk files in %s\n', numel(files), data_dir);
    fprintf('Chunk range: %d -> %d\n', chunk_ids(1), chunk_ids(end));
end

first_data = load_variable(files(1).fullpath, params.variable_name);
EDMD_outputs = first_data.(params.variable_name);

all_fields = fieldnames(EDMD_outputs);
concat_fields = validate_field_list(params.concat_fields, all_fields, 'concat_fields');
copied_fields = setdiff(all_fields, concat_fields, 'stable');
equal_fields = intersect(params.required_equal_fields, copied_fields, 'stable');

n_chunks = numel(files);
n_concat_fields = numel(concat_fields);

chunk_lengths = zeros(n_chunks, 1);
concat_lengths_by_field = zeros(n_chunks, n_concat_fields);
ref_values = cell(n_concat_fields, 1);
ref_sizes = cell(n_concat_fields, 1);
final_sizes = cell(n_concat_fields, 1);
estimated_bytes = zeros(n_concat_fields, 1);

for i = 1:n_concat_fields
    field_name = concat_fields{i};
    ref_values{i} = EDMD_outputs.(field_name);
    ref_sizes{i} = size(ref_values{i});
    validate_concat_value(ref_values{i}, field_name);
    concat_lengths_by_field(1, i) = get_concat_length(ref_values{i}, params.concat_dim);
end

chunk_lengths(1) = validate_chunk_lengths(concat_lengths_by_field(1, :), files(1).name);

for k = 2:n_chunks
    if should_print_progress(k, n_chunks, params.progress_every) && params.verbose
        fprintf('[scan %d/%d] %s\n', k, n_chunks, files(k).name);
    end

    current_data = load_variable(files(k).fullpath, params.variable_name);
    current_outputs = current_data.(params.variable_name);

    check_equal_fields(EDMD_outputs, current_outputs, equal_fields, files(k).name);

    for i = 1:n_concat_fields
        field_name = concat_fields{i};
        current_value = read_required_field(current_outputs, field_name, files(k).name);
        validate_concat_value(current_value, field_name);
        validate_concat_shape(ref_sizes{i}, size(current_value), params.concat_dim, field_name, files(k).name);
        concat_lengths_by_field(k, i) = get_concat_length(current_value, params.concat_dim);
    end

    chunk_lengths(k) = validate_chunk_lengths(concat_lengths_by_field(k, :), files(k).name);
end

total_length = sum(chunk_lengths);

for i = 1:n_concat_fields
    final_sizes{i} = build_final_size(ref_sizes{i}, params.concat_dim, sum(concat_lengths_by_field(:, i)));
    estimated_bytes(i) = prod(double(final_sizes{i})) * class_nbytes(class(ref_values{i}));
end

if params.verbose
    fprintf('Estimated concatenated array size: %.3f GiB\n', sum(estimated_bytes) / 1024^3);
end

for i = 1:n_concat_fields
    field_name = concat_fields{i};
    EDMD_outputs.(field_name) = zeros(final_sizes{i}, 'like', ref_values{i});
end

chunk_start_idx = zeros(n_chunks, 1);
chunk_end_idx = zeros(n_chunks, 1);
cursor = 1;

for k = 1:n_chunks
    if should_print_progress(k, n_chunks, params.progress_every) && params.verbose
        fprintf('[load %d/%d] %s\n', k, n_chunks, files(k).name);
    end

    if k == 1
        current_outputs = first_data.(params.variable_name);
    else
        current_data = load_variable(files(k).fullpath, params.variable_name);
        current_outputs = current_data.(params.variable_name);
    end

    chunk_start_idx(k) = cursor;
    chunk_end_idx(k) = cursor + chunk_lengths(k) - 1;

    for i = 1:n_concat_fields
        field_name = concat_fields{i};
        EDMD_outputs.(field_name) = assign_concat_block( ...
            EDMD_outputs.(field_name), ...
            current_outputs.(field_name), ...
            params.concat_dim, ...
            cursor);
    end

    cursor = chunk_end_idx(k) + 1;
end

concat_info = struct();
concat_info.data_dir = data_dir;
concat_info.file_pattern = params.filename_pattern;
concat_info.variable_name = params.variable_name;
concat_info.file_prefix = files(1).prefix;
concat_info.n_chunks = n_chunks;
concat_info.chunk_ids = chunk_ids;
concat_info.files = {files.fullpath}';
concat_info.chunk_lengths = chunk_lengths;
concat_info.chunk_start_idx = chunk_start_idx;
concat_info.chunk_end_idx = chunk_end_idx;
concat_info.total_length = total_length;
concat_info.concat_fields = concat_fields(:);
concat_info.concat_dim = params.concat_dim;
concat_info.required_equal_fields = equal_fields(:);
concat_info.copied_from_first_fields = copied_fields(:);
concat_info.estimated_concat_bytes = sum(estimated_bytes);
concat_info.estimated_concat_gib = sum(estimated_bytes) / 1024^3;

for i = 1:n_concat_fields
    field_name = concat_fields{i};
    concat_info.final_sizes.(field_name) = final_sizes{i};
end
end


function params = apply_default_params(params)
if ~isfield(params, 'filename_pattern') || isempty(params.filename_pattern)
    params.filename_pattern = '*_outputs_*.mat';
end

if ~isfield(params, 'variable_name') || isempty(params.variable_name)
    params.variable_name = 'EDMD_outputs';
end

if ~isfield(params, 'concat_fields') || isempty(params.concat_fields)
    params.concat_fields = {'efuns'};
end

if ~isfield(params, 'concat_dim') || isempty(params.concat_dim)
    params.concat_dim = 1;
end

if ~isfield(params, 'required_equal_fields') || isempty(params.required_equal_fields)
    params.required_equal_fields = {'evalues', 'kpm_modes', 'N_dict', 'residual_form'};
end

if ~isfield(params, 'allow_missing_chunks') || isempty(params.allow_missing_chunks)
    params.allow_missing_chunks = false;
end

if ~isfield(params, 'verbose') || isempty(params.verbose)
    params.verbose = true;
end

if ~isfield(params, 'progress_every') || isempty(params.progress_every)
    params.progress_every = 50;
end

params.concat_fields = normalize_cellstr(params.concat_fields, 'params.concat_fields');
params.required_equal_fields = normalize_cellstr(params.required_equal_fields, ...
    'params.required_equal_fields');

if ~isscalar(params.concat_dim) || params.concat_dim < 1 || params.concat_dim ~= floor(params.concat_dim)
    error('params.concat_dim must be a positive integer.');
end

if ~isscalar(params.progress_every) || params.progress_every < 1 || ...
        params.progress_every ~= floor(params.progress_every)
    error('params.progress_every must be a positive integer.');
end
end


function files = collect_chunk_files(data_dir, filename_pattern)
L = dir(fullfile(data_dir, filename_pattern));

if isempty(L)
    files = struct('name', {}, 'fullpath', {}, 'chunk_id', {}, 'prefix', {});
    return;
end

files = repmat(struct('name', '', 'fullpath', '', 'chunk_id', [], 'prefix', ''), numel(L), 1);
n_keep = 0;
ignored_names = cell(0, 1);

for i = 1:numel(L)
    tokens = regexp(L(i).name, '^(.*)_outputs_(\d+)\.mat$', 'tokens', 'once');
    if isempty(tokens)
        ignored_names{end+1, 1} = L(i).name; %#ok<AGROW>
        continue;
    end

    n_keep = n_keep + 1;
    files(n_keep).name = L(i).name;
    files(n_keep).fullpath = fullfile(L(i).folder, L(i).name);
    files(n_keep).chunk_id = str2double(tokens{2});
    files(n_keep).prefix = tokens{1};
end

files = files(1:n_keep);

if isempty(files)
    if isempty(ignored_names)
        error('No files matching %s were found in %s.', filename_pattern, data_dir);
    end

    error(['No files in %s matched the required *_outputs_<chunk>.mat naming pattern. ', ...
        'Ignored examples: %s'], data_dir, strjoin(ignored_names(1:min(end, 3)), ', '));
end

[~, order] = sort([files.chunk_id]);
files = files(order);

prefixes = {files.prefix};
if numel(unique(prefixes)) > 1
    error('Multiple EDMD output prefixes were found in %s. Narrow filename_pattern first.', data_dir);
end
end


function data = load_variable(file_path, variable_name)
data = load(file_path, variable_name);

if ~isfield(data, variable_name)
    error('File %s does not contain variable %s.', file_path, variable_name);
end
end


function check_chunk_sequence(chunk_ids, allow_missing_chunks)
expected = (chunk_ids(1):chunk_ids(end))';
missing = setdiff(expected, chunk_ids);

if isempty(missing)
    return;
end

msg = summarize_integer_list(missing, 12);

if allow_missing_chunks
    warning('Missing chunk IDs detected: %s', msg);
else
    error('Missing chunk IDs detected: %s', msg);
end
end


function fields = validate_field_list(field_list, all_fields, list_name)
fields = normalize_cellstr(field_list, list_name);

missing = setdiff(fields, all_fields);
if ~isempty(missing)
    error('%s contains unknown fields: %s', list_name, strjoin(missing, ', '));
end
end


function out = normalize_cellstr(value, value_name)
if ischar(value) || isstring(value)
    out = cellstr(value);
elseif iscell(value)
    out = value(:)';
else
    error('%s must be a char, string, or cell array of text.', value_name);
end

for i = 1:numel(out)
    if ~(ischar(out{i}) || isstring(out{i}))
        error('%s must contain only text values.', value_name);
    end
    out{i} = char(out{i});
end
end


function check_equal_fields(ref_struct, current_struct, field_names, file_name)
for i = 1:numel(field_names)
    field_name = field_names{i};
    ref_value = read_required_field(ref_struct, field_name, file_name);
    current_value = read_required_field(current_struct, field_name, file_name);

    if ~isequaln(ref_value, current_value)
        error('Field %s is not identical in %s.', field_name, file_name);
    end
end
end


function value = read_required_field(S, field_name, file_name)
if ~isfield(S, field_name)
    error('Field %s is missing from %s.', field_name, file_name);
end

value = S.(field_name);
end


function validate_concat_value(value, field_name)
if isstruct(value) || iscell(value)
    error('Concat field %s must be a plain MATLAB array, not a struct or cell.', field_name);
end
end


function L = get_concat_length(value, concat_dim)
sz = size(value);

if numel(sz) < concat_dim
    sz(end+1:concat_dim) = 1;
end

L = sz(concat_dim);
end


function out_size = build_final_size(ref_size, concat_dim, total_concat_length)
out_size = ref_size;

if numel(out_size) < concat_dim
    out_size(end+1:concat_dim) = 1;
end

out_size(concat_dim) = total_concat_length;
end


function validate_concat_shape(ref_size, current_size, concat_dim, field_name, file_name)
n_dims = max([numel(ref_size), numel(current_size), concat_dim]);
ref_size(end+1:n_dims) = 1;
current_size(end+1:n_dims) = 1;

ref_other = ref_size;
current_other = current_size;
ref_other(concat_dim) = [];
current_other(concat_dim) = [];

if ~isequal(ref_other, current_other)
    error(['Field %s in %s is not compatible for concatenation along dimension %d. ', ...
        'Reference size = [%s], current size = [%s].'], ...
        field_name, file_name, concat_dim, num2str(ref_size), num2str(current_size));
end
end


function chunk_length = validate_chunk_lengths(lengths, file_name)
chunk_length = lengths(1);

if all(lengths == chunk_length)
    return;
end

error('Concatenated fields in %s do not share the same time length.', file_name);
end


function A = assign_concat_block(A, block, concat_dim, start_idx)
idx = repmat({':'}, 1, max(ndims(A), concat_dim));
stop_idx = start_idx + size_with_padding(block, concat_dim) - 1;
idx{concat_dim} = start_idx:stop_idx;
A(idx{:}) = block;
end


function n = size_with_padding(value, dim_idx)
sz = size(value);
if numel(sz) < dim_idx
    sz(end+1:dim_idx) = 1;
end
n = sz(dim_idx);
end


function tf = should_print_progress(k, n_total, progress_every)
tf = (k == 1) || (k == n_total) || (mod(k, progress_every) == 0);
end


function msg = summarize_integer_list(values, max_items)
values = values(:)';
show_values = values(1:min(numel(values), max_items));
msg = strjoin(cellstr(string(show_values)), ', ');

if numel(values) > max_items
    msg = sprintf('%s, ...', msg);
end
end


function nbytes = class_nbytes(class_name)
switch class_name
    case {'double', 'int64', 'uint64'}
        nbytes = 8;
    case {'single', 'int32', 'uint32'}
        nbytes = 4;
    case {'int16', 'uint16', 'char'}
        nbytes = 2;
    case {'int8', 'uint8', 'logical'}
        nbytes = 1;
    otherwise
        error('Unsupported class %s for memory estimation.', class_name);
end
end
