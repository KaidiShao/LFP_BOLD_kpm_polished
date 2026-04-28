function [EDMD_outputs, concat_info] = load_and_concat_edmd_output_chunks(data_dir, params)
% Load EDMD output chunk files and concatenate selected fields in time.
%
% This loader uses a simple two-pass flow:
%   1. scan all chunks to validate fields and measure lengths
%   2. preallocate the output arrays and copy each chunk into place
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
%                                        'N_dict','residual_form',...
%                                        'observable_tag','observable_mode'})
%            .allow_missing_chunks   (default false)
%            .verbose                (default true)
%
% Outputs
%   EDMD_outputs Concatenated output struct
%   concat_info  Basic metadata for the concatenation

if nargin < 1 || isempty(data_dir)
    error('data_dir must be provided.');
end

if nargin < 2
    params = struct();
end

if isfield(params, 'scan_mode') && ~isempty(params.scan_mode)
    scan_mode = char(lower(string(params.scan_mode)));
    if ~strcmp(scan_mode, 'full')
        warning(['Ignoring params.scan_mode = %s. ', ...
            'This loader now always uses a full two-pass scan.'], scan_mode);
    end
end

params = apply_default_params(params);

% Resolve chunk files
files = collect_chunk_files(data_dir, params.filename_pattern);

if isempty(files)
    error('No files matching %s were found in %s.', params.filename_pattern, data_dir);
end

chunk_ids = [files.chunk_id]';
expected_chunk_ids = (chunk_ids(1):chunk_ids(end))';
missing_chunk_ids = setdiff(expected_chunk_ids, chunk_ids);
if ~isempty(missing_chunk_ids)
    shown_missing = missing_chunk_ids(1:min(numel(missing_chunk_ids), 12));
    missing_msg = strjoin(cellstr(string(shown_missing(:).')), ', ');
    if numel(missing_chunk_ids) > 12
        missing_msg = sprintf('%s, ...', missing_msg);
    end

    if params.allow_missing_chunks
        warning('Missing chunk IDs detected: %s', missing_msg);
    else
        error('Missing chunk IDs detected: %s', missing_msg);
    end
end

if params.verbose
    fprintf('Found %d EDMD chunk files in %s\n', numel(files), data_dir);
    fprintf('Chunk range: %d -> %d\n', chunk_ids(1), chunk_ids(end));
end

% Load the first chunk and decide which fields to concatenate
first_data = load(files(1).fullpath, params.variable_name);
if ~isfield(first_data, params.variable_name)
    error('File %s does not contain variable %s.', files(1).fullpath, params.variable_name);
end
ref_outputs = first_data.(params.variable_name);
EDMD_outputs = ref_outputs;

all_fields = fieldnames(ref_outputs);
concat_fields = params.concat_fields;
missing_concat_fields = setdiff(concat_fields, all_fields);
if ~isempty(missing_concat_fields)
    error('params.concat_fields contains unknown fields: %s', ...
        strjoin(missing_concat_fields, ', '));
end
copied_fields = setdiff(all_fields, concat_fields, 'stable');
equal_fields = intersect(params.required_equal_fields, copied_fields, 'stable');

n_chunks = numel(files);
n_concat_fields = numel(concat_fields);

chunk_lengths = zeros(n_chunks, 1);
concat_lengths_by_field = zeros(n_chunks, n_concat_fields);
ref_values = cell(n_concat_fields, 1);
ref_sizes = cell(n_concat_fields, 1);
final_sizes = cell(n_concat_fields, 1);

% First pass: validate chunk contents and measure concatenated lengths
for i = 1:n_concat_fields
    field_name = concat_fields{i};
    ref_values{i} = ref_outputs.(field_name);
    ref_sizes{i} = size(ref_values{i});
    if isstruct(ref_values{i}) || iscell(ref_values{i})
        error('Concat field %s must be a plain MATLAB array, not a struct or cell.', field_name);
    end

    sz = ref_sizes{i};
    if numel(sz) < params.concat_dim
        sz(end+1:params.concat_dim) = 1;
    end
    concat_lengths_by_field(1, i) = sz(params.concat_dim);
end

chunk_lengths(1) = concat_lengths_by_field(1, 1);
if any(concat_lengths_by_field(1, :) ~= chunk_lengths(1))
    error('Concatenated fields in %s do not share the same time length.', files(1).name);
end

for k = 2:n_chunks
    current_data = load(files(k).fullpath, params.variable_name);
    if ~isfield(current_data, params.variable_name)
        error('File %s does not contain variable %s.', files(k).fullpath, params.variable_name);
    end
    current_outputs = current_data.(params.variable_name);
    concat_lengths_by_field(k, :) = local_measure_chunk_lengths( ...
        ref_outputs, current_outputs, equal_fields, concat_fields, ...
        ref_sizes, params.concat_dim, files(k).name);
    chunk_lengths(k) = concat_lengths_by_field(k, 1);
end

total_length = sum(chunk_lengths);

% Second pass: preallocate output arrays and copy chunks into place
for i = 1:n_concat_fields
    field_name = concat_fields{i};
    final_sizes{i} = ref_sizes{i};
    if numel(final_sizes{i}) < params.concat_dim
        final_sizes{i}(end+1:params.concat_dim) = 1;
    end
    final_sizes{i}(params.concat_dim) = sum(concat_lengths_by_field(:, i));
    EDMD_outputs.(field_name) = zeros(final_sizes{i}, 'like', ref_values{i});
end

chunk_start_idx = zeros(n_chunks, 1);
chunk_end_idx = zeros(n_chunks, 1);
cursor = 1;

for k = 1:n_chunks
    if k == 1
        current_outputs = ref_outputs;
    else
        current_data = load(files(k).fullpath, params.variable_name);
        if ~isfield(current_data, params.variable_name)
            error('File %s does not contain variable %s.', files(k).fullpath, params.variable_name);
        end
        current_outputs = current_data.(params.variable_name);
    end

    chunk_start_idx(k) = cursor;
    chunk_end_idx(k) = cursor + chunk_lengths(k) - 1;

    EDMD_outputs = local_assign_concat_blocks( ...
        EDMD_outputs, current_outputs, concat_fields, params.concat_dim, cursor);

    cursor = chunk_end_idx(k) + 1;
end

% Assemble a thin concat summary that downstream code can still use
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

if params.verbose
    fprintf('Concatenated %d chunks into %d samples.\n', n_chunks, total_length);
    for i = 1:n_concat_fields
        field_name = concat_fields{i};
        fprintf('  %s size = [%s]\n', field_name, ...
            num2str(size(EDMD_outputs.(field_name))));
    end
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
    params.required_equal_fields = { ...
        'evalues', 'kpm_modes', 'N_dict', 'residual_form', ...
        'observable_tag', 'observable_mode', ...
        'dt', 'dx', 'sampling_period', 'sample_period', 'fs', 'sampling_frequency', ...
        'session_dx', 'session_fs', 'session_ids', 'session_lengths', ...
        'session_start_idx', 'session_end_idx', 'border_idx'};
end

if ~isfield(params, 'allow_missing_chunks') || isempty(params.allow_missing_chunks)
    params.allow_missing_chunks = false;
end

if ~isfield(params, 'verbose') || isempty(params.verbose)
    params.verbose = true;
end

if ischar(params.concat_fields) || isstring(params.concat_fields)
    params.concat_fields = cellstr(params.concat_fields);
elseif iscell(params.concat_fields)
    params.concat_fields = params.concat_fields(:)';
else
    error('params.concat_fields must be a char, string, or cell array of text.');
end
for i = 1:numel(params.concat_fields)
    if ~(ischar(params.concat_fields{i}) || isstring(params.concat_fields{i}))
        error('params.concat_fields must contain only text values.');
    end
    params.concat_fields{i} = char(params.concat_fields{i});
end

if ischar(params.required_equal_fields) || isstring(params.required_equal_fields)
    params.required_equal_fields = cellstr(params.required_equal_fields);
elseif iscell(params.required_equal_fields)
    params.required_equal_fields = params.required_equal_fields(:)';
else
    error('params.required_equal_fields must be a char, string, or cell array of text.');
end
for i = 1:numel(params.required_equal_fields)
    if ~(ischar(params.required_equal_fields{i}) || isstring(params.required_equal_fields{i}))
        error('params.required_equal_fields must contain only text values.');
    end
    params.required_equal_fields{i} = char(params.required_equal_fields{i});
end

params.allow_missing_chunks = logical(params.allow_missing_chunks);
params.verbose = logical(params.verbose);

if ~isscalar(params.concat_dim) || params.concat_dim < 1 || ...
        params.concat_dim ~= floor(params.concat_dim)
    error('params.concat_dim must be a positive integer.');
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


function lengths = local_measure_chunk_lengths( ...
    ref_outputs, current_outputs, equal_fields, concat_fields, ref_sizes, concat_dim, file_name)
for i = 1:numel(equal_fields)
    field_name = equal_fields{i};
    if ~isfield(current_outputs, field_name)
        error('Field %s is missing from %s.', field_name, file_name);
    end
    if ~isequaln(ref_outputs.(field_name), current_outputs.(field_name))
        error('Field %s is not identical in %s.', field_name, file_name);
    end
end

lengths = zeros(1, numel(concat_fields));
for i = 1:numel(concat_fields)
    field_name = concat_fields{i};
    if ~isfield(current_outputs, field_name)
        error('Field %s is missing from %s.', field_name, file_name);
    end

    value = current_outputs.(field_name);
    if isstruct(value) || iscell(value)
        error('Concat field %s must be a plain MATLAB array, not a struct or cell.', field_name);
    end

    current_size = size(value);
    n_dims = max([numel(ref_sizes{i}), numel(current_size), concat_dim]);
    ref_size = ref_sizes{i};
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

    lengths(i) = current_size(concat_dim);
end

if any(lengths ~= lengths(1))
    error('Concatenated fields in %s do not share the same time length.', file_name);
end
end


function EDMD_outputs = local_assign_concat_blocks( ...
    EDMD_outputs, current_outputs, concat_fields, concat_dim, start_idx)
for i = 1:numel(concat_fields)
    field_name = concat_fields{i};
    block = current_outputs.(field_name);
    block_size = size(block);
    if numel(block_size) < concat_dim
        block_size(end+1:concat_dim) = 1;
    end

    idx = repmat({':'}, 1, max(ndims(EDMD_outputs.(field_name)), concat_dim));
    idx{concat_dim} = start_idx:(start_idx + block_size(concat_dim) - 1);
    EDMD_outputs.(field_name)(idx{:}) = block;
end
end
