function runs = filter_invalid_blp_mlp_runs(runs, params)
%FILTER_INVALID_BLP_MLP_RUNS Remove BLP MLP runs marked invalid/stale.
%
% The registry lets current P5/P6-style discovery downgrade a bad P4 source
% without deleting files on disk. Set params.include_invalid_runs = true only
% for explicit forensic/history runs.

if nargin < 2
    params = struct();
end
if isempty(runs)
    return;
end
if isfield(params, 'include_invalid_runs') && ...
        ~isempty(params.include_invalid_runs) && logical(params.include_invalid_runs)
    return;
end

registry_file = local_registry_file(params);
if exist(registry_file, 'file') ~= 2
    return;
end

try
    T = readtable(registry_file, 'TextType', 'string');
catch ME
    warning('Could not read invalid BLP MLP run registry %s: %s', ...
        registry_file, ME.message);
    return;
end
if isempty(T)
    return;
end

status = lower(strtrim(local_table_col(T, 'status', "invalid")));
invalid_statuses = local_invalid_statuses(params);
active_rows = ismember(status, invalid_statuses);
if ~any(active_rows)
    return;
end

keep = true(numel(runs), 1);
reason_col = local_table_col(T, 'reason', "");
for i = 1:numel(runs)
    for j = find(active_rows).'
        if local_row_matches_run(T(j, :), runs(i))
            keep(i) = false;
            fprintf(2, ['Skipping invalid/stale BLP MLP run: %s | %s | %s | %s ' ...
                '(%s: %s)\n'], runs(i).dataset, runs(i).observable_mode, ...
                runs(i).residual_form, runs(i).run_name, status(j), reason_col(j));
            break;
        end
    end
end
runs = runs(keep);
end


function registry_file = local_registry_file(params)
if isfield(params, 'invalid_run_registry_file') && ...
        ~isempty(params.invalid_run_registry_file)
    registry_file = char(string(params.invalid_run_registry_file));
    return;
end

this_file = mfilename('fullpath');
repo_root = fileparts(fileparts(fileparts(this_file)));
registry_file = fullfile(repo_root, 'configs', 'current_invalid_blp_mlp_runs.csv');
end


function values = local_table_col(T, name, default_value)
names = string(T.Properties.VariableNames);
idx = find(strcmpi(names, string(name)), 1);
if isempty(idx)
    values = repmat(string(default_value), height(T), 1);
else
    values = string(T.(char(names(idx))));
end
values = values(:);
end


function statuses = local_invalid_statuses(params)
if isfield(params, 'invalid_run_statuses') && ...
        ~isempty(params.invalid_run_statuses)
    statuses = lower(strtrim(string(params.invalid_run_statuses(:))));
else
    statuses = ["invalid"; "stale"; "legacy"];
end
end


function tf = local_row_matches_run(row, run_info)
dataset = local_table_col(row, 'dataset', "");
observable_mode = local_table_col(row, 'observable_mode', "");
residual_form = local_table_col(row, 'residual_form', "");
run_name = local_table_col(row, 'run_name', "");
output_dir = local_table_col(row, 'output_dir', "");

tf = local_field_matches(run_info.dataset, dataset) && ...
    local_field_matches(run_info.observable_mode, observable_mode) && ...
    local_field_matches(run_info.residual_form, residual_form) && ...
    local_field_matches(run_info.run_name, run_name) && ...
    local_path_matches(run_info.output_dir, output_dir);
end


function tf = local_field_matches(actual, expected)
expected = strtrim(string(expected));
if expected == "" || expected == "*"
    tf = true;
else
    tf = strcmpi(strtrim(string(actual)), expected);
end
end


function tf = local_path_matches(actual, expected)
expected = strtrim(string(expected));
if expected == "" || expected == "*"
    tf = true;
else
    tf = strcmpi(local_normalize_path(actual), local_normalize_path(expected));
end
end


function out = local_normalize_path(value)
out = lower(strtrim(string(value)));
out = replace(out, '/', '\');
out = regexprep(out, '[\\]+$', '');
end
