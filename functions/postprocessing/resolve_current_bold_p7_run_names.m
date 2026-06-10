function run_names = resolve_current_bold_p7_run_names(params)
%RESOLVE_CURRENT_BOLD_P7_RUN_NAMES Return P7 run names selected from BOLD P4.
%
% This helper keeps downstream P8/P10 tied to the current P7 mainline:
% one best-validation BOLD checkpoint per dataset x observable x residual.

if nargin < 1 || isempty(params)
    params = struct();
end

if isfield(params, 'current_p7_run_names') && ~isempty(params.current_p7_run_names)
    run_names = local_unique_cellstr(params.current_p7_run_names);
    return;
end

p7 = build_bold_reskoopnet_postprocessing_params();

copy_fields = {'processed_root', 'dataset_stems', 'exclude_dataset_stems', ...
    'observable_modes', 'residual_forms', 'run_name_filter', ...
    'run_name_contains'};
for i = 1:numel(copy_fields)
    name = copy_fields{i};
    if isfield(params, name) && ~isempty(params.(name))
        p7.(name) = params.(name);
    end
end

if isfield(params, 'current_p7_autodl_roots') && ~isempty(params.current_p7_autodl_roots)
    p7.autodl_roots = params.current_p7_autodl_roots;
elseif isfield(params, 'autodl_roots') && ~isempty(params.autodl_roots)
    p7.autodl_roots = params.autodl_roots;
end
if isfield(params, 'current_p7_allow_summary_only_outputs') && ...
        ~isempty(params.current_p7_allow_summary_only_outputs)
    p7.allow_summary_only_outputs = logical(params.current_p7_allow_summary_only_outputs);
elseif isfield(params, 'allow_summary_only_outputs') && ~isempty(params.allow_summary_only_outputs)
    p7.allow_summary_only_outputs = logical(params.allow_summary_only_outputs);
end

p7.dedupe_by_condition = true;
p7.run_selection_mode = 'best_val_metric';
p7.include_smoke = false;
p7.require_summary_file = true;

runs = discover_completed_bold_reskoopnet_runs(p7);
if isempty(runs)
    run_names = {};
    return;
end

run_names = local_unique_cellstr({runs.run_name});
end


function values = local_unique_cellstr(values)
values = cellstr(string(values(:)).');
values = values(~cellfun(@isempty, values));
[~, idx] = unique(lower(string(values)), 'stable');
values = values(sort(idx));
end
