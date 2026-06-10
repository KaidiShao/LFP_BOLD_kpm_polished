function candidates = discover_completed_bold_eigenfunction_reduction_runs(params)
%DISCOVER_COMPLETED_BOLD_EIGENFUNCTION_REDUCTION_RUNS Find P7 BOLD_POST inputs for P9.

if nargin < 1 || isempty(params)
    params = build_bold_eigenfunction_reduction_params();
else
    params = local_apply_defaults(params);
end

discover_params = struct();
copy_fields = {'processed_root', 'dataset_stems', 'exclude_dataset_stems', ...
    'observable_modes', 'residual_forms', 'run_name_filter', 'run_name_contains', ...
    'current_best_p7_only', 'current_p7_run_names', 'current_p7_autodl_roots', ...
    'current_p7_allow_summary_only_outputs'};
for i = 1:numel(copy_fields)
    name = copy_fields{i};
    if isfield(params, name)
        discover_params.(name) = params.(name);
    end
end

candidates = discover_completed_bold_cross_modal_coupling_runs(discover_params);
if isempty(candidates)
    return;
end

for i = 1:numel(candidates)
    candidates(i).p9_root = fullfile( ...
        io_project.get_pipeline_stage_dir(params.processed_root, ...
        candidates(i).dataset_stem, 9, 'bold_eigenfunction_reduction'), ...
        candidates(i).run_tag);
end

if ~isempty(params.max_runs)
    n_keep = min(numel(candidates), double(params.max_runs));
    candidates = candidates(1:n_keep);
end
end


function params = local_apply_defaults(params)
defaults = build_bold_eigenfunction_reduction_params();
params = local_merge_defaults(defaults, params);
params.dataset_stems = cellstr(string(params.dataset_stems(:)).');
params.exclude_dataset_stems = cellstr(string(params.exclude_dataset_stems(:)).');
params.observable_modes = cellstr(string(params.observable_modes(:)).');
params.residual_forms = cellstr(string(params.residual_forms(:)).');
params.run_name_filter = cellstr(string(params.run_name_filter(:)).');
params.run_name_contains = cellstr(string(params.run_name_contains(:)).');
end


function out = local_merge_defaults(defaults, overrides)
out = defaults;
names = fieldnames(overrides);
for i = 1:numel(names)
    name = names{i};
    value = overrides.(name);
    if isstruct(value) && isfield(defaults, name) && isstruct(defaults.(name)) && ...
            isscalar(value) && isscalar(defaults.(name))
        out.(name) = local_merge_defaults(defaults.(name), value);
    elseif ~isempty(value)
        out.(name) = value;
    else
        out.(name) = defaults.(name);
    end
end
end
