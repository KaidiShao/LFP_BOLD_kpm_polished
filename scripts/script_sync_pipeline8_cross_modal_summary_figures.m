% Canonical batch sync for pipeline 8 xcorr summary figures.
%
% Usage from MATLAB:
%   run('scripts/script_sync_pipeline8_cross_modal_summary_figures.m');
%
% Optional controls, set before run(...) if needed:
%   dataset_stems = {'e10gb1', 'e10gh1', 'e10fV1'};
%   observable_modes = {'global_svd100', 'svd', 'roi_mean'};
%   residual_forms = {'projected_kv', 'projected_vlambda'};
%   run_name_filter = {};
%   run_name_contains = {};
%   processed_root = '';
%   max_runs = [];

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('processed_root', 'var') || isempty(processed_root)
    processed_root = io_project.get_project_processed_root();
end
if ~exist('dataset_stems', 'var') || isempty(dataset_stems)
    dataset_stems = {};
end
if ~exist('observable_modes', 'var') || isempty(observable_modes)
    observable_modes = {};
end
if ~exist('residual_forms', 'var') || isempty(residual_forms)
    residual_forms = {};
end
if ~exist('run_name_filter', 'var') || isempty(run_name_filter)
    run_name_filter = {};
end
if ~exist('run_name_contains', 'var') || isempty(run_name_contains)
    run_name_contains = {};
end
if ~exist('max_runs', 'var') || isempty(max_runs)
    max_runs = [];
end

params = build_bold_cross_modal_coupling_params();
params.processed_root = processed_root;
params.dataset_stems = cellstr(string(dataset_stems(:)).');
params.observable_modes = cellstr(string(observable_modes(:)).');
params.residual_forms = cellstr(string(residual_forms(:)).');
params.run_name_filter = cellstr(string(run_name_filter(:)).');
params.run_name_contains = cellstr(string(run_name_contains(:)).');
params.max_runs = max_runs;

candidates = discover_completed_bold_cross_modal_coupling_runs(params);
if isempty(candidates)
    error('No pipeline 8 candidates were found for summary sync.');
end
if ~isempty(max_runs)
    candidates = candidates(1:min(numel(candidates), max_runs));
end

sync_rows = repmat(local_empty_row(), numel(candidates), 1);
fprintf('Syncing pipeline 8 summary figures for %d run(s).\n', numel(candidates));

for i_run = 1:numel(candidates)
    candidate = candidates(i_run);
    fprintf('\n[%d/%d] %s | %s | %s\n', ...
        i_run, numel(candidates), candidate.dataset_stem, ...
        candidate.observable_mode, candidate.residual_form);
    fprintf('Run name:\n  %s\n', candidate.run_name);

    sync_result = sync_one_pipeline8_cross_modal_summary_figures(candidate, candidate, processed_root);
    sync_rows(i_run) = local_make_row(sync_result);

    fprintf('Summary sync status: %s\n', sync_result.status);
    if ~isempty(sync_result.summary_dir)
        fprintf('Summary dir:\n  %s\n', sync_result.summary_dir);
    end
    fprintf('Copied figure count: %d\n', sync_result.n_copied_files);
end

sync_table = struct2table(sync_rows);
fprintf('\nPipeline 8 summary sync finished.\n');


function row = local_empty_row()
row = struct('dataset_stem', '', 'run_name', '', 'observable_mode', '', ...
    'residual_form', '', 'scope', '', 'status', '', ...
    'summary_dir', '', 'n_copied_files', 0);
end


function row = local_make_row(sync_result)
row = local_empty_row();
row.dataset_stem = sync_result.dataset_stem;
row.run_name = sync_result.run_name;
row.observable_mode = sync_result.observable_mode;
row.residual_form = sync_result.residual_form;
row.scope = sync_result.scope;
row.status = sync_result.status;
row.summary_dir = sync_result.summary_dir;
row.n_copied_files = sync_result.n_copied_files;
end
