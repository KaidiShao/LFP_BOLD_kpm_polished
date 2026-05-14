% Convenience wrapper: run pipeline 7 for one completed BOLD ResKoopNet run.
%
% Canonical interactive entry:
%   scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m
%
% This wrapper is still useful when you already know the exact completed
% AutoDL output folder and want to postprocess just that one run.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

warning(['script_run_one_bold_reskoopnet_postprocess is now a convenience ', ...
    'wrapper. Use script_run_one_cfg_bold_reskoopnet_postprocessing.m ', ...
    'for the canonical step-by-step pipeline 7 workflow.']);

%% -------------------- user settings --------------------
result_dir = ['E:\autodl_results_local\bold_wsl\e10gb1\mlp\outputs\', ...
    'mlp_obs_bold_wsl_20260423_e10gb1_projected_kv_HP_svd100'];

dataset_stem = 'e10gb1';
dataset_id = '';
observable_mode = 'HP_svd100';

skip_existing = false;
make_main_plot = true;
make_deconv_plot = true;
make_timescale_plot = true;

%% -------------------- run pipeline 7 --------------------
if exist(result_dir, 'dir') ~= 7
    error('result_dir does not exist: %s', result_dir);
end

[~, run_name] = fileparts(result_dir);
if isempty(dataset_id)
    dataset_id = io_project.get_dataset_id_from_stem(dataset_stem);
end

run_info = struct();
run_info.dataset_stem = dataset_stem;
run_info.dataset_id = dataset_id;
run_info.run_name = run_name;
run_info.output_dir = result_dir;
run_info.observable_mode = observable_mode;

params = build_bold_reskoopnet_postprocessing_params();
params.processed_root = io_project.get_project_processed_root();
params.skip_existing = logical(skip_existing);
params.make_main_plot = logical(make_main_plot);
params.make_deconv_plot = logical(make_deconv_plot);
params.make_timescale_plot = logical(make_timescale_plot);
params.headless = false;
params.close_figures = false;

result = run_one_bold_reskoopnet_postprocessing_core(run_info, params);

fprintf('\nFinished pipeline 7 for one BOLD ResKoopNet run.\n');
fprintf('Post MAT:\n  %s\n', result.post_file);
fprintf('Main figure:\n  %s\n', result.main_png);
fprintf('Deconv figure:\n  %s\n', result.deconv_png);
fprintf('Timescale figure:\n  %s\n', result.timescale_png);
if isfield(result, 'intrinsic_activation_dir') && ~isempty(result.intrinsic_activation_dir)
    fprintf('Intrinsic activation maps:\n  %s\n', result.intrinsic_activation_dir);
end
if isfield(result, 'intrinsic_roi_summary_png') && ~isempty(result.intrinsic_roi_summary_png)
    fprintf('Intrinsic ROI summary:\n  %s\n', result.intrinsic_roi_summary_png);
end
