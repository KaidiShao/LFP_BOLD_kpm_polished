% Worker branch for the current standardized complex-split k09:k16 extension.
%
% This worker only fills P5 reduction + adaptive RMS-envelope density for the
% K13 standardized complex-split datasets.  The primary k09:k16 process is
% still responsible for the final P8/P10 refresh and Python summary.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

target_dataset_stems = {'k13m17', 'k13m18', 'k13m19', ...
    'k13m20', 'k13m21', 'k13m23'}; %#ok<NASGU>

run_reduction = true; %#ok<NASGU>
run_density = true; %#ok<NASGU>
run_peak_stats = false; %#ok<NASGU>
skip_existing_reduction = true; %#ok<NASGU>
skip_existing_density = true; %#ok<NASGU>
component_count_sweep = 9:16; %#ok<NASGU>
time_component_count_sweep = 9:16; %#ok<NASGU>
spectrum_component_count_sweep = 9:16; %#ok<NASGU>
method_filter = {'svd', 'nmf', 'mds', 'umap'}; %#ok<NASGU>
allow_summary_only_outputs = false; %#ok<NASGU>
autodl_root_override = 'E:\DataPons_processed\derived_autodl_results_standardize'; %#ok<NASGU>

fprintf('K13 worker for standardized complex-split P5 k09:k16 extension\n');
fprintf('Datasets: %s\n', strjoin(target_dataset_stems, ', '));

run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p5_for_p12.m'));

fprintf('\nFinished K13 worker for standardized complex-split P5 k09:k16 extension.\n');
