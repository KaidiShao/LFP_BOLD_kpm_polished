% Resume F12 SOP after P5 reduction/P7/P9 by filling P5 density and rerunning P8/P10.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
cd(repo_root);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

target_dataset_stems = {'f12m02', 'f12m03', 'f12m05'}; %#ok<NASGU>
cfg_names = {'F12m02', 'F12m03', 'F12m05'}; %#ok<NASGU>
roi_modes = {'roi_mean'}; %#ok<NASGU>

fprintf('Resume F12 standardized-csplit density + P8/P10\n');
fprintf('  datasets: %s\n', strjoin(target_dataset_stems, ', '));

%% Fill P5 adaptive RMS-envelope density only.
run_reduction = false; %#ok<NASGU>
run_density = true; %#ok<NASGU>
run_peak_stats = false; %#ok<NASGU>
skip_existing_density = false; %#ok<NASGU>
method_filter = {'svd', 'nmf', 'mds', 'umap'}; %#ok<NASGU>
component_count_sweep = 3:8; %#ok<NASGU>
time_component_count_sweep = 3:8; %#ok<NASGU>
spectrum_component_count_sweep = 3:8; %#ok<NASGU>
autodl_root_override = 'E:\DataPons_processed\derived_autodl_results_standardize'; %#ok<NASGU>
allow_summary_only_outputs = false; %#ok<NASGU>

run(fullfile(repo_root, 'scripts', 'script_run_standardized_csplit_p5_for_p12.m'));

%% Rerun P8/P10 xcorr so raw/dimred density sources are included.
main_observable_modes = roi_modes; %#ok<NASGU>
main_feature_names = {'efun_real', 'deconv_real'}; %#ok<NASGU>
main_method_tags = {'svd_k05', 'svd_k06', 'svd_k07', 'svd_k08', ...
    'nmf_k05', 'nmf_k06', 'nmf_k07', 'nmf_k08', ...
    'mds_k05', 'mds_k06', 'mds_k07', 'mds_k08', ...
    'umap_k05', 'umap_k06', 'umap_k07', 'umap_k08'}; %#ok<NASGU>
current_best_roots = {'E:\DataPons_processed\derived_autodl_results_bold_full_export'}; %#ok<NASGU>
force_recompute_standardized_csplit = true; %#ok<NASGU>
residual_forms = {}; %#ok<NASGU>
run_name_filter = {}; %#ok<NASGU>
run_name_contains = {}; %#ok<NASGU>
max_runs = []; %#ok<NASGU>

fprintf('\n[P8 roi_mean standardized csplit density refresh]\n');
run(fullfile(repo_root, 'scripts', 'script_run_standardized_csplit_p8_for_p12.m'));

fprintf('\n[P10 roi_mean standardized csplit density refresh]\n');
run(fullfile(repo_root, 'scripts', 'script_run_standardized_csplit_p10_for_p12.m'));

fprintf('\nFinished F12 density + P8/P10 refresh.\n');
