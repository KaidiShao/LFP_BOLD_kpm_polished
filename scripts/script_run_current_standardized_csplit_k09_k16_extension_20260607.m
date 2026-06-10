% Extend current standardized complex-split BLP dimred grid to k09:k16.
%
% Scope:
%   - Current E10/F12/K13 standardized complex-split BLP datasets.
%   - P5 reduction and adaptive RMS-envelope density for missing k09:k16.
%   - Refresh minimal-SOP roi_mean P8/P10 coupling with BLP density sources
%     expanded to k03:k16.
%
% BOLD-side P9 method tags remain the current minimal-SOP k05:k08 grid; this
% script expands the BLP dimred density/component axis, not the BOLD-side
% reduction axis.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

target_dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'e10gw1', ...
    'f12m01', 'f12m02', 'f12m03', 'f12m05', ...
    'k13m17', 'k13m18', 'k13m19', 'k13m20', 'k13m21', 'k13m23'};
cfg_names = {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', ...
    'F12m01', 'F12m02', 'F12m03', 'F12m05', ...
    'K13m17', 'K13m18', 'K13m19', 'K13m20', 'K13m21', 'K13m23'};

fprintf('Current standardized complex-split k09:k16 extension\n');
fprintf('Datasets: %s\n', strjoin(target_dataset_stems, ', '));

%% P5: only the missing higher-k reduction/density branch.
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

run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p5_for_p12.m'));

%% P8/P10: refresh roi_mean current-SOP coupling with BLP k03:k16 density.
clear run_name_filter run_name_contains residual_forms max_runs current_p7_run_names
component_count_sweep = 3:16; %#ok<NASGU>
main_observable_modes = {'roi_mean'}; %#ok<NASGU>
current_best_roots = {'E:\autodl_results_local\bold_wsl', ...
    'E:\DataPons_processed\derived_autodl_results_bold_full_export'}; %#ok<NASGU>
force_recompute_standardized_csplit = true; %#ok<NASGU>

fprintf('\nRefreshing P8 roi_mean xcorr with BLP dimred k03:k16 density sources.\n');
run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p8_for_p12.m'));

clear run_name_filter run_name_contains residual_forms max_runs current_p7_run_names
component_count_sweep = 3:16; %#ok<NASGU>
main_observable_modes = {'roi_mean'}; %#ok<NASGU>
main_feature_names = {'efun_real', 'deconv_real'}; %#ok<NASGU>
current_best_roots = {'E:\autodl_results_local\bold_wsl', ...
    'E:\DataPons_processed\derived_autodl_results_bold_full_export'}; %#ok<NASGU>
force_recompute_standardized_csplit = true; %#ok<NASGU>

fprintf('\nRefreshing P10 roi_mean xcorr with BLP dimred k03:k16 density sources.\n');
run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p10_for_p12.m'));

fprintf('\nFinished current standardized complex-split k09:k16 extension.\n');
