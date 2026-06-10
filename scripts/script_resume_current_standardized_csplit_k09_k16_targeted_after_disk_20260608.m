% Targeted resume after the 2026-06-08 disk-full interruption.
%
% Completed before this script:
%   - P5 k09:k16 reductions for all current datasets.
%   - P5 adaptive density for all except k13m20/k13m21 k09:k16.
%   - P8 for most datasets, but k13m20/k13m21 lacked dimred density sources.
%   - P10 through f12m03; f12m05 failed while writing deconv_real svd_k08.
%
% This script intentionally does only the missing work:
%   1. P5 adaptive dimred density k09:k16 for k13m20/k13m21.
%   2. P8 roi_mean rerun for k13m20/k13m21.
%   3. P10 roi_mean rerun from f12m05 onward.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

fprintf('Targeted standardized complex-split k03:k16 resume after disk-full interruption\n');

%% 1. Missing P5 adaptive density for k13m20/k13m21 k09:k16.
target_dataset_stems = {'k13m20', 'k13m21'}; %#ok<NASGU>
run_reduction = false; %#ok<NASGU>
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

fprintf('\n[1/3] Backfilling P5 adaptive density for k13m20/k13m21 k09:k16.\n');
run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p5_for_p12.m'));

%% 2. P8 rerun for k13m20/k13m21 so dimred sources are included.
clear target_dataset_stems run_reduction run_density run_peak_stats
clear skip_existing_reduction skip_existing_density time_component_count_sweep
clear spectrum_component_count_sweep method_filter allow_summary_only_outputs
clear autodl_root_override

cfg_names = {'K13m20', 'K13m21'}; %#ok<NASGU>
current_p7_run_names = { ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_k13m20_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_k13m21_projected_vlambda_roi_mean'}; %#ok<NASGU>
component_count_sweep = 3:16; %#ok<NASGU>
main_observable_modes = {'roi_mean'}; %#ok<NASGU>
force_recompute_standardized_csplit = true; %#ok<NASGU>

fprintf('\n[2/3] Rerunning P8 for k13m20/k13m21 with k03:k16 dimred density sources.\n');
run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p8_for_p12.m'));

%% 3. P10 rerun from f12m05 onward.
clear cfg_names current_p7_run_names run_name_filter run_name_contains
clear residual_forms max_runs

cfg_names = {'F12m05', 'K13m17', 'K13m18', 'K13m19', ...
    'K13m20', 'K13m21', 'K13m23'}; %#ok<NASGU>
current_p7_run_names = { ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_f12m05_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_20260521_allbold_k13m17_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_k13m18_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_k13m19_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_k13m20_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allbold_k13m21_projected_vlambda_roi_mean', ...
    'mlp_obs_bold_batch4_roi_mean_stdT_l1e4_r1e3_b2000_i2_20260521_allbold_k13m23_projected_vlambda_roi_mean'}; %#ok<NASGU>
component_count_sweep = 3:16; %#ok<NASGU>
main_observable_modes = {'roi_mean'}; %#ok<NASGU>
main_feature_names = {'efun_real', 'deconv_real'}; %#ok<NASGU>
force_recompute_standardized_csplit = true; %#ok<NASGU>

fprintf('\n[3/3] Rerunning P10 from f12m05 onward with BLP k03:k16 density sources.\n');
run(fullfile(repo_root, 'scripts', ...
    'script_run_standardized_csplit_p10_for_p12.m'));

fprintf('\nFinished targeted standardized complex-split k03:k16 density/P8/P10 resume after disk-full interruption.\n');
