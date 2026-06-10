% Targeted backfill for K13m23 standardized complex-split P8/P10 xcorr.
%
% Current audit found K13m23 missing only the global_svd100 and gsvd100_ds
% BOLD observable branches.  Keep numeric-first outputs only.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg_names = {'K13m23'};
main_observable_modes = {'global_svd100', 'gsvd100_ds'};

run(fullfile(repo_root, 'scripts', 'script_run_standardized_csplit_p8_for_p12.m'));
run(fullfile(repo_root, 'scripts', 'script_run_standardized_csplit_p10_for_p12.m'));
