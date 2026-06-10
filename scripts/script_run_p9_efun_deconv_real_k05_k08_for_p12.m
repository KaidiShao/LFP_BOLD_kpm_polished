% Backfill P9 BOLD efun/deconv_real reductions needed by P12/P10.
%
% Numeric-first: summary plots are off. Existing MAT files are skipped.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    cfg_names = {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', 'F12m01', ...
        'K13m17', 'K13m23'};
end
if ~exist('main_observable_modes', 'var') || isempty(main_observable_modes)
    main_observable_modes = {'HP_svd100', 'global_svd100', 'gsvd100_ds', 'roi_mean'};
end
if ~exist('current_best_roots', 'var') || isempty(current_best_roots)
    current_best_roots = {'E:\autodl_results_local\bold_wsl', ...
        'E:\DataPons_processed\derived_autodl_results_bold_full_export'};
end
if ~exist('main_feature_names', 'var') || isempty(main_feature_names)
    main_feature_names = {'efun_real', 'deconv_real'};
end
if ~exist('main_force_recompute', 'var') || isempty(main_force_recompute)
    main_force_recompute = false;
end

for i_cfg = 1:numel(cfg_names)
    cfg_name = cfg_names{i_cfg}; %#ok<NASGU>
    observable_modes = main_observable_modes; %#ok<NASGU>
    current_best_p7_only = true; %#ok<NASGU>
    current_p7_autodl_roots = current_best_roots; %#ok<NASGU>
    feature_names = main_feature_names; %#ok<NASGU>
    method_filter = {'svd', 'nmf', 'mds', 'umap'}; %#ok<NASGU>
    component_count_sweep = 5:8; %#ok<NASGU>
    force_recompute = main_force_recompute; %#ok<NASGU>
    make_summary_plot = false; %#ok<NASGU>
    continue_on_error = true; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_bold_eigenfunction_reduction.m'));
end
