% Run P10 xcorr for standardized complex-split RMS-envelope BLP density.
%
% Uses P9 efun_real/deconv_real BOLD reductions and compares them against:
%   event density, raw csplit density, dimred csplit density.
%
% Numeric-first: xcorr overview figures, activation maps, and ROI summaries
% are disabled. Existing xcorr MAT/CSV outputs are skipped.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    cfg_names = {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', 'F12m01', ...
        'F12m02', 'F12m03', 'F12m05', ...
        'K13m17', 'K13m18', 'K13m19', 'K13m20', 'K13m21', 'K13m23'};
end
if ~exist('main_observable_modes', 'var') || isempty(main_observable_modes)
    main_observable_modes = {'HP_svd100', 'global_svd100', 'gsvd100_ds', 'roi_mean'};
end
if ~exist('main_feature_names', 'var') || isempty(main_feature_names)
    main_feature_names = {'efun_real', 'deconv_real'};
end
if ~exist('current_best_roots', 'var') || isempty(current_best_roots)
    current_best_roots = {'E:\autodl_results_local\bold_wsl', ...
        'E:\DataPons_processed\derived_autodl_results_bold_full_export'};
end
if ~exist('main_method_tags', 'var') || isempty(main_method_tags)
    main_method_tags = {'svd_k05', 'svd_k06', 'svd_k07', 'svd_k08', ...
        'nmf_k05', 'nmf_k06', 'nmf_k07', 'nmf_k08', ...
        'mds_k05', 'mds_k06', 'mds_k07', 'mds_k08', ...
        'umap_k05', 'umap_k06', 'umap_k07', 'umap_k08'};
end
if ~exist('force_recompute_standardized_csplit', 'var') || isempty(force_recompute_standardized_csplit)
    force_recompute_standardized_csplit = false;
end
if ~exist('component_count_sweep', 'var') || isempty(component_count_sweep)
    component_count_sweep = 3:8;
end
if ~exist('blp_dimred_method_tags', 'var') || isempty(blp_dimred_method_tags)
    blp_dimred_method_tags = local_make_method_tags({'svd', 'nmf', 'mds', 'umap'}, ...
        component_count_sweep);
end

for i_cfg = 1:numel(cfg_names)
    cfg_name = cfg_names{i_cfg}; %#ok<NASGU>
    observable_modes = main_observable_modes; %#ok<NASGU>
    feature_names = main_feature_names; %#ok<NASGU>
    method_tags = main_method_tags; %#ok<NASGU>
    density_source_kinds = {'event_density', ...
        'raw_complex_split_density', 'dimred_complex_split_density'}; %#ok<NASGU>
    blp_dimred_method_tags = cellstr(string(blp_dimred_method_tags(:)).'); %#ok<NASGU>
    blp_density_condition_suffix = 'standardize_rmsenv_adaptive'; %#ok<NASGU>
    xcorr_save_tag = 'dimred_xcorr_csplit_standardize_rmsenv_adaptive'; %#ok<NASGU>
    force_recompute = logical(force_recompute_standardized_csplit); %#ok<NASGU>
    current_best_p7_only = true; %#ok<NASGU>
    current_p7_autodl_roots = current_best_roots; %#ok<NASGU>
    current_p7_allow_summary_only_outputs = true; %#ok<NASGU>
    output_mode = 'separate'; %#ok<NASGU>
    top_n = 50; %#ok<NASGU>
    make_xcorr_figures = false; %#ok<NASGU>
    make_activation_maps = false; %#ok<NASGU>
    make_roi_summaries = false; %#ok<NASGU>
    require_all_density_sources = false; %#ok<NASGU>
    continue_on_error = true; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_bold_dimred_cross_modal_coupling.m'));
end


function tags = local_make_method_tags(methods, counts)
methods = cellstr(string(methods(:)).');
counts = double(counts(:)).';
tags = cell(1, numel(methods) * numel(counts));
idx = 0;
for i_method = 1:numel(methods)
    for k = counts
        idx = idx + 1;
        tags{idx} = sprintf('%s_k%02d', methods{i_method}, k);
    end
end
end
