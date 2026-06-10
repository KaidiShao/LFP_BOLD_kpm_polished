% Batch wrapper for P10 xcorr against P5 *_rmsenv_adaptive density sources.
%
% Numeric-first: activation maps, ROI summaries, and browse figures are off.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

cfg_names = {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', 'F12m01'};
main_observable_modes = {'HP_svd100', 'global_svd100', 'gsvd100_ds', 'roi_mean'};
main_feature_names = {'efun_real', 'deconv_real'};
main_method_tags = {'svd_k05', 'svd_k08', ...
    'nmf_k05', 'nmf_k08', ...
    'mds_k05', 'mds_k08', ...
    'umap_k05', 'umap_k08'};

for i_cfg = 1:numel(cfg_names)
    cfg_name = cfg_names{i_cfg}; %#ok<NASGU>
    blp_density_condition_suffix = 'rmsenv_adaptive'; %#ok<NASGU>
    xcorr_save_tag = 'dimred_xcorr_rmsenv_adaptive'; %#ok<NASGU>
    observable_modes = main_observable_modes; %#ok<NASGU>
    feature_names = main_feature_names; %#ok<NASGU>
    method_tags = main_method_tags; %#ok<NASGU>
    force_recompute = false; %#ok<NASGU>
    current_best_p7_only = true; %#ok<NASGU>
    output_mode = 'separate'; %#ok<NASGU>
    top_n = 5; %#ok<NASGU>
    make_xcorr_figures = false; %#ok<NASGU>
    make_activation_maps = false; %#ok<NASGU>
    make_roi_summaries = false; %#ok<NASGU>
    require_all_density_sources = false; %#ok<NASGU>
    continue_on_error = true; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_bold_dimred_cross_modal_coupling.m'));
end
