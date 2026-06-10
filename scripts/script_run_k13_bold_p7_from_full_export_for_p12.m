% Run P7 for K13 BOLD runs exported from summary-only checkpoints.
%
% This consumes derived full EDMD chunks under:
%   E:\DataPons_processed\derived_autodl_results_bold_full_export
%
% It keeps figures off; P12 needs BOLD_POST numeric outputs and deconv efuns.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    cfg_names = {'K13m17', 'K13m23'};
end
if ~exist('main_observable_modes', 'var') || isempty(main_observable_modes)
    main_observable_modes = {'HP_svd100', 'global_svd100', 'gsvd100_ds', 'roi_mean'};
end
if ~exist('derived_bold_root', 'var') || isempty(derived_bold_root)
    derived_bold_root = 'E:\DataPons_processed\derived_autodl_results_bold_full_export';
end

for i_cfg = 1:numel(cfg_names)
    cfg_name = cfg_names{i_cfg}; %#ok<NASGU>
    observable_modes = main_observable_modes; %#ok<NASGU>
    autodl_roots = {derived_bold_root}; %#ok<NASGU>
    force_recompute = false; %#ok<NASGU>
    make_main_plot = false; %#ok<NASGU>
    compute_deconv = true; %#ok<NASGU>
    make_deconv_plot = false; %#ok<NASGU>
    make_timescale_plot = false; %#ok<NASGU>
    make_intrinsic_activation_maps = false; %#ok<NASGU>
    make_intrinsic_roi_summary = false; %#ok<NASGU>
    headless = true; %#ok<NASGU>
    close_figures_after_each_run = true; %#ok<NASGU>
    continue_on_error = true; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_bold_reskoopnet_postprocessing.m'));
end
