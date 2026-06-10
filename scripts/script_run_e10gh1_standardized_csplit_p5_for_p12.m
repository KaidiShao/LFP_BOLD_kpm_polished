% Run E10gH1 standardized complex-split P5 products needed by Pipeline12.
%
% Inputs are full EDMD efun chunks exported under:
%
%   E:\DataPons_processed\derived_autodl_results_standardize\e10gh1\mlp\outputs\
%
% This script is numeric-first/headless and only targets:
%
%   complex_split_projected_vlambda_standardize
%   complex_split_projected_vlambda_standardize_rmsenv_adaptive

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

if ~exist('run_reduction', 'var') || isempty(run_reduction)
    run_reduction = true;
end
if ~exist('run_density', 'var') || isempty(run_density)
    run_density = true;
end
if ~exist('run_peak_stats', 'var') || isempty(run_peak_stats)
    run_peak_stats = true;
end
if ~exist('force_recompute_peak_stats', 'var') || isempty(force_recompute_peak_stats)
    force_recompute_peak_stats = false;
end

if ~exist('method_filter', 'var') || isempty(method_filter)
    method_filter = {'svd', 'nmf', 'mds', 'umap'};
end
if ~exist('component_count_sweep', 'var') || isempty(component_count_sweep)
    component_count_sweep = 3:8;
end
if ~exist('time_component_count_sweep', 'var') || isempty(time_component_count_sweep)
    time_component_count_sweep = component_count_sweep;
end
if ~exist('spectrum_component_count_sweep', 'var') || isempty(spectrum_component_count_sweep)
    spectrum_component_count_sweep = component_count_sweep;
end

autodl_root = 'E:\DataPons_processed\derived_autodl_results_standardize';
run_name = ['mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_', ...
    '20260522_pat40_allblp_e10gh1_seed1234_projected_vlambda_complex_split'];

% Keep P5 headless.  P12 is summary-first and does not need P5 plots here.
make_overview_plot = false;
make_state_space_plot = false;
make_consensus_state_space_plot = false;
make_spectrum_diagnostics = false;
make_top30_window_plots = false;
make_thresholded_density = false;
make_thresholded_events = false;
make_dimred_thresholded_density = false;
make_dimred_thresholded_events = false;

density_value_transform = 'abs';
lfp_activity_transform = 'rms_envelope';
lfp_activity_window_policy = 'eigenvalue_adaptive_rms_envelope';
envelope_enable = true;
envelope_policy = 'eigenvalue_adaptive_rms';
envelope_alpha = 0.35;
envelope_min_window_sec = 0.03;
envelope_max_window_sec = 1.0;
envelope_fallback_window_sec = 0.10;

fprintf('E10gH1 standardized complex-split P5 for P12\n');
fprintf('  run_reduction : %d\n', logical(run_reduction));
fprintf('  run_density   : %d\n', logical(run_density));
fprintf('  run_peak_stats: %d\n', logical(run_peak_stats));
fprintf('  autodl root   : %s\n', autodl_root);
fprintf('  run name      : %s\n', run_name);

if run_reduction
    cfg_name = 'E10gH1';
    condition_key_filter = {'e10gh1|complex_split|projected_vlambda'};
    run_name_filter = {run_name};
    condition_output_tag = 'complex_split_projected_vlambda_standardize';
    condition_tag_mode = 'condition';
    allow_summary_only_outputs = false;
    continue_on_error = true;

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_blp_eigenfunction_reduction.m'));
end

if run_density
    dataset_filter = {'e10gh1'};
    condition_filter = {'complex_split'};
    source_condition_extra_suffix = 'standardize';
    output_condition_extra_suffix = 'standardize_rmsenv_adaptive';
    activity_suffix = 'rmsenv_adaptive';
    autodl_root_override = autodl_root;
    run_name_filter_override = {run_name};
    run_raw_density = true;
    run_dimred_density = true;
    skip_existing_density = false;
    dry_run = false;

    run(fullfile(repo_root, 'scripts', ...
        'script_run_current_p5_adaptive_envelope_density_from_existing_reductions_unblocked.m'));
end

if run_peak_stats
    dataset_stems = {'e10gh1'};
    variant_filter = {'complex_split_projected_vlambda_standardize'};
    target_variant_suffix = 'rmsenv_adaptive';
    peak_mode = 'max_abs';
    peak_activity_suffix = 'rmsenv_adaptive';
    force_recompute = logical(force_recompute_peak_stats);
    save_figures = false;
    dry_run = false;
    max_tasks = [];

    run(fullfile(repo_root, 'scripts', ...
        'script_run_current_pipeline5_peak_state_all_datasets.m'));
end

fprintf('\nFinished E10gH1 standardized complex-split P5 for P12.\n');
