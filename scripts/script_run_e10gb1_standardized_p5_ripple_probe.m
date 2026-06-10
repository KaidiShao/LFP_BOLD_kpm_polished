% Run an e10gb1 standardized-BLP P5 ripple-selectivity probe.
%
% This is an experimental branch, not the current P5 mainline.  It only uses
% the two standardized BLP P4 runs identified in the 2026-05-23 raw-vs-std
% monitor report and writes separate P5 condition tags:
%
%   abs_projected_vlambda_standardize
%   complex_split_projected_vlambda_standardize
%
% Peak statistics are written with the current RMS-envelope activity policy
% under:
%
%   pipeline5_eigenfunction_peaks_by_state_rmsenv_adaptive\...

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

if ~exist('run_reduction', 'var') || isempty(run_reduction)
    run_reduction = true;
end
if ~exist('run_peak_stats', 'var') || isempty(run_peak_stats)
    run_peak_stats = true;
end
if ~exist('force_recompute_reduction', 'var') || isempty(force_recompute_reduction)
    force_recompute_reduction = false; %#ok<NASGU> % reserved for future use
end
if ~exist('force_recompute_peak_stats', 'var') || isempty(force_recompute_peak_stats)
    force_recompute_peak_stats = false;
end

method_filter = {'svd', 'nmf', 'mds', 'umap'};
component_count_sweep = 3:8;
time_component_count_sweep = 3:8;
spectrum_component_count_sweep = 3:8;

% Keep this probe numeric-first and headless.
make_overview_plot = false;
make_state_space_plot = false;
make_consensus_state_space_plot = false;
make_spectrum_diagnostics = false;
make_top30_window_plots = false;

% Reduction-stage density outputs are not needed to answer the immediate
% ripple-selectivity question; peak statistics below recompute activity from
% the reduction MAT files using RMS-envelope activity.
make_thresholded_density = false;
make_thresholded_events = false;
make_dimred_thresholded_density = false;
make_dimred_thresholded_events = false;

% Keep all activity metadata explicit for provenance.
density_value_transform = 'abs';
lfp_activity_transform = 'rms_envelope';
lfp_activity_window_policy = 'eigenvalue_adaptive_rms_envelope';
envelope_enable = true;
envelope_policy = 'eigenvalue_adaptive_rms';
envelope_alpha = 0.35;
envelope_min_window_sec = 0.03;
envelope_max_window_sec = 1.0;
envelope_fallback_window_sec = 0.10;

autodl_root = 'E:\DataPons_processed\derived_autodl_results_standardize';
allow_summary_only_outputs = false;
cfg_name = 'E10gb1';
continue_on_error = true;

run_specs = struct( ...
    'condition_key', { ...
        'e10gb1|abs|projected_vlambda', ...
        'e10gb1|complex_split|projected_vlambda'}, ...
    'run_name', { ...
        'mlp_obs_blp_vlambda_abs_stdT_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gb1_seed1234_projected_vlambda_abs', ...
        'mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_b2000_i2_pat40_20260522_pat40_allblp_e10gb1_seed1234_projected_vlambda_complex_split'}, ...
    'output_tag', { ...
        'abs_projected_vlambda_standardize', ...
        'complex_split_projected_vlambda_standardize'} ...
    );

fprintf('e10gb1 standardized P5 ripple probe\n');
fprintf('  run_reduction : %d\n', logical(run_reduction));
fprintf('  run_peak_stats: %d\n', logical(run_peak_stats));
fprintf('  methods       : %s\n', strjoin(method_filter, ', '));
fprintf('  k sweep       : %s\n', mat2str(component_count_sweep));

if run_reduction
    for i_spec = 1:numel(run_specs)
        spec = run_specs(i_spec);
        fprintf('\n[P5 reduction probe %d/%d] %s\n', ...
            i_spec, numel(run_specs), spec.output_tag);
        fprintf('  source run: %s\n', spec.run_name);

        condition_key_filter = {spec.condition_key};
        run_name_filter = {spec.run_name};
        condition_output_tag = spec.output_tag;
        condition_tag_mode = 'condition';

        run(fullfile(repo_root, 'scripts', ...
            'script_run_one_cfg_blp_eigenfunction_reduction.m'));
    end
end

if run_peak_stats
    fprintf('\n[P5 peak-stat probe] standardized variants -> RMS-envelope stats\n');
    dataset_stems = {'e10gb1'};
    variant_filter = { ...
        'abs_projected_vlambda_standardize', ...
        'complex_split_projected_vlambda_standardize'};
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

fprintf('\nFinished e10gb1 standardized P5 ripple probe.\n');
