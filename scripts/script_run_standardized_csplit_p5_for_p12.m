% Run standardized complex-split P5 products needed by Pipeline12.
%
% This is the generic version of the earlier E10gH1-only entry.  It expects
% full EDMD output chunks exported under:
%
%   E:\DataPons_processed\derived_autodl_results_standardize\<dataset>\mlp\outputs\
%
% Default target datasets are the newly stable standardized complex-split
% runs.  Override before run(...) if needed:
%
%   target_dataset_stems = {'e10fV1', 'e10gw1'};
%   run_reduction = true;
%   run_density = true;
%   run_peak_stats = true;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

if ~exist('target_dataset_stems', 'var') || isempty(target_dataset_stems)
    target_dataset_stems = {'e10fV1', 'e10gw1'};
end
if ~exist('run_reduction', 'var') || isempty(run_reduction)
    run_reduction = true;
end
if ~exist('run_density', 'var') || isempty(run_density)
    run_density = true;
end
if ~exist('run_peak_stats', 'var') || isempty(run_peak_stats)
    run_peak_stats = true;
end
if ~exist('skip_existing_reduction', 'var') || isempty(skip_existing_reduction)
    skip_existing_reduction = true;
end
if ~exist('skip_existing_density', 'var') || isempty(skip_existing_density)
    skip_existing_density = true;
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

if ~exist('autodl_root_override', 'var') || isempty(autodl_root_override)
    autodl_root_override = 'E:\DataPons_processed\derived_autodl_results_standardize';
end
if ~exist('allow_summary_only_outputs', 'var') || isempty(allow_summary_only_outputs)
    allow_summary_only_outputs = false;
end

autodl_root = char(string(autodl_root_override));
source_condition_tag = 'complex_split_projected_vlambda_standardize';
target_condition_tag = 'complex_split_projected_vlambda_standardize_rmsenv_adaptive';

all_specs = local_standardized_csplit_specs();
specs = local_filter_specs(all_specs, target_dataset_stems);
if isempty(specs)
    error('No standardized complex-split specs matched target_dataset_stems.');
end

% Keep P5 headless.  P12 is numeric/summary-first.
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

fprintf('Standardized complex-split P5 for P12\n');
fprintf('  run_reduction          : %d\n', logical(run_reduction));
fprintf('  run_density            : %d\n', logical(run_density));
fprintf('  run_peak_stats         : %d\n', logical(run_peak_stats));
fprintf('  skip_existing_reduction: %d\n', logical(skip_existing_reduction));
fprintf('  skip_existing_density  : %d\n', logical(skip_existing_density));
fprintf('  summary-only outputs   : %d\n', logical(allow_summary_only_outputs));
fprintf('  autodl root            : %s\n', autodl_root);
fprintf('  datasets               : %s\n', strjoin({specs.dataset}, ', '));

if run_reduction
    processed_root = io_project.get_project_processed_root();
    for i_spec = 1:numel(specs)
        spec = specs(i_spec);
        if skip_existing_reduction && local_has_complete_reduction( ...
                processed_root, spec.dataset, source_condition_tag, ...
                method_filter, component_count_sweep)
            fprintf('\n[reduction skip] %s already has all requested method-k MAT files.\n', ...
                spec.dataset);
            continue;
        end

        fprintf('\n[reduction] %s\n', spec.dataset);
        cfg_name = spec.cfg_name; %#ok<NASGU>
        condition_key_filter = {lower(sprintf('%s|complex_split|projected_vlambda', ...
            spec.dataset))}; %#ok<NASGU>
        run_name_filter = {spec.run_name}; %#ok<NASGU>
        condition_output_tag = source_condition_tag; %#ok<NASGU>
        condition_tag_mode = 'condition'; %#ok<NASGU>
        allow_summary_only_outputs = logical(allow_summary_only_outputs); %#ok<NASGU>
        continue_on_error = true; %#ok<NASGU>

        run(fullfile(repo_root, 'scripts', ...
            'script_run_one_cfg_blp_eigenfunction_reduction.m'));
    end
end

if run_density
    fprintf('\n[density] %s\n', strjoin({specs.dataset}, ', '));
    dataset_filter = {specs.dataset}; %#ok<NASGU>
    condition_filter = {'complex_split'}; %#ok<NASGU>
    source_condition_extra_suffix = 'standardize'; %#ok<NASGU>
    output_condition_extra_suffix = 'standardize_rmsenv_adaptive'; %#ok<NASGU>
    activity_suffix = 'rmsenv_adaptive'; %#ok<NASGU>
    autodl_root_override = autodl_root; %#ok<NASGU>
    run_name_filter_override = {specs.run_name}; %#ok<NASGU>
    dataset_specs_override = repmat(struct('dataset', '', 'cfg_name', '', ...
        'conditions', {{'complex_split'}}), 1, numel(specs)); %#ok<NASGU>
    for i_override = 1:numel(specs)
        dataset_specs_override(i_override).dataset = specs(i_override).dataset;
        dataset_specs_override(i_override).cfg_name = specs(i_override).cfg_name;
        dataset_specs_override(i_override).conditions = {'complex_split'};
    end
    run_raw_density = true; %#ok<NASGU>
    run_dimred_density = true; %#ok<NASGU>
    allow_summary_only_outputs = logical(allow_summary_only_outputs); %#ok<NASGU>
    skip_existing_density = logical(skip_existing_density); %#ok<NASGU>
    dry_run = false; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_current_p5_adaptive_envelope_density_from_existing_reductions_unblocked.m'));
end

if run_peak_stats
    fprintf('\n[peak stats] %s\n', strjoin({specs.dataset}, ', '));
    dataset_stems = {specs.dataset}; %#ok<NASGU>
    variant_filter = {source_condition_tag}; %#ok<NASGU>
    target_variant_suffix = 'rmsenv_adaptive'; %#ok<NASGU>
    peak_mode = 'max_abs'; %#ok<NASGU>
    peak_activity_suffix = 'rmsenv_adaptive'; %#ok<NASGU>
    force_recompute = logical(force_recompute_peak_stats); %#ok<NASGU>
    save_figures = false; %#ok<NASGU>
    dry_run = false; %#ok<NASGU>
    max_tasks = []; %#ok<NASGU>
    continue_on_error = true; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_current_pipeline5_peak_state_all_datasets.m'));
end

fprintf('\nFinished standardized complex-split P5 for P12.\n');


function specs = local_standardized_csplit_specs()
template = ['mlp_obs_blp_vlambda_complex_split_stdComplexPair_l1e4_r1e3_', ...
    'b2000_i2_pat40_20260522_pat40_allblp_%s_seed1234_projected_vlambda_complex_split'];
datasets = {'e10gb1', 'e10fV1', 'e10gh1', 'e10gw1', ...
    'f12m01', 'f12m02', 'f12m03', 'f12m05', ...
    'k13m17', 'k13m18', 'k13m19', 'k13m20', 'k13m21', 'k13m23'};
cfg_names = {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', ...
    'F12m01', 'F12m02', 'F12m03', 'F12m05', ...
    'K13m17', 'K13m18', 'K13m19', 'K13m20', 'K13m21', 'K13m23'};
specs = repmat(struct('dataset', '', 'cfg_name', '', 'run_name', ''), ...
    1, numel(datasets));
for i = 1:numel(datasets)
    specs(i).dataset = datasets{i};
    specs(i).cfg_name = cfg_names{i};
    specs(i).run_name = sprintf(template, datasets{i});
end
end


function specs = local_filter_specs(all_specs, target_dataset_stems)
targets = lower(string(cellstr(string(target_dataset_stems(:)).')));
keep = false(size(all_specs));
for i = 1:numel(all_specs)
    keep(i) = any(lower(string(all_specs(i).dataset)) == targets);
end
specs = all_specs(keep);
end


function tf = local_has_complete_reduction(processed_root, dataset, condition_tag, ...
        method_filter, component_count_sweep)
tf = true;
stage_dir = fullfile(processed_root, dataset, ...
    io_project.get_pipeline_stage_name(5, 'eigenfunction_reduction'), ...
    condition_tag);
if exist(stage_dir, 'dir') ~= 7
    tf = false;
    return;
end

methods = cellstr(string(method_filter(:)).');
for i_method = 1:numel(methods)
    method = lower(char(string(methods{i_method})));
    for k = component_count_sweep(:).'
        method_tag = sprintf('%s_k%02d', method, k);
        mat_dir = fullfile(stage_dir, method_tag, 'mat');
        if exist(mat_dir, 'dir') ~= 7 || isempty(dir(fullfile(mat_dir, '*.mat')))
            tf = false;
            return;
        end
    end
end
end
