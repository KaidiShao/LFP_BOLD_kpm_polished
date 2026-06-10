% Run current P5 activity-magnitude outputs for datasets not blocked by P4.
%
% This runner intentionally processes one dataset at a time.  It skips datasets
% whose current BLP P4 raw/standardize outputs are known to be missing, and it
% skips the invalid/stale e10fV1 complex_split run until P4 is repaired.
%
% Optional controls, set before run(...):
%   dataset_filter = {'e10gb1','f12m01'};
%   condition_filter = {'abs','complex_split'};
%   method_filter = {'svd','nmf','mds','umap'};
%   component_count_sweep = 3:8;
%   make_state_space_plot = true;
%   make_consensus_state_space_plot = true;
%   make_spectrum_diagnostics = true;
%   make_top30_window_plots = false;
%   dry_run = true;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('dry_run', 'var') || isempty(dry_run)
    dry_run = false;
end
if ~exist('dataset_filter', 'var')
    dataset_filter = {};
end
if ~exist('condition_filter', 'var')
    condition_filter = {};
end
if ~exist('method_filter', 'var') || isempty(method_filter)
    method_filter = {'svd', 'nmf', 'mds', 'umap'};
end
if ~exist('component_count_sweep', 'var') || isempty(component_count_sweep)
    component_count_sweep = 3:8;
end
if ~exist('make_state_space_plot', 'var') || isempty(make_state_space_plot)
    make_state_space_plot = true;
end
if ~exist('make_consensus_state_space_plot', 'var') || isempty(make_consensus_state_space_plot)
    make_consensus_state_space_plot = true;
end
if ~exist('make_spectrum_diagnostics', 'var') || isempty(make_spectrum_diagnostics)
    make_spectrum_diagnostics = true;
end
if ~exist('make_top30_window_plots', 'var') || isempty(make_top30_window_plots)
    make_top30_window_plots = false;
end

dataset_specs = struct( ...
    'dataset', {'e10gb1', 'e10fV1', 'e10gh1', 'e10gw1', 'f12m01'}, ...
    'cfg_name', {'E10gb1', 'E10fV1', 'E10gH1', 'E10gW1', 'F12m01'}, ...
    'conditions', {{'abs', 'complex_split'}, {'abs'}, ...
        {'abs', 'complex_split'}, {'abs', 'complex_split'}, ...
        {'abs', 'complex_split'}} ...
    );

dataset_filter = cellstr(string(dataset_filter(:)).');
condition_filter = cellstr(string(condition_filter(:)).');
dataset_filter_lc = lower(string(dataset_filter));
condition_filter_lc = lower(string(condition_filter));

fprintf('Current P5 activity-magnitude rerun (unblocked subset)\n');
fprintf('  autodl_root: E:\\autodl_results_new\n');
fprintf('  methods    : %s\n', strjoin(cellstr(string(method_filter)), ', '));
fprintf('  k sweep    : %s\n', mat2str(component_count_sweep));
fprintf('  top30 figs : %d\n', logical(make_top30_window_plots));
fprintf('  dry_run    : %d\n', logical(dry_run));

for i_dataset = 1:numel(dataset_specs)
    spec = dataset_specs(i_dataset);
    dataset_lc = lower(string(spec.dataset));
    if ~isempty(dataset_filter_lc) && ~any(dataset_filter_lc == dataset_lc)
        continue;
    end

    conditions = spec.conditions;
    if ~isempty(condition_filter_lc)
        keep = false(size(conditions));
        for i_cond = 1:numel(conditions)
            keep(i_cond) = any(condition_filter_lc == lower(string(conditions{i_cond})));
        end
        conditions = conditions(keep);
    end
    if isempty(conditions)
        continue;
    end

    condition_key_filter = cell(1, numel(conditions)); %#ok<NASGU>
    for i_cond = 1:numel(conditions)
        condition_key_filter{i_cond} = sprintf('%s|%s|projected_vlambda', ...
            spec.dataset, conditions{i_cond});
    end

    fprintf('\n[P5] %s (%s)\n', spec.dataset, spec.cfg_name);
    fprintf('  condition_key_filter:\n');
    for i_cond = 1:numel(condition_key_filter)
        fprintf('    - %s\n', condition_key_filter{i_cond});
    end

    if dry_run
        continue;
    end

    cfg_name = spec.cfg_name; %#ok<NASGU>
    autodl_root = 'E:\autodl_results_new'; %#ok<NASGU>
    condition_tag_mode = 'condition'; %#ok<NASGU>
    density_value_transform = 'abs'; %#ok<NASGU>
    lfp_activity_transform = 'abs_magnitude'; %#ok<NASGU>
    lfp_activity_window_policy = 'samplewise_abs_no_envelope'; %#ok<NASGU>
    envelope_enable = false; %#ok<NASGU>
    envelope_policy = 'none'; %#ok<NASGU>
    component_timescale_weight_transform = 'abs'; %#ok<NASGU>
    component_timescale_top_modes = 5; %#ok<NASGU>
    make_state_space_plot = logical(make_state_space_plot); %#ok<NASGU>
    make_consensus_state_space_plot = logical(make_consensus_state_space_plot); %#ok<NASGU>
    make_spectrum_diagnostics = logical(make_spectrum_diagnostics); %#ok<NASGU>
    make_top30_window_plots = logical(make_top30_window_plots); %#ok<NASGU>
    make_thresholded_density = true; %#ok<NASGU>
    make_dimred_thresholded_density = true; %#ok<NASGU>
    make_thresholded_events = false; %#ok<NASGU>
    make_dimred_thresholded_events = false; %#ok<NASGU>
    threshold_mode = 'quantile'; %#ok<NASGU>
    threshold_ratio = 0.7; %#ok<NASGU>
    threshold_ratio_sweep = []; %#ok<NASGU>
    continue_on_error = true; %#ok<NASGU>

    run(fullfile(repo_root, 'scripts', ...
        'script_run_one_cfg_blp_eigenfunction_reduction.m'));
end

fprintf('\nFinished current P5 activity-magnitude rerun script.\n');
