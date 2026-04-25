% Analyze reduced temporal-component peaks inside e10gb1 consensus-state windows.
%
% This script expects that script_run_eigenfunction_reduction_minimal.m has
% already saved a compact eigenfunction-reduction result.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));
set(groot, 'defaultFigureVisible', 'off');

cfg = cfg_eigenfunction_reduction_minimal();
output_root = cfg.dataset.processed_root;

result_file = local_find_latest_result_file(cfg.save.dir);
fprintf('Loading eigenfunction reduction result:\n  %s\n', result_file);
S = load(result_file, 'result');
if ~isfield(S, 'result')
    error('Result file %s does not contain variable result.', result_file);
end
result = S.result;

consensus_loader_cfg = struct();
consensus_loader_cfg.file_stem = cfg.dataset.name;
[C, source_consensus_file] = load_consensus_state_results( ...
    consensus_loader_cfg, output_root, []);
fprintf('Loading consensus states:\n  %s\n', source_consensus_file);

params = struct();
params.component_source = 'smooth_if_available';
params.state_start_idx = cfg.viz.state_space_consensus.state_start_idx;
params.peak_mode = 'max';       % 'max' | 'max_abs'
params.baseline_mode = 'matched';
params.baseline_code = 0;
params.baseline_label = 'baseline';
params.min_window_samples = 1;
params.baseline_random_seed = 1;
params.alpha = 0.05;
params.save_dir = fullfile(cfg.output.root, 'peaks');
params.save_tag = sprintf('%s_peaks', cfg.dataset.name);
params.save_results = true;
params.write_csv = true;
params.save_figures = true;
params.close_figures = true;

A = analyze_eigenfunction_component_peaks_by_consensus_state( ...
    result, C, params, result_file, source_consensus_file);

fprintf('\nSaved component-peak analysis:\n');
if isfield(A.save_paths, 'main_mat')
    fprintf('  MAT   : %s\n', A.save_paths.main_mat);
end
if isfield(A.save_paths, 'stats_csv')
    fprintf('  Stats : %s\n', A.save_paths.stats_csv);
end
if isfield(A.save_paths, 'event_peak_csv')
    fprintf('  Peaks : %s\n', A.save_paths.event_peak_csv);
end
if isfield(A.save_paths, 'peak_distribution_png')
    fprintf('  Fig 1 : %s\n', A.save_paths.peak_distribution_png);
end
if isfield(A.save_paths, 'mean_peak_heatmap_png')
    fprintf('  Fig 2 : %s\n', A.save_paths.mean_peak_heatmap_png);
end
if isfield(A.save_paths, 'baseline_effect_heatmap_png')
    fprintf('  Fig 3 : %s\n', A.save_paths.baseline_effect_heatmap_png);
end

disp(A.stats_table);


function result_file = local_find_latest_result_file(result_dir)
if exist(result_dir, 'dir') ~= 7
    error(['Result directory does not exist:\n  %s\n', ...
        'Run script_run_eigenfunction_reduction_minimal.m first.'], result_dir);
end

L = dir(fullfile(result_dir, '*.mat'));
if isempty(L)
    error(['No eigenfunction reduction MAT files were found in:\n  %s\n', ...
        'Run script_run_eigenfunction_reduction_minimal.m first.'], result_dir);
end

[~, idx] = max([L.datenum]);
result_file = fullfile(L(idx).folder, L(idx).name);
end
