% Canonical single-dataset state-diversity analysis stage entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_analyze_one_cfg_blp_consensus_state_diversity_windows.m');
%
% Optional controls, set before run(...) if needed:
%   force_window_recompute = false;
%   window_length_samples = 6000;
%   top_k = 30;

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_fun = ['cfg_' char(string(cfg_name))];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
output_root = io_project.get_project_processed_root();

[C, source_consensus_file] = io_results.load_consensus_state_results(cfg, output_root, []);

params = build_blp_consensus_state_pipeline_params();
if exist('force_window_recompute', 'var') && ~isempty(force_window_recompute)
    params.window_params.force_recompute = logical(force_window_recompute);
end
if exist('window_length_samples', 'var') && ~isempty(window_length_samples)
    params.window_params.window_length_samples = double(window_length_samples);
end
if exist('top_k', 'var') && ~isempty(top_k)
    params.window_params.top_k = double(top_k);
end

W = analyze_blp_consensus_state_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    params.window_params, ...
    source_consensus_file);

disp(W.top_windows_table);
fprintf('Saved per-window consensus-state-diversity summary to:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('Saved all-window CSV to:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('Saved top-window CSV to:\n  %s\n', W.top_csv_file);
end
