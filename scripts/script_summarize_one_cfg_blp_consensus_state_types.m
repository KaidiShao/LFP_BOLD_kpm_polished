% Canonical single-dataset consensus-state summary stage entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_summarize_one_cfg_blp_consensus_state_types.m');
%
% Optional controls, set before run(...) if needed:
%   force_summary_recompute = false;
%   save_summary_csv = true;

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
if exist('force_summary_recompute', 'var') && ~isempty(force_summary_recompute)
    params.summary_params.force_recompute = logical(force_summary_recompute);
end
if exist('save_summary_csv', 'var') && ~isempty(save_summary_csv)
    params.summary_params.save_csv = logical(save_summary_csv);
end

S = summarize_blp_consensus_state_types( ...
    cfg, ...
    output_root, ...
    C, ...
    params.summary_params, ...
    source_consensus_file);

disp(S.summary_table);
fprintf('Saved consensus-state type summary to:\n  %s\n', S.save_file);
if ~isempty(S.csv_file)
    fprintf('Saved CSV summary to:\n  %s\n', S.csv_file);
end
