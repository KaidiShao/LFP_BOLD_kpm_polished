% Convenience local chain for E10gW1 using the current cfg_E10gW1.m.
%
% Runs:
%   1. BLP streamed observables / dictionaries
%   2. BLP consensus-state top-window pipeline
%   3. BOLD observable pipeline

this_script_dir = fileparts(mfilename('fullpath'));
scripts_dir = fileparts(this_script_dir);
repo_root = fileparts(scripts_dir);
addpath(genpath(repo_root));

cfg_name = 'E10gW1';
run(fullfile(scripts_dir, 'script_preprocess_one_cfg_to_observables_streamed.m'));

cfg_name = 'E10gW1';
run(fullfile(scripts_dir, 'script_run_one_cfg_to_consensus_state_top_windows.m'));

cfg_name = 'E10gW1';
run(fullfile(scripts_dir, 'script_build_one_cfg_bold_observables.m'));
