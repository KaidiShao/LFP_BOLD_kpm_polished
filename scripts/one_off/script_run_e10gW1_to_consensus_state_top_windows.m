% One-off fixed-dataset wrapper for the consensus-state top-window pipeline.
this_script_dir = fileparts(mfilename('fullpath'));
scripts_dir = fileparts(this_script_dir);
repo_root = fileparts(scripts_dir);
addpath(genpath(repo_root));

cfg_name = 'E10gW1';
run(fullfile(scripts_dir, 'script_run_one_cfg_to_consensus_state_top_windows.m'));
