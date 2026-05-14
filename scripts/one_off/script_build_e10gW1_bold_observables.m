% One-off fixed-dataset wrapper for the BOLD observable pipeline.
this_script_dir = fileparts(mfilename('fullpath'));
scripts_dir = fileparts(this_script_dir);
repo_root = fileparts(scripts_dir);
addpath(genpath(repo_root));

cfg_name = 'E10gW1';
run(fullfile(scripts_dir, 'script_build_one_cfg_bold_observables.m'));
