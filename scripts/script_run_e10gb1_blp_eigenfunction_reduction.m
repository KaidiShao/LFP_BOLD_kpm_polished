% Fixed-dataset example wrapper for pipeline 5.

this_script_dir = fileparts(mfilename('fullpath'));
cfg_name = 'E10gb1';
run(fullfile(this_script_dir, 'script_run_one_cfg_blp_eigenfunction_reduction.m'));
