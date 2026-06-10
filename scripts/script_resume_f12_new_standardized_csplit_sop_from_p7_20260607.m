% Resume the F12 standardized-csplit SOP after successful P5 prerequisites.
%
% This entry intentionally skips BLP full export and P5 reduction/density.
% It exists to avoid stale workspace filters leaking from earlier run(...)
% calls and to resume only P7/P9/P8/P10 after P5 is already complete.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
cd(repo_root);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

run_p5_prereqs = false; %#ok<NASGU>
run(fullfile(repo_root, 'scripts', 'script_run_f12_new_standardized_csplit_sop_20260607.m'));
