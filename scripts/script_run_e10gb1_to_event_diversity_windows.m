% Example single-dataset event-diversity quick-check script.
%
% This is not the reusable entry. For running the same quick-check branch on
% one arbitrary cfg_*.m definition, use:
%   scripts/script_run_one_cfg_to_event_diversity_windows.m
%
% For batch execution across cfg_*.m entries, use:
%   scripts/script_run_cfgs_to_event_diversity_windows.m

cfg_name = 'E10gb1';
run('scripts/script_run_one_cfg_to_event_diversity_windows.m');
