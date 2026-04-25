this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

cfg = cfg_E10fV1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = get_project_processed_root();
params = struct();
params.window_params = struct('top_k', 30, 'window_length_samples', 6000, ...
    'window_mode', 'global', 'keep_partial_window', false, 'save_csv', true, 'force_recompute', false);
params.plot_params = struct();

out = run_blp_pipeline_to_state_diversity_top_windows(cfg, output_root, params);
disp(out.W.top_windows_table);
