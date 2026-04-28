% Canonical e10gb1 event-diversity plotting wrapper.
% This script always uses the 6000-sample global-window diversity result and
% delegates figure export to the canonical event-diversity plotting function.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

cfg = cfg_E10gb1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = io_project.get_project_processed_root();
window_result_file = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', ...
    [cfg.file_stem, '_event_diversity_windows_6000samp_globalwin.mat']);

S = load(window_result_file, 'W');
W = S.W;
local_validate_e10gb1_diversity_result(W, window_result_file);
params = struct();
params.save_dir = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', 'top_window_plots');
params.save_png = true;
params.close_after_save = true;
params.skip_existing = false;
params.fallback_size_px = [4979, 2888];
params.reference_dir = fullfile(params.save_dir, 'ref_figs');
params.freq_range_to_plot = [0, 250];
params.color_limits = [];
params.show_consensus = true;

if exist('othercolor', 'file') == 2
    params.spec_colormap = flipud(othercolor('Spectral10'));
else
    params.spec_colormap = flipud(turbo(256));
end

P = export_top_consensus_event_diversity_window_plots(cfg, output_root, W, params);
fprintf('Saved %d top-window plots to:\n  %s\n', height(W.top_windows_table), P.save_dir);


function local_validate_e10gb1_diversity_result(W, window_result_file)
if ~isfield(W, 'window_mode') || ~strcmpi(char(string(W.window_mode)), 'global')
    error(['%s is not a global-window diversity result. ', ...
        'Use the canonical 6000-sample global-window file for e10gb1.'], window_result_file);
end

if ~isfield(W, 'window_length_samples') || double(W.window_length_samples) ~= 6000
    error(['%s does not use 6000-sample windows. ', ...
        'Use the canonical 6000-sample global-window file for e10gb1.'], window_result_file);
end

required_vars = {'global_window_idx', 'start_session_id', 'end_session_id'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, W.top_windows_table.Properties.VariableNames)
        error('Top-window table in %s is missing required column "%s".', ...
            window_result_file, required_vars{i});
    end
end
end
