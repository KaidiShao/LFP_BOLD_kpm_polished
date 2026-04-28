% Canonical interactive event-diversity plotting entry.
%
% Usage from MATLAB:
%   cfg_name = 'F12m01';
%   run('scripts/script_plot_top_consensus_event_diversity_windows.m');
%
% Optional controls, set before run(...) if needed:
%   window_length_samples = 6000;
%   window_mode = 'global';
%   window_result_file = '';
%   skip_existing = true;
%   plot_show_consensus = false;
%
% Role:
%   - canonical single-dataset plotting entry for event-diversity windows
%   - loads a saved event-diversity result, then delegates figure export to
%     the canonical event-diversity plotting function
%
% For the fixed e10gb1 example, use:
%   scripts/script_plot_top_consensus_event_diversity_windows_e10gb1.m

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''F12m01'';']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

othercolor_root = 'D:\Onedrive\util_functions\othercolor\';
if exist(othercolor_root, 'dir') == 7
    addpath(genpath(othercolor_root));
end

close all force;

cfg_fun = ['cfg_' char(string(cfg_name))];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
cfg.plot.trace_scale = 0.18;
cfg.plot.trace_clip = 4;
cfg.plot.within_gap = 1.4;
cfg.plot.between_gap = 2.2;
cfg.plot.trace_linewidth = 0.4;

output_root = io_project.get_project_processed_root();
branch_params = build_blp_event_diversity_pipeline_params();
if exist('window_length_samples', 'var') && ~isempty(window_length_samples)
    branch_params.window_params.window_length_samples = window_length_samples;
end
if exist('window_mode', 'var') && ~isempty(window_mode)
    branch_params.window_params.window_mode = char(string(window_mode));
end

if ~exist('window_result_file', 'var') || isempty(window_result_file)
    save_tag = local_build_save_tag( ...
        cfg.file_stem, ...
        branch_params.window_params.window_length_samples, ...
        branch_params.window_params.window_mode);
    window_result_file = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', [save_tag, '.mat']);
end

if exist(window_result_file, 'file') ~= 2
    error('Event-diversity result file does not exist:\n  %s', window_result_file);
end

S = load(window_result_file, 'W');
W = S.W;

params = struct();
params.save_dir = fullfile(output_root, cfg.file_stem, 'event_diversity_windows', 'top_window_plots');
params.save_png = true;
params.close_after_save = true;
params.skip_existing = true;
params.freq_range_to_plot = [0, 250];
params.color_limits = [];
params.show_consensus = false;
if exist('skip_existing', 'var') && ~isempty(skip_existing)
    params.skip_existing = logical(skip_existing);
end
if exist('plot_show_consensus', 'var') && ~isempty(plot_show_consensus)
    params.show_consensus = logical(plot_show_consensus);
end

if exist('othercolor', 'file') == 2
    params.spec_colormap = flipud(othercolor('Spectral10'));
else
    params.spec_colormap = flipud(turbo(256));
end

P = export_top_consensus_event_diversity_window_plots(cfg, output_root, W, params);
fprintf('Saved %d top-window plots to:\n  %s\n', height(W.top_windows_table), P.save_dir);


function tag = local_build_save_tag(file_stem, window_length_samples, window_mode)
tag = sprintf('%s_event_diversity_windows_%dsamp', file_stem, window_length_samples);
if strcmpi(char(string(window_mode)), 'global')
    tag = sprintf('%s_globalwin', tag);
end
end
