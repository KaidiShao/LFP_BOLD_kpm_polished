% Canonical e10gb1 diversity-reference export entry point.
% This script always reads:
%   e10gb1_event_diversity_windows_6000samp_globalwin_top.csv
% and exports matching reference PNGs into:
%   E:\DataPons_processed\e10gb1\event_diversity_windows\top_window_plots\ref_figs

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

output_root = get_project_processed_root();

src_dir = ['D:\Onedrive\ICPBR\Alberta\koopman_events\code\test\', ...
    'test63_event_labeling_4th_round\figs_new'];
top_csv = fullfile(output_root, 'e10gb1', 'event_diversity_windows', ...
    'e10gb1_event_diversity_windows_6000samp_globalwin_top.csv');
ref_dir = fullfile(output_root, 'e10gb1', 'event_diversity_windows', ...
    'top_window_plots', 'ref_figs');

if exist(src_dir, 'dir') ~= 7
    error('Source fig directory does not exist:\n  %s', src_dir);
end

if exist(top_csv, 'file') ~= 2
    error('Top diversity CSV does not exist:\n  %s', top_csv);
end

if exist(ref_dir, 'dir') ~= 7
    mkdir(ref_dir);
end

T = readtable(top_csv, 'TextType', 'string');
if isempty(T)
    error('Top diversity CSV is empty:\n  %s', top_csv);
end
local_validate_top_csv(T, top_csv);

png_files = strings(height(T), 1);
source_figs = strings(height(T), 1);
global_window_idx = zeros(height(T), 1);

for i = 1:height(T)
    if ismember('global_window_idx', T.Properties.VariableNames)
        gidx = double(T.global_window_idx(i));
    else
        error('Top diversity CSV does not contain global_window_idx. Recompute with window_mode=''global'' first.');
    end
    global_window_idx(i) = gidx;

    fig_name = sprintf('uncorrected_signal_wavelet_koopman_%d_realimag_updated.fig', gidx);
    fig_path = fullfile(src_dir, fig_name);
    if exist(fig_path, 'file') ~= 2
        error('Missing reference fig for rank %d:\n  %s', i, fig_path);
    end

    out_name = sprintf('rank_%02d_globalwin_%03d_ref.png', ...
        double(T.diversity_rank(i)), gidx);
    out_path = fullfile(ref_dir, out_name);

    h = openfig(fig_path, 'invisible');
    drawnow;
    exportgraphics(h, out_path, 'Resolution', 220);
    close(h);

    source_figs(i) = string(fig_path);
    png_files(i) = string(out_path);
    fprintf('[%d/%d] Exported %s\n', i, height(T), out_name);
end

manifest = T;
manifest.global_window_idx = global_window_idx;
manifest.source_fig = source_figs;
manifest.ref_png = png_files;
writetable(manifest, fullfile(ref_dir, 'ref_manifest.csv'));

fprintf('Saved %d reference PNGs to:\n  %s\n', height(T), ref_dir);


function local_validate_top_csv(T, top_csv)
required_vars = {'diversity_rank', 'global_window_idx', 'global_start_idx', 'global_end_idx'};
for i = 1:numel(required_vars)
    if ~ismember(required_vars{i}, T.Properties.VariableNames)
        error('Top diversity CSV %s is missing required column "%s".', top_csv, required_vars{i});
    end
end
end
