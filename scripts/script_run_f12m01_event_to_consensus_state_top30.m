% Run F12m01 BLP event pipeline through consensus-state top30 windows.
%
% This wrapper is cache-first:
%   1. reuse existing event_detection result when present
%   2. update/reuse event_density
%   3. update/reuse consensus_states
%   4. update/reuse consensus_state_summary
%   5. recompute the canonical consensus-state-diversity top30 windows
%
% Canonical top-window settings:
%   window_length_samples = 6000
%   window_mode = global
%   top_k = 30

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
    % Older MATLAB releases do not expose this graphics default.
end
set(groot, 'defaultFigureVisible', 'off');

cfg = cfg_F12m01();
output_root = get_project_processed_root();

opts = struct();
opts.force_event_recompute = false;
opts.force_density_recompute = false;
opts.force_consensus_recompute = false;
opts.force_summary_recompute = true;
opts.force_window_recompute = true;
opts.compute_density = true;

if exist('f12m01_event_top30_override', 'var') == 1 && isstruct(f12m01_event_top30_override)
    opts = local_merge_struct(opts, f12m01_event_top30_override);
end

fprintf('Running F12m01 event-to-consensus-state-top30 pipeline.\n');
fprintf('Dataset: %s | file_stem: %s\n', cfg.dataset_id, cfg.file_stem);
fprintf('Output root:\n  %s\n', output_root);

event_params = struct();
event_params.passband = [2, 15; 30, 90; 90, 190];
event_params.band_labels = {'theta', 'gamma', 'ripple'};
event_params.L_start_range = [151, 101, 51];
event_params.L_extract_range = [301, 201, 101];
event_params.ThresRatio_range = [3.5, 4, 4];
event_params.input_normalization = 'zscore_per_channel';
event_params.force_recompute = logical(opts.force_event_recompute);

if opts.force_event_recompute
    [R, source_event_file] = local_compute_events_from_raw(cfg, output_root, event_params);
else
    try
        [R, source_event_file] = load_event_results(cfg, output_root, []);
        fprintf('[stage] Reusing event detection:\n  %s\n', source_event_file);
    catch ME
        fprintf('[stage] Existing event detection not found or not loadable: %s\n', ME.message);
        [R, source_event_file] = local_compute_events_from_raw(cfg, output_root, event_params);
    end
end

E = [];
if opts.compute_density
    fprintf('[stage] Updating event density...\n');
    density_params = struct();
    density_params.bin_sec = 2;
    density_params.smooth_sigma_sec = 2;
    density_params.force_recompute = logical(opts.force_density_recompute);

    E = compute_blp_event_density(cfg, output_root, R, density_params, source_event_file);
    fprintf('[stage] Event density ready:\n  %s\n', E.save_file);
end

fprintf('[stage] Updating consensus states...\n');
consensus_params = struct();
consensus_params.min_channel_count = [];
consensus_params.require_region_presence = false;
consensus_params.required_regions = {'hp', 'pl'};
consensus_params.force_recompute = logical(opts.force_consensus_recompute);

C = compute_blp_consensus_states(cfg, output_root, R, consensus_params, source_event_file);
fprintf('[stage] Consensus states ready:\n  %s\n', C.save_file);

fprintf('[stage] Updating consensus-state summary...\n');
summary_params = struct();
summary_params.save_csv = true;
summary_params.force_recompute = logical(opts.force_summary_recompute);

S = summarize_blp_consensus_state_types(cfg, output_root, C, summary_params, C.save_file);
fprintf('[stage] Consensus-state summary ready:\n  %s\n', S.save_file);
if ~isempty(S.csv_file)
    fprintf('[stage] Consensus-state summary CSV:\n  %s\n', S.csv_file);
end

fprintf('[stage] Computing consensus-state-diversity top30 windows...\n');
window_params = struct();
window_params.window_length_samples = 6000;
window_params.window_mode = 'global';
window_params.keep_partial_window = false;
window_params.top_k = 30;
window_params.save_csv = true;
window_params.force_recompute = logical(opts.force_window_recompute);

W = analyze_blp_consensus_state_diversity_windows( ...
    cfg, ...
    output_root, ...
    C, ...
    window_params, ...
    C.save_file);

fprintf('[stage] State-diversity windows ready:\n  %s\n', W.save_file);
if ~isempty(W.csv_file)
    fprintf('[stage] All-window CSV:\n  %s\n', W.csv_file);
end
if ~isempty(W.top_csv_file)
    fprintf('[stage] Top30 CSV:\n  %s\n', W.top_csv_file);
end

disp('Top consensus-state-diversity windows:');
disp(W.top_windows_table);

out = struct();
out.cfg = cfg;
out.output_root = output_root;
out.R = R;
out.E = E;
out.C = C;
out.S = S;
out.W = W;
out.opts = opts;

fprintf('F12m01 event-to-consensus-state-top30 pipeline finished.\n');

clear f12m01_event_top30_override


function [R, source_event_file] = local_compute_events_from_raw(cfg, output_root, event_params)
fprintf('[stage] Loading raw BLP for event detection...\n');
D = load_blp_dataset(cfg);
fprintf('[stage] Loaded concatenated raw data: [%d, %d]\n', size(D.data, 1), size(D.data, 2));

fprintf('[stage] Computing event detection...\n');
R = compute_blp_bandpass_events(D, cfg, output_root, event_params);
source_event_file = R.save_file;
fprintf('[stage] Event detection ready:\n  %s\n', source_event_file);
end


function dst = local_merge_struct(dst, src)
names = fieldnames(src);
for i = 1:numel(names)
    dst.(names{i}) = src.(names{i});
end
end
