% Maintenance runner for the currently missing completed-model conditions.
%
% This intentionally skips the old event-diversity branch. It runs the
% top-state-diversity-window postprocessing stages:
%   main postprocess, timescale diagnostics, deconv residuals,
%   window-normalized deconv residuals, and full-time SPKT/MUA cross-correlation.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_names = {'E10gH1', 'F12m01'};

% Remaining conditions after the 2026-04-23 E10/F12 postprocessing runs.
% E10gH1 abs-vlambda was exported later, so it still needs this stage.
% Keep abs-kv in the list: its top-window plots already exist, but the
% SPKT/MUA cross-correlation stages still need to be generated. Existing
% window PNGs are skipped by the canonical single-cfg script.
% Format: dataset|observable_mode|residual_form
condition_key_filter = { ...
    'e10gh1|abs|projected_vlambda', ...
    'f12m01|abs|projected_kv', ...
    'f12m01|abs|projected_vlambda', ...
    'f12m01|complex_split|projected_kv', ...
    'f12m01|complex_split|projected_vlambda'};

fprintf('Running missing completed-model postprocessing conditions.\n');
fprintf('Condition filter:\n');
for i = 1:numel(condition_key_filter)
    fprintf('  - %s\n', condition_key_filter{i});
end

run(fullfile(this_script_dir, 'script_run_cfgs_blp_eigenfunction_postprocessing.m'));

n_window_rows = 0;
n_spkt_rows = 0;
n_mua_rows = 0;
for i_cfg = 1:numel(batch_manifests)
    if ~isfield(batch_manifests(i_cfg), 'manifest') || isempty(batch_manifests(i_cfg).manifest)
        continue;
    end
    manifest_i = batch_manifests(i_cfg).manifest;
    if isfield(manifest_i, 'table') && ~isempty(manifest_i.table)
        n_window_rows = n_window_rows + height(manifest_i.table);
    end
    if isfield(manifest_i, 'spkt_table') && ~isempty(manifest_i.spkt_table)
        n_spkt_rows = n_spkt_rows + height(manifest_i.spkt_table);
    end
    if isfield(manifest_i, 'mua_table') && ~isempty(manifest_i.mua_table)
        n_mua_rows = n_mua_rows + height(manifest_i.mua_table);
    end
end

fprintf('\nFinished missing completed-model postprocessing.\n');
fprintf('Window rows: %d\n', n_window_rows);
fprintf('SPKT rows: %d\n', n_spkt_rows);
fprintf('MUA rows: %d\n', n_mua_rows);
