% Maintenance batch runner for completed E10 BLP eigenfunction postprocessing.
%
% Windows are always taken from each dataset's saved consensus-state-diversity
% top-30 global windows. Outputs are written under each dataset's
% pipeline6_* folders. F12m01 is intentionally excluded.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_names = {'E10fV1', 'E10gb1', 'E10gH1'};
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

fprintf('\nFinished completed MLP top-window postprocessing.\n');
fprintf('Window rows: %d\n', n_window_rows);
fprintf('SPKT rows: %d\n', n_spkt_rows);
fprintf('MUA rows: %d\n', n_mua_rows);
