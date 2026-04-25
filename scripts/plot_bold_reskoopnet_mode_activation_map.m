%% Plot a BOLD ResKoopNet Koopman-mode activation map for one eigenbasis
% Edit the user settings below, then run this file in MATLAB.

clear; close all; clc;

%% -------------------- user settings --------------------
result_mat = ['E:\autodl_results_local\bold_wsl\e10gb1\mlp\outputs\', ...
    'mlp_obs_bold_wsl_20260423_gpu_e10gb1_projected_vlambda_HP\', ...
    'E10.gb1_bold_observables_HP_Python_resdmd_Layer_100_Ndict_1205_outputs_1.mat'];

roi_ts_mat = 'E:\DataPons\E10.gb1\roits\e10gb1_0001_roits.mat';

roi_index = 29;          % E10.gb1 roiTs{29} = HP, matching the HP observable.
basis_index = 2;         % Eigenbasis/mode index to plot.
index_mode = 'sorted';   % 'sorted': basis_index after sorting by |lambda| desc; 'raw': file order.
mode_values = 'abs';     % 'abs' or 'real'. Use 'real' for signed maps.

slice_list = [];         % [] = all slices containing ROI voxels.
marker_size = 18;
marker_alpha = 0.85;
save_figures = true;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
out_dir = fullfile(project_root, 'analysis_outputs', 'bold_reskoopnet_activation_maps');

%% -------------------- load outputs and ROI template --------------------
S = load(result_mat, 'EDMD_outputs');
EDMD_outputs = S.EDMD_outputs;

R = load(roi_ts_mat, 'roiTs');
roi = R.roiTs{1, roi_index};

if ~isfield(EDMD_outputs, 'kpm_modes')
    error('EDMD_outputs.kpm_modes is missing in %s.', result_mat);
end
if ~isfield(EDMD_outputs, 'evalues')
    error('EDMD_outputs.evalues is missing in %s.', result_mat);
end
if ~isfield(EDMD_outputs, 'efuns')
    warning('EDMD_outputs.efuns is missing. The activation map can still be plotted.');
end

evalues = EDMD_outputs.evalues(:);
kpm_modes = EDMD_outputs.kpm_modes;
if size(kpm_modes, 1) ~= numel(evalues)
    error('Dimension mismatch: kpm_modes has %d rows but evalues has %d entries.', ...
        size(kpm_modes, 1), numel(evalues));
end

coords = roi.coords;
ana = roi.ana;
if size(coords, 2) ~= 3
    error('roi.coords must be [N_voxel x 3].');
end

if size(kpm_modes, 2) ~= size(coords, 1)
    error(['Mode dimension (%d) does not match ROI coordinate count (%d). ', ...
        'Use a direct spatial observable such as HP/eleHP, or back-project SVD modes first.'], ...
        size(kpm_modes, 2), size(coords, 1));
end

%% -------------------- choose eigenbasis --------------------
switch lower(index_mode)
    case 'sorted'
        [~, order] = sort(abs(evalues), 'descend');
        raw_index = order(basis_index);
    case 'raw'
        raw_index = basis_index;
    otherwise
        error('index_mode must be ''sorted'' or ''raw''.');
end

if raw_index < 1 || raw_index > size(kpm_modes, 1)
    error('Selected raw_index=%d is outside [1, %d].', raw_index, size(kpm_modes, 1));
end

lambda = evalues(raw_index);
mode_vec = kpm_modes(raw_index, :);

switch lower(mode_values)
    case 'abs'
        vals = abs(mode_vec(:));
        vals = local_minmax(vals);
        value_label = '|mode weight|';
        cmap = turbo(256);
        clim_use = [0 1];
    case 'real'
        vals = real(mode_vec(:));
        max_abs = max(abs(vals), [], 'omitnan');
        if max_abs > 0
            vals = vals ./ max_abs;
        end
        value_label = 'real(mode weight)';
        cmap = local_redblue(256);
        clim_use = [-1 1];
    otherwise
        error('mode_values must be ''abs'' or ''real''.');
end

if isempty(slice_list)
    slice_list = unique(coords(:, 3))';
    slice_list = slice_list(slice_list >= 1 & slice_list <= size(ana, 3));
end

if isempty(slice_list)
    error('No valid slices found in roi.coords for anatomy with %d slices.', size(ana, 3));
end

%% -------------------- plot slice mosaic --------------------
fig_mosaic = figure('Color', 'w', 'Name', 'BOLD ResKoopNet mode activation map');
n_slice = numel(slice_list);
n_col = min(8, n_slice);
n_row = ceil(n_slice / n_col);
tiledlayout(n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

for ii = 1:n_slice
    z = slice_list(ii);
    nexttile;

    img = squeeze(ana(:, :, z))';
    imagesc(img);
    colormap(gca, gray);
    axis image xy off;
    hold on;

    idx = coords(:, 3) == z;
    if any(idx)
        scatter(coords(idx, 1), coords(idx, 2), marker_size, vals(idx), ...
            'filled', 'MarkerEdgeColor', 'none', ...
            'MarkerFaceAlpha', marker_alpha);
        colormap(gca, cmap);
        clim(clim_use);
    end
    title(sprintf('S%d', z), 'FontWeight', 'normal');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = value_label;

sgtitle(sprintf('%s | %s index %d -> raw %d | |lambda|=%.4g', ...
    local_basename(result_mat), index_mode, basis_index, raw_index, abs(lambda)), ...
    'Interpreter', 'none');

%% -------------------- plot 3D scatter summary --------------------
fig_3d = figure('Color', 'w', 'Name', 'BOLD ResKoopNet mode activation map 3D');
scatter3(coords(:, 1), coords(:, 2), coords(:, 3), marker_size, vals, ...
    'filled', 'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', marker_alpha);
axis equal; grid on;
xlabel('x'); ylabel('y'); zlabel('slice');
title(sprintf('ROI %s | %s index %d -> raw %d | |lambda|=%.4g', ...
    roi.name, index_mode, basis_index, raw_index, abs(lambda)), 'Interpreter', 'none');
colormap(cmap);
clim(clim_use);
cb = colorbar;
cb.Label.String = value_label;
view(3);

%% -------------------- optional eigenfunction trace --------------------
fig_trace = [];
if isfield(EDMD_outputs, 'efuns') && size(EDMD_outputs.efuns, 2) >= raw_index
    fig_trace = figure('Color', 'w', 'Name', 'Selected eigenfunction trace');
    efun = EDMD_outputs.efuns(:, raw_index);
    plot(abs(efun), 'k', 'LineWidth', 1);
    xlabel('Time index');
    ylabel('|eigenfunction|');
    title(sprintf('Eigenfunction trace | %s index %d -> raw %d', ...
        index_mode, basis_index, raw_index), 'Interpreter', 'none');
    box off;
end

%% -------------------- save --------------------
if save_figures
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    tag = sprintf('%s_%sIdx%03d_raw%04d_%s', roi.name, index_mode, basis_index, raw_index, mode_values);
    exportgraphics(fig_mosaic, fullfile(out_dir, ['activation_mosaic_' tag '.png']), 'Resolution', 200);
    savefig(fig_mosaic, fullfile(out_dir, ['activation_mosaic_' tag '.fig']));
    exportgraphics(fig_3d, fullfile(out_dir, ['activation_3d_' tag '.png']), 'Resolution', 200);
    savefig(fig_3d, fullfile(out_dir, ['activation_3d_' tag '.fig']));
    if ~isempty(fig_trace)
        exportgraphics(fig_trace, fullfile(out_dir, ['eigenfunction_trace_' tag '.png']), 'Resolution', 200);
        savefig(fig_trace, fullfile(out_dir, ['eigenfunction_trace_' tag '.fig']));
    end
    fprintf('Saved activation figures to:\n  %s\n', out_dir);
end

fprintf('Selected %s basis index %d (raw index %d): lambda = %.6g%+.6gi, |lambda| = %.6g\n', ...
    index_mode, basis_index, raw_index, real(lambda), imag(lambda), abs(lambda));

%% -------------------- local helpers --------------------
function y = local_minmax(x)
    x = double(x);
    mn = min(x, [], 'omitnan');
    mx = max(x, [], 'omitnan');
    if ~isfinite(mn) || ~isfinite(mx) || mx <= mn
        y = zeros(size(x));
    else
        y = (x - mn) ./ (mx - mn);
    end
end

function cmap = local_redblue(n)
    if nargin < 1 || isempty(n)
        n = 256;
    end
    x = linspace(-1, 1, n)';
    cmap = zeros(n, 3);
    cmap(:, 1) = max(0, x);
    cmap(:, 3) = max(0, -x);
    cmap(:, 2) = 1 - abs(x);
    cmap = 0.15 + 0.85 * cmap;
    cmap = min(max(cmap, 0), 1);
end

function name = local_basename(path_str)
    [~, name, ext] = fileparts(path_str);
    name = [name ext];
end
