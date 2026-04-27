%% Reference-style BOLD ResKoopNet Koopman-mode activation map
% This version makes a slice mosaic like the manuscript/reference panels:
% anatomical slice in gray/black, activation values overlaid as a colored image.

clear; close all; clc;

%% -------------------- user settings --------------------
result_mat = ['E:\autodl_results_local\bold_wsl\e10gb1\mlp\outputs\', ...
    'mlp_obs_bold_wsl_20260423_cpu2_e10gb1_projected_vlambda_global_svd100\', ...
    'E10.gb1_bold_observables_global_svd100_Python_resdmd_Layer_100_Ndict_201_outputs_1.mat'];

observable_mat = 'E:\DataPons_processed\bold_observables\E10.gb1\E10.gb1_bold_observables_global_svd100.mat';
roi_ts_mat = 'E:\DataPons\E10.gb1\roits\e10gb1_0001_roits.mat';

basis_index = 2;         % eigenbasis / Koopman mode index to plot
index_mode = 'sorted';   % 'sorted' = sorted by |lambda| desc; 'raw' = file order
value_mode = 'abs';      % 'abs' for 0-1 activation; 'real' for signed mode

slice_list = 1:20;       % reference figure style: S1-S20
tiles_per_row = 10;
overlay_alpha = 0.86;
background_limits = [];  % [] = auto robust grayscale limits

save_figures = true;
script_dir = fileparts(mfilename('fullpath'));
project_root = fileparts(script_dir);
out_dir = fullfile(project_root, 'analysis_outputs', 'bold_reskoopnet_activation_maps_reference_style');

%% -------------------- load ResKoopNet and observable model --------------------
S = load(result_mat, 'EDMD_outputs');
EDMD_outputs = S.EDMD_outputs;

Obs = load(observable_mat, 'O', 'params');
O = Obs.O;

if ~isfield(EDMD_outputs, 'kpm_modes')
    error('EDMD_outputs.kpm_modes is missing in %s.', result_mat);
end
if ~isfield(EDMD_outputs, 'evalues')
    error('EDMD_outputs.evalues is missing in %s.', result_mat);
end

evalues = EDMD_outputs.evalues(:);
kpm_modes = EDMD_outputs.kpm_modes;

if size(kpm_modes, 1) ~= numel(evalues)
    error('Dimension mismatch: kpm_modes has %d rows but evalues has %d entries.', ...
        size(kpm_modes, 1), numel(evalues));
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
mode_obs = double(kpm_modes(raw_index, :));

%% -------------------- project mode into voxel space --------------------
% Direct voxel observables have kpm_modes rows already in voxel space.
% SVD observables have kpm_modes rows in score space; use O.model.coeff to
% map score-space mode weights back to original voxel variables.
if isfield(O, 'model') && isfield(O.model, 'coeff') && size(O.model.coeff, 2) == numel(mode_obs)
    voxel_values = mode_obs * double(O.model.coeff).';
    map_space = sprintf('%s back-projected to voxels', char(string(Obs.params.observable_branch)));
else
    voxel_values = mode_obs;
    map_space = 'direct voxel mode';
end

%% -------------------- load first-session ROI coordinates and anatomy --------------------
R = load(roi_ts_mat, 'roiTs');
[coords, region_labels, n_var_region, ana] = local_collect_roi_coords_and_anatomy(R.roiTs);

if numel(voxel_values) ~= size(coords, 1)
    error(['Voxel-space mode has %d values but first-session ROI coordinates have %d voxels. ', ...
        'Check that result_mat and observable_mat describe the same observable branch/dataset.'], ...
        numel(voxel_values), size(coords, 1));
end

vals = local_prepare_values(voxel_values(:), value_mode);
[cmap, clim_use, value_label] = local_colormap_and_limits(value_mode);
act_vol = local_values_to_volume(vals, coords, size(ana));

if isempty(background_limits)
    background_limits = local_robust_limits(double(ana(:)));
end

%% -------------------- plot reference-style slice mosaic --------------------
n_slice = numel(slice_list);
n_col = tiles_per_row;
n_row = ceil(n_slice / n_col);

fig = figure('Color', 'w', 'Name', 'Reference-style BOLD ResKoopNet activation map');
fig.Position = [80 80 1450 170 * n_row + 110];
tiledlayout(n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

for ii = 1:n_slice
    z = slice_list(ii);
    nexttile;

    if z < 1 || z > size(ana, 3)
        axis off;
        title(sprintf('S%d', z));
        continue;
    end

    img = squeeze(ana(:, :, z))';
    overlay = squeeze(act_vol(:, :, z))';
    alpha_data = overlay_alpha * double(~isnan(overlay));

    image(local_gray_rgb(img, background_limits));
    axis image xy off;
    set(gca, 'Color', 'k');
    hold on;

    h = imagesc(overlay);
    set(h, 'AlphaData', alpha_data);
    colormap(gca, cmap);
    clim(clim_use);

    title(sprintf('S%d', z), 'FontSize', 9, 'FontWeight', 'bold');
end

cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = value_label;

sgtitle(sprintf('%s | %s index %d -> raw %d | |lambda|=%.4g | %s', ...
    local_basename(result_mat), index_mode, basis_index, raw_index, abs(lambda), map_space), ...
    'Interpreter', 'none');

%% -------------------- save --------------------
if save_figures
    if ~exist(out_dir, 'dir')
        mkdir(out_dir);
    end
    tag = sprintf('%sIdx%03d_raw%04d_%s', index_mode, basis_index, raw_index, value_mode);
    png_file = fullfile(out_dir, ['activation_reference_mosaic_' tag '.png']);
    fig_file = fullfile(out_dir, ['activation_reference_mosaic_' tag '.fig']);
    exportgraphics(fig, png_file, 'Resolution', 220);
    savefig(fig, fig_file);
    fprintf('Saved reference-style activation map:\n  %s\n', png_file);
end

fprintf('Selected %s basis index %d (raw index %d): lambda = %.6g%+.6gi, |lambda| = %.6g\n', ...
    index_mode, basis_index, raw_index, real(lambda), imag(lambda), abs(lambda));
fprintf('Collected %d voxel coordinates from %d ROI(s): %s\n', ...
    size(coords, 1), numel(region_labels), strjoin(region_labels(:).', ', '));
fprintf('n_var_region first entries: %s\n', mat2str(n_var_region(1:min(10,end)).'));

%% -------------------- local helpers --------------------
function [coords, region_labels, n_var_region, ana] = local_collect_roi_coords_and_anatomy(roiTs)
    if ~iscell(roiTs)
        error('Expected roiTs to be a cell array of ROI structs.');
    end

    coords_cells = {};
    region_labels = {};
    n_var_region = [];
    ana = [];

    for i_region = 1:size(roiTs, 2)
        R = roiTs{1, i_region};
        if ~isstruct(R) || ~isfield(R, 'coords') || isempty(R.coords)
            continue;
        end
        if ~isfield(R, 'dat') || size(R.dat, 2) ~= size(R.coords, 1)
            continue;
        end
        coords_cells{end+1, 1} = double(R.coords); %#ok<AGROW>
        n_var_region(end+1, 1) = size(R.coords, 1); %#ok<AGROW>
        region_labels{end+1, 1} = local_roi_label(R, i_region); %#ok<AGROW>
        if isempty(ana) && isfield(R, 'ana') && ~isempty(R.ana)
            ana = R.ana;
        end
    end

    if isempty(coords_cells)
        error('No ROI coordinates were found in roiTs.');
    end
    if isempty(ana)
        error('No anatomy field .ana was found in roiTs.');
    end
    coords = cat(1, coords_cells{:});
end

function label = local_roi_label(R, i_region)
    if ~isfield(R, 'name') || isempty(R.name)
        error('roiTs{1,%d} is missing required ROI label field .name.', i_region);
    end
    if ischar(R.name) || isstring(R.name)
        label = char(string(R.name));
        return;
    end
    if iscell(R.name) && numel(R.name) == 1 && (ischar(R.name{1}) || isstring(R.name{1}))
        label = char(string(R.name{1}));
        return;
    end
    error('roiTs{1,%d}.name must be a string-like scalar.', i_region);
end

function vals = local_prepare_values(raw_vals, value_mode)
    raw_vals = double(raw_vals);
    switch lower(value_mode)
        case 'abs'
            vals = abs(raw_vals);
            vals = local_minmax(vals);
        case 'real'
            vals = real(raw_vals);
            max_abs = max(abs(vals), [], 'omitnan');
            if isfinite(max_abs) && max_abs > 0
                vals = vals ./ max_abs;
            end
        otherwise
            error('value_mode must be ''abs'' or ''real''.');
    end
end

function vol = local_values_to_volume(vals, coords, vol_size)
    coords = round(double(coords));
    in_bounds = coords(:, 1) >= 1 & coords(:, 1) <= vol_size(1) & ...
        coords(:, 2) >= 1 & coords(:, 2) <= vol_size(2) & ...
        coords(:, 3) >= 1 & coords(:, 3) <= vol_size(3);
    if ~all(in_bounds)
        warning('Dropping %d out-of-bounds coordinate(s).', sum(~in_bounds));
    end

    vol = nan(vol_size);
    idx = sub2ind(vol_size, coords(in_bounds, 1), coords(in_bounds, 2), coords(in_bounds, 3));
    vol(idx) = vals(in_bounds);
end

function [cmap, clim_use, value_label] = local_colormap_and_limits(value_mode)
    switch lower(value_mode)
        case 'abs'
            cmap = turbo(256);
            clim_use = [0 1];
            value_label = 'normalized |mode weight|';
        case 'real'
            cmap = local_redblue(256);
            clim_use = [-1 1];
            value_label = 'normalized real(mode weight)';
    end
end

function lim = local_robust_limits(x)
    x = x(isfinite(x));
    if isempty(x)
        lim = [0 1];
        return;
    end
    lo = prctile(x, 1);
    hi = prctile(x, 99.5);
    if hi <= lo
        hi = max(x);
        lo = min(x);
    end
    if hi <= lo
        lim = [lo lo + 1];
    else
        lim = [lo hi];
    end
end

function rgb = local_gray_rgb(img, lim)
    img = double(img);
    if nargin < 2 || isempty(lim)
        lim = local_robust_limits(img(:));
    end
    denom = lim(2) - lim(1);
    if denom <= 0 || ~isfinite(denom)
        denom = 1;
    end
    g = (img - lim(1)) ./ denom;
    g = min(max(g, 0), 1);
    rgb = repmat(g, 1, 1, 3);
end

function y = local_minmax(x)
    mn = min(x, [], 'omitnan');
    mx = max(x, [], 'omitnan');
    if ~isfinite(mn) || ~isfinite(mx) || mx <= mn
        y = zeros(size(x));
    else
        y = (x - mn) ./ (mx - mn);
    end
end

function cmap = local_redblue(n)
    x = linspace(-1, 1, n)';
    cmap = zeros(n, 3);
    cmap(:, 1) = max(0, x);
    cmap(:, 3) = max(0, -x);
    cmap(:, 2) = 1 - abs(x);
    cmap = 0.12 + 0.88 * cmap;
    cmap = min(max(cmap, 0), 1);
end

function name = local_basename(path_str)
    [~, name, ext] = fileparts(path_str);
    name = [name ext];
end
