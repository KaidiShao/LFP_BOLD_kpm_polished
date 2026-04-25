% Batch-generate missing BOLD activation maps from existing BOLD post + xcorr.
%
% This script only processes runs that already have both:
%   1) E:\DataPons_processed\<dataset>\bold_reskoopnet_postprocessing\<run>\mat\*_bold_post.mat
%   2) E:\DataPons_processed\<dataset>\bold_reskoopnet_postprocessing\<run>\density_cross_correlation\bold_efun_density_xcorr.mat
%
% Outputs:
%   E:\DataPons_processed\<dataset>\bold_reskoopnet_single_run_analysis\<run>\fig\activation_maps_top5

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');
try
    opengl('software');
catch
end

if ~exist('params', 'var') || ~isstruct(params)
    params = struct();
end
params = local_apply_defaults(params);

manifest = local_generate_activation_maps(params);

fprintf('\nFinished BOLD activation-map batch.\n');
fprintf('Manifest CSV:\n  %s\n', manifest.csv_file);


function manifest = local_generate_activation_maps(params)
candidates = local_discover_candidates(params);
if isempty(candidates)
    error('No BOLD post + xcorr pairs were found for the requested filters.');
end

rows = repmat(local_empty_row(), numel(candidates), 1);
row_count = 0;

fprintf('Discovered %d activation-map candidate run(s).\n', numel(candidates));

for i_run = 1:numel(candidates)
    cand = candidates(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s\n', i_run, numel(candidates), ...
        cand.dataset_stem, cand.run_name);
    fprintf('BOLD_POST:\n  %s\n', cand.bold_post_file);
    fprintf('XCORR:\n  %s\n', cand.xcorr_file);

    t_run = tic;
    try
        S_xcorr = load(cand.xcorr_file, 'out');
        if ~isfield(S_xcorr, 'out')
            error('XCORR variable ''out'' missing in %s.', cand.xcorr_file);
        end
        xcorr_out = S_xcorr.out;
        if ~isfield(xcorr_out, 'top_table') || isempty(xcorr_out.top_table)
            row_count = row_count + 1;
            rows(row_count) = local_make_row(cand, 'no_top_table', ...
                'XCORR top_table is empty.', '', 0, 0, toc(t_run));
            fprintf('XCORR top_table is empty; skipping.\n');
            continue;
        end

        result_dir = local_resolve_result_dir(cand, params);
        raw_outputs = local_load_first_edmd_output(result_dir);
        observable_file = local_resolve_observable_file(cand, params);
        Obs = local_load_observable_struct(observable_file);
        roi_ts_file = local_resolve_roi_ts_file(cand, params);
        plot_ctx = local_prepare_plot_context(raw_outputs, Obs, roi_ts_file, ...
            params.value_mode);

        out_root = fullfile(params.processed_root, cand.dataset_stem, ...
            params.output_folder_name, cand.run_name);
        mat_dir = fullfile(out_root, 'mat');
        fig_dir = fullfile(out_root, 'fig');
        act_dir = fullfile(fig_dir, 'activation_maps_top5');
        if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
        if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end
        if exist(act_dir, 'dir') ~= 7, mkdir(act_dir); end

        raw_indices = local_top_raw_indices(xcorr_out.top_table, ...
            raw_outputs, params.top_n);
        raw_indices = unique(raw_indices(:).', 'stable');
        raw_indices = raw_indices(1:min(numel(raw_indices), params.top_n));
        if isempty(raw_indices)
            row_count = row_count + 1;
            rows(row_count) = local_make_row(cand, 'no_raw_indices', ...
                'XCORR top_table did not map to any raw Koopman basis index.', ...
                act_dir, 0, 0, toc(t_run));
            fprintf('No raw Koopman basis indices were resolved; skipping.\n');
            continue;
        end

        n_saved = 0;
        n_existing = 0;
        for i_idx = 1:numel(raw_indices)
            raw_idx = raw_indices(i_idx);
            png_file = fullfile(act_dir, sprintf( ...
                'top%02d_raw%04d_%s_activation_reference.png', ...
                i_idx, raw_idx, params.value_mode));
            fig_file = fullfile(act_dir, sprintf( ...
                'top%02d_raw%04d_%s_activation_reference.fig', ...
                i_idx, raw_idx, params.value_mode));

            if params.skip_existing && exist(png_file, 'file') == 2 && ...
                    (~params.save_fig || exist(fig_file, 'file') == 2)
                n_existing = n_existing + 1;
                fprintf('Existing activation map found; skipping raw index %d.\n', raw_idx);
                continue;
            end

            meta = local_top_meta_for_raw_index( ...
                xcorr_out.top_table, raw_outputs, raw_idx);
            fig = local_plot_reference_activation_map_ctx( ...
                plot_ctx, raw_idx, params.slice_list, params.tiles_per_row, meta);

            if params.save_png
                exportgraphics(fig, png_file, 'Resolution', 220);
            end
            if params.save_fig
                savefig(fig, fig_file);
            end
            close(fig);
            drawnow limitrate;
            n_saved = n_saved + 1;
        end

        activation_info = struct();
        activation_info.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
        activation_info.dataset_stem = cand.dataset_stem;
        activation_info.run_name = cand.run_name;
        activation_info.observable_mode = cand.observable_mode;
        activation_info.residual_form = cand.residual_form;
        activation_info.bold_post_file = cand.bold_post_file;
        activation_info.xcorr_file = cand.xcorr_file;
        activation_info.observable_file = observable_file;
        activation_info.roi_ts_file = roi_ts_file;
        activation_info.result_dir = result_dir;
        activation_info.value_mode = params.value_mode;
        activation_info.slice_list = params.slice_list;
        activation_info.tiles_per_row = params.tiles_per_row;
        activation_info.raw_indices = raw_indices;
        activation_info.top_table = xcorr_out.top_table( ...
            1:min(height(xcorr_out.top_table), params.top_n), :);
        save(fullfile(mat_dir, 'activation_map_batch_info.mat'), ...
            'activation_info', '-v7.3');

        if n_saved > 0
            status = 'ok';
            message = sprintf('Saved %d map(s); %d already existed.', ...
                n_saved, n_existing);
        else
            status = 'skipped_existing';
            message = sprintf('All %d requested map(s) already existed.', ...
                numel(raw_indices));
        end

        row_count = row_count + 1;
        rows(row_count) = local_make_row(cand, status, message, act_dir, ...
            numel(raw_indices), n_saved, toc(t_run));
        clear S_xcorr xcorr_out raw_outputs Obs plot_ctx
    catch ME
        message = getReport(ME, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s failed:\n%s\n', cand.run_name, message);
        row_count = row_count + 1;
        rows(row_count) = local_make_row(cand, 'error', message, '', ...
            0, 0, toc(t_run));
        if ~params.continue_on_error
            rethrow(ME);
        end
    end
    fclose('all');
    close all force;
end

rows = rows(1:row_count);
T = struct2table(rows);
tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest_dir = fullfile(params.processed_root, 'postprocessing_manifests', ...
    'bold_activation_maps');
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end

manifest = struct();
manifest.params = params;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(manifest_dir, ...
    ['bold_activation_maps_manifest_', tag, '.csv']);
manifest.mat_file = fullfile(manifest_dir, ...
    ['bold_activation_maps_manifest_', tag, '.mat']);

writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');
end


function candidates = local_discover_candidates(params)
dataset_stems = cellstr(string(params.dataset_stems(:)).');
observable_modes = cellstr(string(params.observable_modes(:)).');
residual_forms = cellstr(string(params.residual_forms(:)).');
run_name_contains = cellstr(string(params.run_name_contains(:)).');

candidates = repmat(local_empty_candidate(), 0, 1);
dataset_dirs = dir(params.processed_root);
dataset_dirs = dataset_dirs([dataset_dirs.isdir]);
dataset_dirs = dataset_dirs(~ismember({dataset_dirs.name}, {'.', '..'}));

for i_ds = 1:numel(dataset_dirs)
    dataset_stem = dataset_dirs(i_ds).name;
    if ~isempty(dataset_stems) && ...
            ~any(strcmpi(dataset_stem, dataset_stems))
        continue;
    end

    L = dir(fullfile(params.processed_root, dataset_stem, ...
        params.post_folder_name, '*', 'mat', '*_bold_post.mat'));
    for i_file = 1:numel(L)
        bold_post_file = fullfile(L(i_file).folder, L(i_file).name);
        [~, base_name] = fileparts(bold_post_file);
        run_name = char(erase(string(base_name), "_bold_post"));
        observable_mode = local_parse_observable_mode(run_name);
        residual_form = local_parse_residual_form(run_name);
        if ~isempty(observable_modes) && ...
                ~any(strcmpi(observable_mode, observable_modes))
            continue;
        end
        if ~isempty(residual_forms) && ...
                ~any(strcmpi(residual_form, residual_forms))
            continue;
        end
        if ~local_matches_all(run_name, run_name_contains)
            continue;
        end

        run_root = fileparts(fileparts(bold_post_file));
        xcorr_file = fullfile(run_root, params.xcorr_folder_name, ...
            params.xcorr_file_name);
        if exist(xcorr_file, 'file') ~= 2
            continue;
        end

        cand = local_empty_candidate();
        cand.dataset_stem = dataset_stem;
        cand.run_name = run_name;
        cand.observable_mode = observable_mode;
        cand.residual_form = residual_form;
        cand.bold_post_file = bold_post_file;
        cand.xcorr_file = xcorr_file;
        candidates(end + 1, 1) = cand; %#ok<AGROW>
    end
end

if isempty(candidates)
    return;
end

[~, order] = sortrows([{candidates.dataset_stem}.', {candidates.run_name}.'], [1 2]);
candidates = candidates(order);
end


function tf = local_matches_all(text_value, patterns)
tf = true;
for i = 1:numel(patterns)
    if ~contains(text_value, patterns{i}, 'IgnoreCase', true)
        tf = false;
        return;
    end
end
end


function params = local_apply_defaults(params)
params = local_set_default(params, 'processed_root', get_project_processed_root());
params = local_set_default(params, 'datapons_root', 'E:\DataPons');
params = local_set_default(params, 'autodl_roots', ...
    {'E:\autodl_results_local\bold_wsl', 'E:\autodl_results\bold'});

params = local_set_default(params, 'dataset_stems', {});
params = local_set_default(params, 'observable_modes', {});
params = local_set_default(params, 'residual_forms', {});
params = local_set_default(params, 'run_name_contains', {});

params = local_set_default(params, 'post_folder_name', ...
    'bold_reskoopnet_postprocessing');
params = local_set_default(params, 'output_folder_name', ...
    'bold_reskoopnet_single_run_analysis');
params = local_set_default(params, 'xcorr_folder_name', ...
    'density_cross_correlation');
params = local_set_default(params, 'xcorr_file_name', ...
    'bold_efun_density_xcorr.mat');

params = local_set_default(params, 'top_n', 5);
params = local_set_default(params, 'value_mode', 'abs');
params = local_set_default(params, 'slice_list', 1:20);
params = local_set_default(params, 'tiles_per_row', 10);
params = local_set_default(params, 'skip_existing', true);
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'continue_on_error', true);
end


function result_dir = local_resolve_result_dir(cand, params)
result_dir = '';

for i_root = 1:numel(params.autodl_roots)
    candidate = fullfile(params.autodl_roots{i_root}, cand.dataset_stem, ...
        'mlp', 'outputs', cand.run_name);
    if exist(candidate, 'dir') == 7
        result_dir = candidate;
        return;
    end
end

error('Could not resolve result_dir for %s.', cand.run_name);
end


function observable_file = local_resolve_observable_file(cand, params)
observable_file = '';
dataset_dir = local_dataset_dir_name(cand.dataset_stem);
candidate = fullfile(params.processed_root, 'bold_observables', dataset_dir, ...
    sprintf('%s_bold_observables_%s.mat', dataset_dir, cand.observable_mode));
if exist(candidate, 'file') == 2
    observable_file = candidate;
    return;
end

error('Could not resolve observable MAT file for %s.', cand.run_name);
end


function Obs = local_load_observable_struct(observable_file)
available = whos('-file', observable_file);
available_names = {available.name};
load_names = intersect({'O', 'cfg', 'params', 'obs_info', ...
    'observable_info', 'observable_labels', 'dataset_id', 'file_stem'}, ...
    available_names, 'stable');
Obs = load(observable_file, load_names{:});
if ~isfield(Obs, 'O')
    Obs.O = struct();
end
if ~isfield(Obs, 'cfg')
    Obs.cfg = struct();
end
if ~isfield(Obs, 'params')
    Obs.params = struct();
end
end


function roi_ts_file = local_resolve_roi_ts_file(cand, params)
dataset_dir = local_dataset_dir_name(cand.dataset_stem);
roi_dir = fullfile(params.datapons_root, dataset_dir, 'roits');
if exist(roi_dir, 'dir') ~= 7
    error('ROI directory not found: %s', roi_dir);
end

L = dir(fullfile(roi_dir, sprintf('%s_*_roits.mat', lower(cand.dataset_stem))));
if isempty(L)
    L = dir(fullfile(roi_dir, '*_roits.mat'));
end
if isempty(L)
    error('No roiTs MAT files found in %s.', roi_dir);
end
names = {L.name};
[~, order] = sort(names);
L = L(order);
roi_ts_file = fullfile(L(1).folder, L(1).name);
end


function raw_outputs = local_load_first_edmd_output(result_dir)
L = dir(fullfile(result_dir, '*_outputs_1.mat'));
if isempty(L)
    L = dir(fullfile(result_dir, '*_outputs_*.mat'));
end
if isempty(L)
    error('No EDMD output chunks found in %s.', result_dir);
end
[~, order] = sort({L.name});
L = L(order);
S = load(fullfile(L(1).folder, L(1).name), 'EDMD_outputs');
if ~isfield(S, 'EDMD_outputs')
    error('EDMD_outputs missing in %s.', fullfile(L(1).folder, L(1).name));
end
raw_outputs = S.EDMD_outputs;
end


function raw_indices = local_top_raw_indices(top_table, EDMD_outputs, top_n)
n = min(height(top_table), top_n);
raw_indices = nan(n, 1);
for ii = 1:n
    mode_idx = top_table.bold_mode_index(ii);
    if isfield(EDMD_outputs, 'idx_final_in_original') && ...
            numel(EDMD_outputs.idx_final_in_original) >= mode_idx
        raw_indices(ii) = EDMD_outputs.idx_final_in_original(mode_idx);
    else
        raw_indices(ii) = mode_idx;
    end
end
raw_indices = raw_indices(isfinite(raw_indices));
raw_indices = round(raw_indices(:));
raw_indices = raw_indices(raw_indices >= 1);
end


function meta = local_top_meta_for_raw_index(top_table, EDMD_outputs, raw_idx)
meta = struct('text', sprintf('raw index %d', raw_idx));
for ii = 1:height(top_table)
    mode_idx = top_table.bold_mode_index(ii);
    if isfield(EDMD_outputs, 'idx_final_in_original') && ...
            numel(EDMD_outputs.idx_final_in_original) >= mode_idx
        row_raw = EDMD_outputs.idx_final_in_original(mode_idx);
    else
        row_raw = mode_idx;
    end
    if row_raw == raw_idx
        meta.text = sprintf('%s | %s | mode %d raw %d | corr %.3f lag %.1fs', ...
            char(top_table.density_name(ii)), char(top_table.bold_feature(ii)), ...
            mode_idx, raw_idx, top_table.peak_corr(ii), top_table.peak_lag_sec(ii));
        return;
    end
end
end


function plot_ctx = local_prepare_plot_context(EDMD_outputs, Obs, roi_ts_mat, value_mode)
plot_ctx = struct();
plot_ctx.evalues = EDMD_outputs.evalues(:);
plot_ctx.kpm_modes = EDMD_outputs.kpm_modes;
plot_ctx.value_mode = value_mode;

O = Obs.O;
plot_ctx.direct_voxel_mode = true;
plot_ctx.coeff_t = [];
if isfield(O, 'model') && isfield(O.model, 'coeff') && ...
        size(O.model.coeff, 2) == size(plot_ctx.kpm_modes, 2)
    plot_ctx.direct_voxel_mode = false;
    plot_ctx.coeff_t = double(O.model.coeff).';
end

R = load(roi_ts_mat, 'roiTs');
if ~isfield(R, 'roiTs')
    error('roiTs variable missing in %s.', roi_ts_mat);
end
target_region_names = local_target_region_names_from_observable(Obs);
[plot_ctx.coords, plot_ctx.region_labels, ~, plot_ctx.ana] = ...
    local_collect_roi_coords_and_anatomy(R.roiTs, target_region_names);
plot_ctx.bg_limits = local_robust_limits(double(plot_ctx.ana(:)));
[plot_ctx.cmap, plot_ctx.clim_use, plot_ctx.value_label] = ...
    local_colormap_and_limits(value_mode);
end


function fig = local_plot_reference_activation_map_ctx(plot_ctx, ...
        raw_index, slice_list, tiles_per_row, meta)
kpm_modes = plot_ctx.kpm_modes;
if raw_index < 1 || raw_index > size(kpm_modes, 1)
    error('raw_index=%d is outside [1, %d].', raw_index, size(kpm_modes, 1));
end
lambda = plot_ctx.evalues(raw_index);
mode_obs = double(kpm_modes(raw_index, :));

if ~plot_ctx.direct_voxel_mode
    voxel_values = mode_obs * plot_ctx.coeff_t;
    map_space = 'SVD back-projected to voxels';
else
    voxel_values = mode_obs;
    map_space = 'direct voxel mode';
end

if numel(voxel_values) ~= size(plot_ctx.coords, 1)
    error('Voxel-space mode has %d values but roiTs has %d voxel coordinates.', ...
        numel(voxel_values), size(plot_ctx.coords, 1));
end

vals = local_prepare_values(voxel_values(:), plot_ctx.value_mode);
act_vol = local_values_to_volume(vals, plot_ctx.coords, size(plot_ctx.ana));

n_slice = numel(slice_list);
n_col = tiles_per_row;
n_row = ceil(n_slice / n_col);
fig = figure('Color', 'w', 'Name', 'Reference-style BOLD activation map');
fig.Position = [80 80 1450 170 * n_row + 130];
tiledlayout(n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

for jj = 1:n_slice
    z = slice_list(jj);
    nexttile;
    if z < 1 || z > size(plot_ctx.ana, 3)
        axis off;
        title(sprintf('S%d', z));
        continue;
    end
    img = squeeze(plot_ctx.ana(:, :, z))';
    overlay = squeeze(act_vol(:, :, z))';
    alpha_data = 0.86 * double(~isnan(overlay));

    image(local_gray_rgb(img, plot_ctx.bg_limits));
    axis image xy off;
    set(gca, 'Color', 'k');
    hold on;
    h = imagesc(overlay);
    set(h, 'AlphaData', alpha_data);
    colormap(gca, plot_ctx.cmap);
    clim(plot_ctx.clim_use);
    title(sprintf('S%d', z), 'FontSize', 9, 'FontWeight', 'bold');
end
cb = colorbar;
cb.Layout.Tile = 'east';
cb.Label.String = plot_ctx.value_label;
sgtitle(sprintf('%s | |lambda|=%.4g | %s | %s | ROI(s): %s', ...
    meta.text, abs(lambda), plot_ctx.value_mode, map_space, ...
    strjoin(plot_ctx.region_labels(:).', ', ')), ...
    'Interpreter', 'none');
end


function target_region_names = local_target_region_names_from_observable(Obs)
target_region_names = {};
if isfield(Obs, 'cfg') && isfield(Obs.cfg, 'bold') && ...
        isfield(Obs.cfg.bold, 'selected_region_names')
    target_region_names = cellstr(string(Obs.cfg.bold.selected_region_names(:)).');
end
if isempty(target_region_names) && isfield(Obs, 'params') && ...
        isfield(Obs.params, 'observable_branch')
    branch = lower(char(string(Obs.params.observable_branch)));
    switch branch
        case {'hp', 'hp_svd100'}
            target_region_names = {'HP'};
        case {'elehp'}
            target_region_names = {'eleHP'};
        case {'global_svd100', 'svd', 'roi_mean', 'slow_band_power'}
            target_region_names = {};
    end
end
end


function [coords, region_labels, n_var_region, ana] = local_collect_roi_coords_and_anatomy(roiTs, target_region_names)
coords_cells = {};
region_labels = {};
n_var_region = [];
ana = [];
target_region_names = cellstr(string(target_region_names(:)).');
for i_region = 1:size(roiTs, 2)
    R = roiTs{1, i_region};
    if ~isstruct(R) || ~isfield(R, 'coords') || isempty(R.coords)
        continue;
    end
    if ~isfield(R, 'dat') || size(R.dat, 2) ~= size(R.coords, 1)
        continue;
    end
    label = local_roi_label(R, i_region);
    if ~isempty(target_region_names) && ...
            ~any(strcmpi(label, target_region_names))
        continue;
    end
    coords_cells{end+1, 1} = double(R.coords); %#ok<AGROW>
    n_var_region(end+1, 1) = size(R.coords, 1); %#ok<AGROW>
    region_labels{end+1, 1} = label; %#ok<AGROW>
    if isempty(ana) && isfield(R, 'ana') && ~isempty(R.ana)
        ana = R.ana;
    end
end
if isempty(coords_cells)
    error('No ROI coordinates were found in roiTs for target regions: %s', ...
        strjoin(target_region_names, ', '));
end
if isempty(ana)
    error('No anatomy field .ana was found in roiTs.');
end
coords = cat(1, coords_cells{:});
end


function label = local_roi_label(R, i_region)
if isfield(R, 'name') && ~isempty(R.name)
    label = char(string(R.name));
else
    label = sprintf('roi%02d', i_region);
end
end


function vals = local_prepare_values(raw_vals, value_mode)
switch lower(value_mode)
    case 'abs'
        vals = local_minmax(abs(double(raw_vals)));
    case 'real'
        vals = real(double(raw_vals));
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


function rgb = local_gray_rgb(img, lim)
x = (double(img) - lim(1)) ./ max(eps, lim(2) - lim(1));
x = min(max(x, 0), 1);
rgb = repmat(x, 1, 1, 3);
end


function lim = local_robust_limits(x)
x = x(isfinite(x));
if isempty(x)
    lim = [0 1];
    return;
end
lim = [prctile(x, 1), prctile(x, 99.5)];
if lim(2) <= lim(1)
    lim = [min(x), max(x)];
end
if lim(2) <= lim(1)
    lim = [0 1];
end
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
cmap = min(max(0.15 + 0.85 * cmap, 0), 1);
end


function mode = local_parse_observable_mode(run_name)
patterns = {'slow_band_power', 'global_svd100', 'HP_svd100', ...
    'roi_mean', 'eleHP', 'HP', 'svd', 'abs', 'complex_split', ...
    'complex', 'identity'};
mode = '';
for i = 1:numel(patterns)
    pat = patterns{i};
    if endsWith(run_name, pat, 'IgnoreCase', true)
        mode = pat;
        return;
    end
end
end


function form = local_parse_residual_form(run_name)
forms = {'projected_kv', 'projected_vlambda'};
form = '';
for i = 1:numel(forms)
    if contains(run_name, forms{i}, 'IgnoreCase', true)
        form = forms{i};
        return;
    end
end
end


function dataset_dir = local_dataset_dir_name(dataset_stem)
switch lower(dataset_stem)
    case 'e10fv1'
        dataset_dir = 'E10.fV1';
    case 'e10gb1'
        dataset_dir = 'E10.gb1';
    case 'e10gh1'
        dataset_dir = 'E10.gH1';
    case 'f12m01'
        dataset_dir = 'F12.m01';
    otherwise
        dataset_dir = dataset_stem;
end
end


function cand = local_empty_candidate()
cand = struct('dataset_stem', '', 'run_name', '', ...
    'observable_mode', '', 'residual_form', '', ...
    'bold_post_file', '', 'xcorr_file', '');
end


function row = local_empty_row()
row = struct('dataset_stem', '', 'run_name', '', ...
    'observable_mode', '', 'residual_form', '', ...
    'status', '', 'message', '', 'bold_post_file', '', ...
    'xcorr_file', '', 'activation_dir', '', ...
    'n_requested_maps', 0, 'n_saved_maps', 0, 'runtime_sec', 0);
end


function row = local_make_row(cand, status, message, activation_dir, ...
        n_requested_maps, n_saved_maps, runtime_sec)
row = local_empty_row();
row.dataset_stem = cand.dataset_stem;
row.run_name = cand.run_name;
row.observable_mode = cand.observable_mode;
row.residual_form = cand.residual_form;
row.status = status;
row.message = message;
row.bold_post_file = cand.bold_post_file;
row.xcorr_file = cand.xcorr_file;
row.activation_dir = activation_dir;
row.n_requested_maps = n_requested_maps;
row.n_saved_maps = n_saved_maps;
row.runtime_sec = runtime_sec;
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end
