% Run one BOLD ResKoopNet result through:
%   1) EDMD/BOLD postprocessing with deconvolved eigenfunctions and timescale plots
%   2) density cross-correlation with explicit session-border masking
%   3) reference-style activation maps for the top correlated BOLD eigenbasis
%
% This script is intentionally single-run and editable. It does not modify cfg.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

%% -------------------- user settings --------------------
result_dir = ['E:\autodl_results_local\bold_wsl\e10gb1\mlp\outputs\', ...
    'mlp_obs_bold_wsl_20260423_e10gb1_projected_kv_HP_svd100'];

dataset_stem = 'e10gb1';
dataset_id = 'E10.gb1';
observable_mode = 'HP_svd100';

observable_mat = fullfile('E:\DataPons_processed\bold_observables', ...
    dataset_id, sprintf('%s_bold_observables_%s.mat', dataset_id, observable_mode));
roi_ts_mat = fullfile('E:\DataPons', dataset_id, 'roits', ...
    sprintf('%s_0001_roits.mat', dataset_stem));

density_sources = struct([]);
density_sources(1).name = 'blp_event_density';
density_sources(1).type = 'event_density';
density_sources(1).file = fullfile(io_project.get_project_processed_root(), ...
    dataset_stem, 'event_density', sprintf('%s_event_density_2s.mat', dataset_stem));

% Optional: add threshold/dimred density files or manifests here.
% density_sources(2).name = 'efun_threshold_density';
% density_sources(2).type = 'threshold_density_scan';
% density_sources(2).manifest_file = 'E:\DataPons_processed\postprocessing_manifests\efun_den\...mat';

max_lag_sec = 120;
border_mask_sec = 120;
top_n = 5;

activation_value_mode = 'abs';    % 'abs' or 'real'
activation_slice_list = 1:20;
activation_tiles_per_row = 10;

%% -------------------- derived paths --------------------
if exist(result_dir, 'dir') ~= 7
    error('result_dir does not exist: %s', result_dir);
end

[~, run_name] = fileparts(result_dir);
out_root = fullfile(io_project.get_project_processed_root(), dataset_stem, ...
    'bold_reskoopnet_single_run_analysis', run_name);
mat_dir = fullfile(out_root, 'mat');
fig_dir = fullfile(out_root, 'fig');
act_dir = fullfile(fig_dir, 'activation_maps_top5');
if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end
if exist(act_dir, 'dir') ~= 7, mkdir(act_dir); end

%% -------------------- load and postprocess --------------------
source_cfg = struct();
source_cfg.mode = 'chunk_dir';
source_cfg.data_dir = result_dir;
source_cfg.concat = struct();
source_cfg.concat.filename_pattern = '*_outputs_*.mat';
source_cfg.concat.variable_name = 'EDMD_outputs';
source_cfg.concat.concat_fields = {'efuns'};
source_cfg.concat.scan_mode = 'uniform_except_last';
source_cfg.concat.verbose = true;

[EDMD_outputs, concat_info, source_info] = io_edmd.load_edmd_source(source_cfg);
Obs = load(observable_mat);

session = local_session_from_observable(Obs, size(EDMD_outputs.efuns, 1));
dt = local_resolve_dt_from_observable(EDMD_outputs, Obs, 2);
t = (0:size(EDMD_outputs.efuns, 1)-1)' * dt;
session_border_t = t(session.border_idx(session.border_idx >= 1 & session.border_idx <= numel(t)));

post_opts = struct();
post_opts.abs_thresh = 0.01;
post_opts.sort_by = 'modulus';
post_opts.sort_dir = 'descend';
post_opts.max_basis = 80;
post_opts.dt = dt;
post_opts.do_plot = true;
post_opts.time_vec = t;
post_opts.session_border = session_border_t;
post_opts.draw_border = true;
[EDMD_outputs, fig_main] = postprocess_EDMD_outputs(EDMD_outputs, post_opts);
exportgraphics(fig_main, fullfile(fig_dir, 'bold_post_eigenfunctions.png'), 'Resolution', 200);

deconv_cfg = struct();
deconv_cfg.method = 'koopman_residual';
deconv_cfg.lambda_source = 'edmd';
deconv_cfg.lambdaType = 'discrete';
deconv_cfg.dt = dt;
deconv_cfg.max_modes_all = 80;
deconv_cfg.max_modes_sel = 40;
deconv_cfg.do_plot = true;
deconv_cfg.time_vec = t;
deconv_cfg.session_border = session_border_t;
deconv_cfg.draw_border = true;
[EDMD_outputs, fig_deconv, deconv] = postprocess_EDMD_outputs_deconv_efuns(EDMD_outputs, deconv_cfg);
exportgraphics(fig_deconv, fullfile(fig_dir, 'bold_post_deconvolved_eigenfunctions.png'), 'Resolution', 200);

timescale_cfg = struct();
timescale_cfg.dt = dt;
timescale_cfg.max_modes_all = 80;
timescale_cfg.max_modes_sel = 40;
timescale_cfg.maxLag = 200;
timescale_cfg.xlim_time = 240;
timescale_cfg.title_prefix = run_name;
[fig_timescale, timescale_info] = postprocess_EDMD_outputs_timescale(EDMD_outputs, timescale_cfg);
exportgraphics(fig_timescale, fullfile(fig_dir, 'bold_post_timescale.png'), 'Resolution', 200);

BOLD_POST = struct();
BOLD_POST.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
BOLD_POST.run_info = struct('dataset_stem', dataset_stem, ...
    'dataset_id', dataset_id, 'run_name', run_name, ...
    'output_dir', result_dir, 'observable_mode', observable_mode);
BOLD_POST.source_info = source_info;
BOLD_POST.concat_info = concat_info;
BOLD_POST.observable_file = observable_mat;
BOLD_POST.observable = Obs;
BOLD_POST.session = session;
BOLD_POST.dt = dt;
BOLD_POST.time_vec = t;
BOLD_POST.session_border_t = session_border_t;
BOLD_POST.EDMD_outputs = EDMD_outputs;
BOLD_POST.deconv = deconv;
BOLD_POST.timescale_info = timescale_info;
BOLD_POST.artifacts = struct();
BOLD_POST.artifacts.post_file = fullfile(mat_dir, [run_name, '_bold_post.mat']);
save(BOLD_POST.artifacts.post_file, 'BOLD_POST', '-v7.3');

%% -------------------- density cross-correlation --------------------
xcorr_params = struct();
xcorr_params.max_lag_sec = max_lag_sec;
xcorr_params.border_mask_sec = border_mask_sec;
xcorr_params.min_valid_samples = 20;
xcorr_params.top_n = top_n;
xcorr_params.feature_names = {'efun_abs', 'efun_real', 'deconv_abs', 'deconv_real'};
xcorr_params.save_dir = fullfile(out_root, 'density_cross_correlation');
xcorr_params.save_tag = 'bold_efun_density_xcorr';
xcorr_params.save_results = true;
xcorr_params.make_figures = true;

xcorr_out = compute_bold_efun_density_cross_correlation( ...
    BOLD_POST, density_sources, xcorr_params);

%% -------------------- activation maps for top correlated eigenbasis --------------------
top_table = xcorr_out.top_table;
if isempty(top_table)
    warning('No top correlated eigenbasis found; skipping activation maps.');
else
    raw_outputs = local_load_first_edmd_output(result_dir);
    raw_indices = local_top_raw_indices(top_table, EDMD_outputs, top_n);
    raw_indices = unique(raw_indices(:).', 'stable');

    fprintf('\nPlotting activation maps for raw eigenbasis indices: %s\n', mat2str(raw_indices));
    for i = 1:numel(raw_indices)
        raw_idx = raw_indices(i);
        meta = local_top_meta_for_raw_index(top_table, EDMD_outputs, raw_idx);
        fig = local_plot_reference_activation_map( ...
            raw_outputs, Obs, roi_ts_mat, raw_idx, activation_value_mode, ...
            activation_slice_list, activation_tiles_per_row, meta);

        png_file = fullfile(act_dir, sprintf('top%02d_raw%04d_%s_activation_reference.png', ...
            i, raw_idx, activation_value_mode));
        fig_file = fullfile(act_dir, sprintf('top%02d_raw%04d_%s_activation_reference.fig', ...
            i, raw_idx, activation_value_mode));
        exportgraphics(fig, png_file, 'Resolution', 220);
        savefig(fig, fig_file);
        close(fig);
    end
end

fprintf('\nFinished one-run BOLD ResKoopNet analysis.\n');
fprintf('Post MAT:\n  %s\n', BOLD_POST.artifacts.post_file);
fprintf('Cross-correlation MAT:\n  %s\n', xcorr_out.save_paths.mat_file);
fprintf('Figures:\n  %s\n', fig_dir);

%% -------------------- local helpers --------------------
function session = local_session_from_observable(Obs, T)
    names = {'session_ids', 'session_lengths', 'session_dx', ...
        'session_start_idx', 'session_end_idx', 'border_idx'};
    session = struct();
    for ii = 1:numel(names)
        name = names{ii};
        if isfield(Obs, name)
            session.(name) = Obs.(name);
        elseif isfield(Obs, 'O') && isfield(Obs.O, name)
            session.(name) = Obs.O.(name);
        else
            session.(name) = [];
        end
    end
    if isempty(session.session_start_idx) || isempty(session.session_end_idx)
        session.session_start_idx = 1;
        session.session_end_idx = T;
        session.session_lengths = T;
        session.session_ids = 1;
    end
    session.session_start_idx = double(session.session_start_idx(:));
    session.session_end_idx = double(session.session_end_idx(:));
    if isempty(session.session_ids)
        session.session_ids = (1:numel(session.session_start_idx)).';
    end
    if isempty(session.border_idx)
        session.border_idx = session.session_end_idx(1:end-1);
    end
end

function dt = local_resolve_dt_from_observable(EDMD_outputs, Obs, fallback)
    dt = [];
    if isfield(EDMD_outputs, 'dt') && ~isempty(EDMD_outputs.dt)
        dt = double(EDMD_outputs.dt);
    elseif isfield(EDMD_outputs, 'dx') && ~isempty(EDMD_outputs.dx)
        dt = double(EDMD_outputs.dx);
    elseif isfield(Obs, 'dx') && ~isempty(Obs.dx)
        dt = double(Obs.dx);
    elseif isfield(Obs, 'O') && isfield(Obs.O, 'dx') && ~isempty(Obs.O.dx)
        dt = double(Obs.O.dx);
    end
    if isempty(dt) || ~isfinite(dt(1)) || dt(1) <= 0
        dt = fallback;
    else
        dt = dt(1);
    end
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

function fig = local_plot_reference_activation_map(EDMD_outputs, Obs, roi_ts_mat, ...
        raw_index, value_mode, slice_list, tiles_per_row, meta)
    evalues = EDMD_outputs.evalues(:);
    kpm_modes = EDMD_outputs.kpm_modes;
    if raw_index < 1 || raw_index > size(kpm_modes, 1)
        error('raw_index=%d is outside [1, %d].', raw_index, size(kpm_modes, 1));
    end
    lambda = evalues(raw_index);
    mode_obs = double(kpm_modes(raw_index, :));

    O = Obs.O;
    if isfield(O, 'model') && isfield(O.model, 'coeff') && ...
            size(O.model.coeff, 2) == numel(mode_obs)
        voxel_values = mode_obs * double(O.model.coeff).';
        map_space = 'SVD back-projected to voxels';
    else
        voxel_values = mode_obs;
        map_space = 'direct voxel mode';
    end

    R = load(roi_ts_mat, 'roiTs');
    target_region_names = local_target_region_names_from_observable(Obs);
    [coords, region_labels, ~, ana] = local_collect_roi_coords_and_anatomy( ...
        R.roiTs, target_region_names);
    if numel(voxel_values) ~= size(coords, 1)
        error('Voxel-space mode has %d values but roiTs has %d voxel coordinates.', ...
            numel(voxel_values), size(coords, 1));
    end

    vals = local_prepare_values(voxel_values(:), value_mode);
    [cmap, clim_use, value_label] = local_colormap_and_limits(value_mode);
    act_vol = local_values_to_volume(vals, coords, size(ana));
    bg_limits = local_robust_limits(double(ana(:)));

    n_slice = numel(slice_list);
    n_col = tiles_per_row;
    n_row = ceil(n_slice / n_col);
    fig = figure('Color', 'w', 'Name', 'Reference-style BOLD activation map');
    fig.Position = [80 80 1450 170 * n_row + 130];
    tiledlayout(n_row, n_col, 'TileSpacing', 'compact', 'Padding', 'compact');

    for jj = 1:n_slice
        z = slice_list(jj);
        nexttile;
        if z < 1 || z > size(ana, 3)
            axis off;
            title(sprintf('S%d', z));
            continue;
        end
        img = squeeze(ana(:, :, z))';
        overlay = squeeze(act_vol(:, :, z))';
        alpha_data = 0.86 * double(~isnan(overlay));

        image(local_gray_rgb(img, bg_limits));
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
    sgtitle(sprintf('%s | |lambda|=%.4g | %s | %s | ROI(s): %s', ...
        meta.text, abs(lambda), value_mode, map_space, strjoin(region_labels(:).', ', ')), ...
        'Interpreter', 'none');
end

function target_region_names = local_target_region_names_from_observable(Obs)
    target_region_names = {};
    if isfield(Obs, 'cfg') && isfield(Obs.cfg, 'bold') && ...
            isfield(Obs.cfg.bold, 'selected_region_names')
        target_region_names = cellstr(string(Obs.cfg.bold.selected_region_names(:)).');
    end
    if isempty(target_region_names)
        target_region_names = local_region_names_from_observable_metadata(Obs);
    end
end

function target_region_names = local_region_names_from_observable_metadata(Obs)
    target_region_names = {};

    obs_info = [];
    if isfield(Obs, 'observable_info') && istable(Obs.observable_info)
        obs_info = Obs.observable_info;
    elseif isfield(Obs, 'obs_info') && istable(Obs.obs_info)
        obs_info = Obs.obs_info;
    end

    if ~isempty(obs_info) && ismember('region_label', obs_info.Properties.VariableNames)
        labels = cellstr(string(obs_info.region_label(:)));
        labels = labels(~cellfun(@isempty, labels));
        if ~isempty(labels)
            target_region_names = unique(labels, 'stable');
            return;
        end
    end

    if isfield(Obs, 'model') && isstruct(Obs.model) && ...
            isfield(Obs.model, 'input_variable_labels') && ~isempty(Obs.model.input_variable_labels)
        target_region_names = local_infer_region_names_from_variable_labels(Obs.model.input_variable_labels);
    end
end

function target_region_names = local_infer_region_names_from_variable_labels(variable_labels)
    labels = cellstr(string(variable_labels(:)));
    region_names = cell(size(labels));
    keep = false(size(labels));

    for i = 1:numel(labels)
        tok = regexp(labels{i}, '^(.*)_v\d+$', 'tokens', 'once');
        if isempty(tok)
            tok = regexp(labels{i}, '^(.*)_mean$', 'tokens', 'once');
        end
        if isempty(tok) || isempty(tok{1})
            continue;
        end
        region_names{i} = tok{1};
        keep(i) = true;
    end

    if any(keep)
        target_region_names = unique(region_names(keep), 'stable');
    else
        target_region_names = {};
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
    if isempty(ana), error('No anatomy field .ana was found in roiTs.'); end
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
    switch lower(value_mode)
        case 'abs'
            vals = local_minmax(abs(double(raw_vals)));
        case 'real'
            vals = real(double(raw_vals));
            max_abs = max(abs(vals), [], 'omitnan');
            if isfinite(max_abs) && max_abs > 0, vals = vals ./ max_abs; end
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
    if isempty(x), lim = [0 1]; return; end
    lim = [prctile(x, 1), prctile(x, 99.5)];
    if lim(2) <= lim(1), lim = [min(x), max(x)]; end
    if lim(2) <= lim(1), lim = [0 1]; end
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
