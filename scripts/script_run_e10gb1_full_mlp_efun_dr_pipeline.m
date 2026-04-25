% Run the complete e10gb1 MLP ResKoopNet output through the efun DR pipeline.
%
% Output layout:
%   E:\DataPons_processed\e10gb1\efun\mlp260415_kvabs\<method>\
%       mat     compact reduction result
%       fig     whole-record figures
%       top30   top consensus-state-diversity window figures
%       peaks   consensus-state peak-distribution analysis
%
% Methods enabled by default: SVD, logSVD, NMF, MDS, UMAP.

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

base_cfg = cfg_eigenfunction_reduction_minimal();

run_opts = struct();
run_opts.dataset_name = 'e10gb1';
run_opts.source_run_name = 'mlp_obs_e10gb1_260415_shuffle_plateau_projected_kv_abs';
run_opts.output_run_name = 'mlp260415_kvabs';
run_opts.source_data_dir = fullfile('E:\autodl_results', ...
    run_opts.dataset_name, 'mlp', 'outputs', run_opts.source_run_name);
run_opts.output_root = fullfile(get_project_processed_root(), ...
    run_opts.dataset_name, 'efun', run_opts.output_run_name);

run_opts.sweep_name = 'dr';
run_opts.continue_on_error = true;
run_opts.skip_unavailable_methods = true;
run_opts.close_figures = true;

run_opts.make_overview = true;
run_opts.make_state_space = true;              % One whole-record ss plot per method.
run_opts.make_consensus_state_space = true;    % Single-patch renderer; safe for long runs.
run_opts.make_spectrum_diagnostics = true;
run_opts.make_top30 = true;
run_opts.make_peak_analysis = true;

run_opts.top30_n_windows = 30;
run_opts.top30_force_recompute_norm_cache = false;

run_opts.peak_component_source = 'smooth_if_available';
run_opts.peak_mode = 'max';                    % 'max' | 'max_abs'
run_opts.peak_baseline_mode = 'matched';
run_opts.peak_random_seed = 1;
run_opts.peak_alpha = 0.05;

method_specs = local_default_method_specs();

base_cfg = local_configure_full_mlp_cfg(base_cfg, run_opts);
base_cfg = local_attach_umap_if_available(base_cfg, repo_root);

local_validate_source_chunks(base_cfg);

if exist(run_opts.output_root, 'dir') ~= 7
    mkdir(run_opts.output_root);
end

consensus_loader_cfg = struct();
consensus_loader_cfg.file_stem = base_cfg.dataset.name;
[C_consensus, source_consensus_file] = load_consensus_state_results( ...
    consensus_loader_cfg, base_cfg.dataset.processed_root, []);

manifest_file = fullfile(run_opts.output_root, ...
    sprintf('%s_%s.csv', base_cfg.dataset.name, run_opts.sweep_name));
manifest_rows = repmat(local_empty_manifest_row(), numel(method_specs), 1);
manifest_count = 0;
source_cache = local_empty_source_cache();

fprintf('Running complete MLP ResKoopNet efun DR pipeline.\n');
fprintf('Source chunks:\n  %s\n', base_cfg.source.data_dir);
fprintf('Output root:\n  %s\n', run_opts.output_root);
if ~isempty(base_cfg.path.spectrum.umap_dir)
    fprintf('UMAP directory:\n  %s\n', base_cfg.path.spectrum.umap_dir);
end
fprintf('Consensus states:\n  %s\n', source_consensus_file);
fprintf('Manifest:\n  %s\n\n', manifest_file);

for i = 1:numel(method_specs)
    spec = method_specs(i);
    method_tag = local_method_tag(spec);
    if ~spec.enabled
        continue;
    end

    [skip_method, skip_reason] = local_should_skip_method( ...
        spec, base_cfg, run_opts.skip_unavailable_methods);
    if skip_method
        fprintf('[SKIP] %s: %s\n', method_tag, skip_reason);
        manifest_count = manifest_count + 1;
        manifest_rows(manifest_count) = local_build_manifest_row( ...
            spec, method_tag, 'skipped', skip_reason, '', '', '', '', ...
            '', '', '', '', '', NaN, NaN);
        local_write_manifest(manifest_rows(1:manifest_count), manifest_file);
        continue;
    end

    method_cfg = local_configure_method_cfg(base_cfg, spec, run_opts);
    method_cfg = local_attach_source_cache(method_cfg, source_cache);
    fprintf('\n[%d/%d] %s\n', i, numel(method_specs), method_tag);
    fprintf('Method output:\n  %s\n', method_cfg.output.root);

    t_start = tic;
    try
        [result, EDMD_outputs_cache, concat_info_cache, source_info_cache] = ...
            run_eigenfunction_reduction_pipeline(method_cfg);
        if ~source_cache.is_ready
            source_cache = local_store_source_cache( ...
                EDMD_outputs_cache, concat_info_cache, source_info_cache);
        end

        plot_paths = local_plot_method_outputs( ...
            result, method_cfg, run_opts, C_consensus);

        top30_manifest_file = '';
        if run_opts.make_top30
            top30_manifest_file = local_run_top30_window_plots( ...
                method_cfg, result.artifacts.result_mat_file, run_opts);
        end

        peak_paths = local_empty_peak_paths();
        if run_opts.make_peak_analysis
            peak_paths = local_run_peak_analysis( ...
                result, method_cfg, run_opts, C_consensus, source_consensus_file);
        end

        runtime_sec = toc(t_start);
        recon_error = local_get_reconstruction_error(result);

        manifest_count = manifest_count + 1;
        manifest_rows(manifest_count) = local_build_manifest_row( ...
            spec, method_tag, 'ok', '', result.artifacts.result_mat_file, ...
            plot_paths.overview_png, plot_paths.state_space_png, ...
            plot_paths.consensus_state_space_png, ...
            plot_paths.spectrum_diagnostics_png, top30_manifest_file, ...
            peak_paths.main_mat, peak_paths.stats_csv, ...
            peak_paths.peak_distribution_png, runtime_sec, recon_error);

        fprintf('[OK] %s finished in %.1f sec.\n', method_tag, runtime_sec);
        fprintf('Result MAT:\n  %s\n', result.artifacts.result_mat_file);

        clear result
    catch ME
        runtime_sec = toc(t_start);
        fprintf(2, '[ERROR] %s failed after %.1f sec:\n  %s\n', ...
            method_tag, runtime_sec, ME.message);

        manifest_count = manifest_count + 1;
        manifest_rows(manifest_count) = local_build_manifest_row( ...
            spec, method_tag, 'error', ME.message, '', '', '', '', ...
            '', '', '', '', '', runtime_sec, NaN);

        if ~run_opts.continue_on_error
            rethrow(ME);
        end
    end

    local_write_manifest(manifest_rows(1:manifest_count), manifest_file);
end

local_write_manifest(manifest_rows(1:manifest_count), manifest_file);

fprintf('\nComplete MLP efun DR pipeline finished.\n');
fprintf('Manifest:\n  %s\n', manifest_file);


function cfg = local_configure_full_mlp_cfg(cfg, run_opts)
cfg.dataset.name = run_opts.dataset_name;
cfg.source.mode = 'chunk_dir';
cfg.source.data_dir = run_opts.source_data_dir;
cfg.source.concat.scan_mode = 'uniform_except_last';
cfg.output.root = run_opts.output_root;
cfg.output.figure_dir = fullfile(run_opts.output_root, 'fig');
cfg.output.result_dir = fullfile(run_opts.output_root, 'mat');
cfg.output.source_run_name = run_opts.source_run_name;
cfg.output.output_run_name = run_opts.output_run_name;
cfg.save.dir = cfg.output.result_dir;
cfg.save.file_stem = sprintf('%s_efun', run_opts.dataset_name);
cfg.save.tag = 'min';
cfg.save.payload = 'compact';
end


function cfg = local_attach_umap_if_available(cfg, repo_root)
umap_dir = local_find_default_umap_dir(repo_root);
if ~isempty(umap_dir) && isempty(cfg.path.spectrum.umap_dir)
    cfg.path.spectrum.umap_dir = umap_dir;
end
end


function local_validate_source_chunks(cfg)
if exist(cfg.source.data_dir, 'dir') ~= 7
    error('Source chunk directory does not exist:\n  %s', cfg.source.data_dir);
end

chunk_files = dir(fullfile(cfg.source.data_dir, cfg.source.concat.filename_pattern));
if isempty(chunk_files)
    error('No EDMD chunk files matched %s in:\n  %s', ...
        cfg.source.concat.filename_pattern, cfg.source.data_dir);
end

fprintf('Found %d EDMD chunk files.\n', numel(chunk_files));
end


function cache = local_empty_source_cache()
cache = struct();
cache.is_ready = false;
cache.EDMD_outputs = [];
cache.concat_info = [];
cache.source_info = [];
end


function cfg = local_attach_source_cache(cfg, cache)
if ~cache.is_ready
    return;
end

cfg.source.mode = 'preloaded';
cfg.source.preloaded_EDMD_outputs = cache.EDMD_outputs;
cfg.source.preloaded_concat_info = cache.concat_info;
cfg.source.preloaded_source_info = cache.source_info;
end


function cache = local_store_source_cache(EDMD_outputs, concat_info, source_info)
cache = local_empty_source_cache();
cache.is_ready = true;
cache.EDMD_outputs = EDMD_outputs;
cache.concat_info = concat_info;
cache.source_info = source_info;
end


function specs = local_default_method_specs()
specs = repmat(local_empty_method_spec(), 5, 1);
specs(1) = local_make_method_spec(true, 'time', 'SVD', struct());
specs(2) = local_make_method_spec(true, 'time', 'logSVD', struct());
specs(3) = local_make_method_spec(true, 'time', 'NMF', struct());
specs(4) = local_make_method_spec(true, 'spectrum', 'MDS', struct());
specs(5) = local_make_method_spec(true, 'spectrum', 'UMAP', struct());
end


function spec = local_empty_method_spec()
spec = struct();
spec.enabled = false;
spec.kind = '';
spec.method = '';
spec.options = struct();
end


function spec = local_make_method_spec(enabled, kind, method, options)
spec = local_empty_method_spec();
spec.enabled = enabled;
spec.kind = char(kind);
spec.method = char(method);
spec.options = options;
end


function umap_dir = local_find_default_umap_dir(repo_root)
umap_dir = '';
icpbr_root = fileparts(fileparts(fileparts(repo_root)));
candidate_dirs = { ...
    fullfile(icpbr_root, 'forNikos', 'umapFileExchange (4.4)', 'umap'), ...
    fullfile(icpbr_root, 'forNikos', 'net_fmri_tutorials', ...
    'umapFileExchange (4.4)', 'umap')};

for i = 1:numel(candidate_dirs)
    candidate = candidate_dirs{i};
    if exist(fullfile(candidate, 'run_umap.m'), 'file') == 2
        umap_dir = candidate;
        return;
    end
end
end


function [skip_method, reason] = local_should_skip_method(spec, cfg, skip_unavailable)
skip_method = false;
reason = '';
if ~skip_unavailable
    return;
end

switch lower(spec.kind)
    case 'time'
        if strcmpi(spec.method, 'NMF') && exist('nnmf', 'file') ~= 2
            skip_method = true;
            reason = 'nnmf is not available on the MATLAB path.';
        end

    case 'spectrum'
        if strcmpi(cfg.path.spectrum.cluster_method, 'kmeans') && ...
                exist('kmeans', 'file') ~= 2
            skip_method = true;
            reason = 'kmeans is not available on the MATLAB path.';
            return;
        end

        switch lower(spec.method)
            case 'mds'
                if exist('mdscale', 'file') ~= 2
                    skip_method = true;
                    reason = 'mdscale is not available on the MATLAB path.';
                end
            case 'umap'
                if ~isempty(cfg.path.spectrum.umap_dir) && ...
                        exist(cfg.path.spectrum.umap_dir, 'dir') == 7
                    addpath(cfg.path.spectrum.umap_dir);
                end
                if exist('run_umap', 'file') ~= 2
                    skip_method = true;
                    reason = 'run_umap is not available; set cfg.path.spectrum.umap_dir first.';
                end
        end

    otherwise
        skip_method = true;
        reason = sprintf('Unknown method kind "%s".', spec.kind);
end
end


function method_cfg = local_configure_method_cfg(base_cfg, spec, run_opts)
method_cfg = base_cfg;
method_tag = local_method_tag(spec);
method_root = fullfile(run_opts.output_root, method_tag);
figure_dir = fullfile(method_root, 'fig');
result_dir = fullfile(method_root, 'mat');

method_cfg.path.kind = spec.kind;
switch lower(spec.kind)
    case 'time'
        method_cfg.path.time.method = spec.method;
        method_cfg = local_apply_options(method_cfg, 'path.time', spec.options);
    case 'spectrum'
        method_cfg.path.spectrum.method = spec.method;
        method_cfg = local_apply_options(method_cfg, 'path.spectrum', spec.options);
    otherwise
        error('Unknown method kind "%s".', spec.kind);
end

method_cfg.output.root = method_root;
method_cfg.output.figure_dir = figure_dir;
method_cfg.output.result_dir = result_dir;
method_cfg.save.dir = result_dir;
method_cfg.save.tag = 'sw';
method_cfg.save.payload = 'compact';

method_cfg.viz.overview.enable = run_opts.make_overview;
method_cfg.viz.overview.save_dir = figure_dir;
method_cfg.viz.overview.save_tag = [method_tag, '_ov'];

method_cfg.viz.state_space.enable = run_opts.make_state_space;
method_cfg.viz.state_space.save_dir = figure_dir;
method_cfg.viz.state_space.save_tag = [method_tag, '_ss'];

method_cfg.viz.state_space_consensus.enable = run_opts.make_consensus_state_space;
method_cfg.viz.state_space_consensus.save_dir = figure_dir;
method_cfg.viz.state_space_consensus.save_tag = [method_tag, '_ssc'];

method_cfg.viz.spectrum.enable = run_opts.make_spectrum_diagnostics;
method_cfg.viz.spectrum.save_dir = figure_dir;
method_cfg.viz.spectrum.save_tag = [method_tag, '_spec'];
end


function cfg = local_apply_options(cfg, dotted_path, options)
if isempty(fieldnames(options))
    return;
end

parts = strsplit(dotted_path, '.');
target = cfg.(parts{1}).(parts{2});
names = fieldnames(options);
for i = 1:numel(names)
    target.(names{i}) = options.(names{i});
end
cfg.(parts{1}).(parts{2}) = target;
end


function plot_paths = local_plot_method_outputs(result, cfg, run_opts, C_consensus)
plot_paths = struct();
plot_paths.overview_png = '';
plot_paths.state_space_png = '';
plot_paths.consensus_state_space_png = '';
plot_paths.spectrum_diagnostics_png = '';

if run_opts.make_overview
    stage_tic = local_stage_start('Plotting overview figure');
    [fig, info] = plot_eigenfunction_component_overview(result, cfg.viz.overview);
    if isfield(info, 'save_path')
        plot_paths.overview_png = info.save_path;
    end
    local_close_if_requested(fig, run_opts);
    local_stage_done('Plotted overview figure', stage_tic, plot_paths.overview_png);
end

if run_opts.make_state_space
    stage_tic = local_stage_start('Plotting state-space trajectory');
    [fig, info] = plot_eigenfunction_state_space_trajectory(result, cfg.viz.state_space);
    if isfield(info, 'save_path')
        plot_paths.state_space_png = info.save_path;
    end
    local_close_if_requested(fig, run_opts);
    local_stage_done('Plotted state-space trajectory', stage_tic, plot_paths.state_space_png);
end

if run_opts.make_consensus_state_space
    stage_tic = local_stage_start('Plotting consensus state-space trajectory');
    [fig, info] = plot_eigenfunction_state_space_consensus_trajectory( ...
        result, C_consensus, cfg.viz.state_space_consensus);
    if isfield(info, 'save_path')
        plot_paths.consensus_state_space_png = info.save_path;
    end
    local_close_if_requested(fig, run_opts);
    local_stage_done('Plotted consensus state-space trajectory', ...
        stage_tic, plot_paths.consensus_state_space_png);
end

if strcmpi(result.meta.path_kind, 'spectrum') && run_opts.make_spectrum_diagnostics
    stage_tic = local_stage_start('Plotting spectrum diagnostics');
    [fig, info] = plot_eigenfunction_spectrum_diagnostics(result, cfg.viz.spectrum);
    if isfield(info, 'save_path')
        plot_paths.spectrum_diagnostics_png = info.save_path;
    end
    local_close_if_requested(fig, run_opts);
    local_stage_done('Plotted spectrum diagnostics', ...
        stage_tic, plot_paths.spectrum_diagnostics_png);
end
end


function top30_manifest_file = local_run_top30_window_plots(method_cfg, result_file, run_opts) %#ok<INUSD>
top30_manifest_file = '';
top30_script = fullfile(method_cfg.repo_root, 'scripts', ...
    'script_plot_e10gb1_efun_svd_top30_state_diversity_windows.m');
if exist(top30_script, 'file') ~= 2
    error('Top30 plotting script was not found:\n  %s', top30_script);
end

cfg = method_cfg; %#ok<NASGU>
top_window_params = struct();
top_window_params.n_top_windows = run_opts.top30_n_windows;
top_window_params.make_overview = true;
top_window_params.make_state_space = false;
top_window_params.make_consensus_state_space = true;
top_window_params.close_figures = run_opts.close_figures;
top_window_params.force_recompute_norm_cache = ...
    run_opts.top30_force_recompute_norm_cache;
top_window_params.save_dir = fullfile(method_cfg.output.root, 'top30');
top_window_params.norm_cache_file = fullfile(top_window_params.save_dir, 'norm.mat'); %#ok<STRNU>

stage_tic = local_stage_start('Running top30 window plots');
run(top30_script);

if exist('manifest_file', 'var') == 1 && ~isempty(manifest_file)
    top30_manifest_file = manifest_file;
end
local_stage_done('Finished top30 window plots', stage_tic, top30_manifest_file);
end


function peak_paths = local_run_peak_analysis( ...
    result, method_cfg, run_opts, C_consensus, source_consensus_file)
stage_tic = local_stage_start('Running consensus-state peak analysis');
peak_params = struct();
peak_params.component_source = run_opts.peak_component_source;
peak_params.state_start_idx = method_cfg.viz.state_space_consensus.state_start_idx;
peak_params.peak_mode = run_opts.peak_mode;
peak_params.baseline_mode = run_opts.peak_baseline_mode;
peak_params.baseline_code = method_cfg.viz.state_space_consensus.baseline_code;
peak_params.baseline_label = method_cfg.viz.state_space_consensus.baseline_label;
peak_params.min_window_samples = 1;
peak_params.baseline_random_seed = run_opts.peak_random_seed;
peak_params.alpha = run_opts.peak_alpha;
peak_params.save_dir = fullfile(method_cfg.output.root, 'peaks');
peak_params.save_tag = sprintf('%s_peaks', method_cfg.dataset.name);
peak_params.save_results = true;
peak_params.write_csv = true;
peak_params.save_figures = true;
peak_params.close_figures = run_opts.close_figures;

A = analyze_eigenfunction_component_peaks_by_consensus_state( ...
    result, C_consensus, peak_params, ...
    result.artifacts.result_mat_file, source_consensus_file);

peak_paths = local_empty_peak_paths();
if isfield(A, 'save_paths')
    fields = fieldnames(peak_paths);
    for i = 1:numel(fields)
        if isfield(A.save_paths, fields{i})
            peak_paths.(fields{i}) = A.save_paths.(fields{i});
        end
    end
end
local_stage_done('Finished consensus-state peak analysis', ...
    stage_tic, peak_paths.main_mat);
end


function peak_paths = local_empty_peak_paths()
peak_paths = struct();
peak_paths.main_mat = '';
peak_paths.stats_csv = '';
peak_paths.peak_distribution_png = '';
end


function stage_tic = local_stage_start(label)
stage_tic = tic;
fprintf('[stage] %s...\n', label);
end


function local_stage_done(label, stage_tic, save_path)
if nargin < 3
    save_path = '';
end

fprintf('[stage] %s in %s.\n', label, local_format_seconds(toc(stage_tic)));
if ~isempty(save_path)
    fprintf('        %s\n', save_path);
end
end


function txt = local_format_seconds(seconds_in)
if ~isfinite(seconds_in)
    txt = '--';
    return;
end

seconds_in = max(0, seconds_in);
if seconds_in < 1
    txt = sprintf('%.3fs', seconds_in);
elseif seconds_in < 60
    txt = sprintf('%.1fs', seconds_in);
else
    hours = floor(seconds_in / 3600);
    minutes = floor(mod(seconds_in, 3600) / 60);
    seconds = floor(mod(seconds_in, 60));
    if hours > 0
        txt = sprintf('%02d:%02d:%02d', hours, minutes, seconds);
    else
        txt = sprintf('%02d:%02d', minutes, seconds);
    end
end
end


function local_close_if_requested(fig, run_opts)
if run_opts.close_figures && ~isempty(fig) && isvalid(fig)
    close(fig);
end
end


function method_tag = local_method_tag(spec)
switch lower(spec.method)
    case {'svd', 'logsvd', 'nmf', 'mds', 'umap'}
        method_tag = lower(spec.method);
    otherwise
        method_tag = sprintf('%s_%s', lower(spec.kind), lower(spec.method));
end
method_tag = regexprep(method_tag, '[^\w\-]+', '_');
end


function recon_error = local_get_reconstruction_error(result)
recon_error = NaN;
if isfield(result, 'quality') && isfield(result.quality, 'reconstruction_error_fro')
    recon_error = double(result.quality.reconstruction_error_fro);
end
end


function row = local_build_manifest_row(spec, method_tag, status, message, ...
    result_mat_file, overview_png, state_space_png, consensus_state_space_png, ...
    spectrum_diagnostics_png, top30_manifest_file, peak_mat_file, ...
    peak_stats_csv, peak_distribution_png, runtime_sec, reconstruction_error_fro)
row = struct();
row.method_tag = string(method_tag);
row.kind = string(spec.kind);
row.method = string(spec.method);
row.status = string(status);
row.message = string(message);
row.result_mat_file = string(result_mat_file);
row.overview_png = string(overview_png);
row.state_space_png = string(state_space_png);
row.consensus_state_space_png = string(consensus_state_space_png);
row.spectrum_diagnostics_png = string(spectrum_diagnostics_png);
row.top30_manifest_file = string(top30_manifest_file);
row.peak_mat_file = string(peak_mat_file);
row.peak_stats_csv = string(peak_stats_csv);
row.peak_distribution_png = string(peak_distribution_png);
row.runtime_sec = runtime_sec;
row.reconstruction_error_fro = reconstruction_error_fro;
end


function row = local_empty_manifest_row()
row = local_build_manifest_row(local_empty_method_spec(), '', '', '', ...
    '', '', '', '', '', '', '', '', '', NaN, NaN);
end


function local_write_manifest(rows, manifest_file)
if isempty(rows)
    return;
end

manifest_dir = fileparts(manifest_file);
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end

T = struct2table(rows);
writetable(T, manifest_file);
end
