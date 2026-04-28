function manifest = postprocess_bold_reskoopnet_results(params)
%POSTPROCESS_BOLD_RESKOOPNET_RESULTS Postprocess completed BOLD ResKoopNet runs.
%
% This wrapper reuses the existing EDMD postprocessing functions for BOLD:
% postprocess_EDMD_outputs, postprocess_EDMD_outputs_deconv_efuns, and
% postprocess_EDMD_outputs_timescale. Spike correlation is intentionally not
% included here.

if nargin < 1 || isempty(params)
    params = struct();
end
params = local_apply_defaults(params);

if params.headless
    set(groot, 'defaultFigureVisible', 'off');
end
try
    set(groot, 'defaultAxesToolbarVisible', 'off');
catch
end

runs = local_discover_runs(params);
if isempty(runs)
    error('No completed BOLD ResKoopNet output folders were found.');
end

if exist(params.manifest_dir, 'dir') ~= 7
    mkdir(params.manifest_dir);
end

rows = repmat(local_empty_row(), numel(runs), 1);
row_count = 0;

fprintf('Discovered %d BOLD ResKoopNet runs.\n', numel(runs));

for i_run = 1:numel(runs)
    run_info = runs(i_run);
    fprintf('\n============================================================\n');
    fprintf('[%d/%d] %s | %s\n', i_run, numel(runs), ...
        run_info.dataset_stem, run_info.run_name);
    fprintf('Source:\n  %s\n', run_info.output_dir);

    t_run = tic;
    try
        out_root = fullfile(params.processed_root, run_info.dataset_stem, ...
            params.output_folder_name, run_info.run_name);
        mat_dir = fullfile(out_root, 'mat');
        fig_dir = fullfile(out_root, 'fig');
        if exist(mat_dir, 'dir') ~= 7, mkdir(mat_dir); end
        if exist(fig_dir, 'dir') ~= 7, mkdir(fig_dir); end

        post_file = fullfile(mat_dir, [run_info.run_name, '_bold_post.mat']);
        if params.skip_existing && exist(post_file, 'file') == 2
            fprintf('Existing postprocessing found; skipping.\n');
            row_count = row_count + 1;
            rows(row_count) = local_make_row(run_info, 'skipped_existing', ...
                '', post_file, '', '', '', toc(t_run));
            continue;
        end

        source_cfg = struct();
        source_cfg.mode = 'chunk_dir';
        source_cfg.data_dir = run_info.output_dir;
        source_cfg.concat = params.concat;
        [EDMD_outputs, concat_info, source_info] = io_edmd.load_edmd_source(source_cfg);

        observable_file = local_find_observable_file(params, run_info, EDMD_outputs);
        obs_meta = local_load_observable_metadata(observable_file);
        session = local_resolve_session_metadata(EDMD_outputs, obs_meta, concat_info);
        dt = local_resolve_dt(EDMD_outputs, obs_meta, params.default_dt);
        [time_vec, session_border_t] = local_time_axis_and_borders( ...
            size(EDMD_outputs.efuns, 1), dt, session);

        post_opts = params.post_opts;
        post_opts.dt = dt;
        post_opts.do_plot = params.make_main_plot;
        post_opts.time_vec = time_vec;
        post_opts.session_border = session_border_t;
        post_opts.draw_border = params.draw_session_borders;
        [EDMD_outputs, fig_main] = postprocess_EDMD_outputs(EDMD_outputs, post_opts);
        main_png = local_save_figure(fig_main, fig_dir, ...
            [run_info.run_name, '_efuns.png'], params);

        deconv_cfg = params.deconv;
        deconv_cfg.dt = dt;
        deconv_cfg.do_plot = params.make_deconv_plot;
        deconv_cfg.time_vec = time_vec;
        deconv_cfg.session_border = session_border_t;
        deconv_cfg.draw_border = params.draw_session_borders;
        deconv_cfg.save_path = '';
        [EDMD_outputs, fig_deconv, deconv] = ...
            postprocess_EDMD_outputs_deconv_efuns(EDMD_outputs, deconv_cfg);
        deconv_png = local_save_figure(fig_deconv, fig_dir, ...
            [run_info.run_name, '_deconv_efuns.png'], params);

        timescale_cfg = params.timescale;
        timescale_cfg.dt = dt;
        timescale_cfg.title_prefix = run_info.run_name;
        [fig_timescale, timescale_info] = ...
            postprocess_EDMD_outputs_timescale(EDMD_outputs, timescale_cfg);
        timescale_png = local_save_figure(fig_timescale, fig_dir, ...
            [run_info.run_name, '_timescale.png'], params);

        BOLD_POST = struct();
        BOLD_POST.created_at = char(datetime('now', 'Format', 'yyyy-MM-dd HH:mm:ss'));
        BOLD_POST.run_info = run_info;
        BOLD_POST.source_info = source_info;
        BOLD_POST.concat_info = concat_info;
        BOLD_POST.observable_file = observable_file;
        BOLD_POST.observable = obs_meta;
        BOLD_POST.session = session;
        BOLD_POST.dt = dt;
        BOLD_POST.time_vec = time_vec;
        BOLD_POST.session_border_t = session_border_t;
        BOLD_POST.EDMD_outputs = EDMD_outputs;
        BOLD_POST.deconv = deconv;
        BOLD_POST.timescale_info = timescale_info;
        BOLD_POST.params = params;
        BOLD_POST.artifacts = struct('post_file', post_file, ...
            'main_png', main_png, 'deconv_png', deconv_png, ...
            'timescale_png', timescale_png);

        save(post_file, 'BOLD_POST', '-v7.3');
        fprintf('Saved BOLD postprocessing:\n  %s\n', post_file);

        row_count = row_count + 1;
        rows(row_count) = local_make_row(run_info, 'ok', '', post_file, ...
            main_png, deconv_png, timescale_png, toc(t_run));
    catch ME
        message = getReport(ME, 'basic', 'hyperlinks', 'off');
        fprintf(2, '[ERROR] %s failed:\n%s\n', run_info.run_name, message);
        row_count = row_count + 1;
        rows(row_count) = local_make_row(run_info, 'error', message, ...
            '', '', '', '', toc(t_run));
        if ~params.continue_on_error
            rethrow(ME);
        end
    end

    if params.close_figures
        close all force;
    end
end

rows = rows(1:row_count);
T = struct2table(rows);
tag = char(datetime('now', 'Format', 'yyyyMMdd_HHmmss'));
manifest = struct();
manifest.params = params;
manifest.runs = runs;
manifest.rows = rows;
manifest.table = T;
manifest.csv_file = fullfile(params.manifest_dir, ...
    ['bold_reskoopnet_postprocessing_manifest_', tag, '.csv']);
manifest.mat_file = fullfile(params.manifest_dir, ...
    ['bold_reskoopnet_postprocessing_manifest_', tag, '.mat']);

writetable(T, manifest.csv_file);
save(manifest.mat_file, 'manifest', '-v7.3');
fprintf('\nBOLD ResKoopNet postprocessing finished.\n');
fprintf('Manifest:\n  %s\n', manifest.csv_file);
end


function params = local_apply_defaults(params)
if ~isfield(params, 'autodl_roots') || isempty(params.autodl_roots)
    params.autodl_roots = {'E:\autodl_results_local\bold_wsl', 'E:\autodl_results\bold'};
end
params.autodl_roots = cellstr(string(params.autodl_roots(:)).');

if ~isfield(params, 'processed_root') || isempty(params.processed_root)
    params.processed_root = io_project.get_project_processed_root();
end
if ~isfield(params, 'dataset_stems') || isempty(params.dataset_stems)
    params.dataset_stems = {};
end
params.dataset_stems = cellstr(string(params.dataset_stems(:)).');
if ~isfield(params, 'observable_modes') || isempty(params.observable_modes)
    params.observable_modes = {};
end
params.observable_modes = cellstr(string(params.observable_modes(:)).');
if ~isfield(params, 'residual_forms') || isempty(params.residual_forms)
    params.residual_forms = {};
end
params.residual_forms = cellstr(string(params.residual_forms(:)).');

if ~isfield(params, 'output_folder_name') || isempty(params.output_folder_name)
    params.output_folder_name = 'bold_reskoopnet_postprocessing';
end
if ~isfield(params, 'manifest_dir') || isempty(params.manifest_dir)
    params.manifest_dir = fullfile(params.processed_root, ...
        'postprocessing_manifests', params.output_folder_name);
end

if ~isfield(params, 'concat') || isempty(params.concat)
    params.concat = struct();
end
params.concat = local_set_default(params.concat, 'filename_pattern', '*_outputs_*.mat');
params.concat = local_set_default(params.concat, 'variable_name', 'EDMD_outputs');
params.concat = local_set_default(params.concat, 'concat_fields', {'efuns'});
params.concat = local_set_default(params.concat, 'scan_mode', 'uniform_except_last');
params.concat = local_set_default(params.concat, 'verbose', true);

if ~isfield(params, 'post_opts') || isempty(params.post_opts)
    params.post_opts = struct();
end
params.post_opts = local_set_default(params.post_opts, 'abs_thresh', 0.01);
params.post_opts = local_set_default(params.post_opts, 'sort_by', 'modulus');
params.post_opts = local_set_default(params.post_opts, 'sort_dir', 'descend');
params.post_opts = local_set_default(params.post_opts, 'max_basis', 80);
params.post_opts = local_set_default(params.post_opts, 'max_plot_samples', 2000);

if ~isfield(params, 'deconv') || isempty(params.deconv)
    params.deconv = struct();
end
params.deconv = local_set_default(params.deconv, 'method', 'koopman_residual');
params.deconv = local_set_default(params.deconv, 'lambda_source', 'edmd');
params.deconv = local_set_default(params.deconv, 'lambdaType', 'discrete');
params.deconv = local_set_default(params.deconv, 'max_modes_all', 80);
params.deconv = local_set_default(params.deconv, 'max_modes_sel', 40);
params.deconv = local_set_default(params.deconv, 'max_plot_samples', 2000);
params.deconv = local_set_default(params.deconv, 'plot_normalize_scope', 'global');

if ~isfield(params, 'timescale') || isempty(params.timescale)
    params.timescale = struct();
end
params.timescale = local_set_default(params.timescale, 'max_modes_all', 80);
params.timescale = local_set_default(params.timescale, 'max_modes_sel', 40);
params.timescale = local_set_default(params.timescale, 'maxLag', 200);
params.timescale = local_set_default(params.timescale, 'xlim_time', 240);

params = local_set_default(params, 'default_dt', 2);
params = local_set_default(params, 'make_main_plot', true);
params = local_set_default(params, 'make_deconv_plot', true);
params = local_set_default(params, 'draw_session_borders', true);
params = local_set_default(params, 'save_png', true);
params = local_set_default(params, 'save_fig', false);
params = local_set_default(params, 'resolution', 180);
params = local_set_default(params, 'skip_existing', true);
params = local_set_default(params, 'require_summary_file', true);
params = local_set_default(params, 'continue_on_error', true);
params = local_set_default(params, 'headless', true);
params = local_set_default(params, 'close_figures', true);
end


function runs = local_discover_runs(params)
runs = repmat(local_empty_run(), 0, 1);
for i_root = 1:numel(params.autodl_roots)
    root_dir = params.autodl_roots{i_root};
    if exist(root_dir, 'dir') ~= 7
        continue;
    end

    dataset_dirs = dir(root_dir);
    dataset_dirs = dataset_dirs([dataset_dirs.isdir]);
    dataset_dirs = dataset_dirs(~ismember({dataset_dirs.name}, {'.', '..'}));
    for i_ds = 1:numel(dataset_dirs)
        dataset_stem = dataset_dirs(i_ds).name;
        if ~isempty(params.dataset_stems) && ...
                ~any(strcmpi(dataset_stem, params.dataset_stems))
            continue;
        end

        output_parent = fullfile(root_dir, dataset_stem, 'mlp', 'outputs');
        if exist(output_parent, 'dir') ~= 7
            continue;
        end

        run_dirs = dir(output_parent);
        run_dirs = run_dirs([run_dirs.isdir]);
        run_dirs = run_dirs(~ismember({run_dirs.name}, {'.', '..'}));
        for i_run = 1:numel(run_dirs)
            run_name = run_dirs(i_run).name;
            output_dir = fullfile(run_dirs(i_run).folder, run_name);
            files = dir(fullfile(output_dir, '*_outputs_*.mat'));
            if isempty(files)
                continue;
            end
            if params.require_summary_file && isempty(dir(fullfile(output_dir, '*_summary.mat')))
                continue;
            end

            observable_mode = local_parse_observable_mode(run_name);
            residual_form = local_parse_residual_form(run_name);
            if ~isempty(params.observable_modes) && ...
                    ~any(strcmpi(observable_mode, params.observable_modes))
                continue;
            end
            if ~isempty(params.residual_forms) && ...
                    ~any(strcmpi(residual_form, params.residual_forms))
                continue;
            end

            run_info = local_empty_run();
            run_info.dataset_stem = dataset_stem;
            run_info.run_name = run_name;
            run_info.output_dir = output_dir;
            run_info.autodl_root = root_dir;
            run_info.observable_mode = observable_mode;
            run_info.residual_form = residual_form;
            run_info.n_output_files = numel(files);
            runs(end + 1, 1) = run_info; %#ok<AGROW>
        end
    end
end
end


function observable_file = local_find_observable_file(params, run_info, EDMD_outputs)
observable_file = '';
if isfield(EDMD_outputs, 'observable_file') && ~isempty(EDMD_outputs.observable_file)
    observable_file = char(string(EDMD_outputs.observable_file));
    if exist(observable_file, 'file') == 2
        return;
    end
end

mode = run_info.observable_mode;
if isempty(mode)
    mode = local_get_field(EDMD_outputs, 'observable_mode', '');
end
if isempty(mode)
    return;
end

dataset_dir = local_dataset_dir_name(run_info.dataset_stem);
candidate = fullfile(params.processed_root, 'bold_observables', dataset_dir, ...
    sprintf('%s_bold_observables_%s.mat', dataset_dir, mode));
if exist(candidate, 'file') == 2
    observable_file = candidate;
end
end


function obs_meta = local_load_observable_metadata(observable_file)
obs_meta = struct();
obs_meta.source_file = observable_file;
if isempty(observable_file) || exist(observable_file, 'file') ~= 2
    return;
end

wanted_names = {'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx', ...
    'dx', 'dt', 'fs', 'obs_info', 'observable_info', ...
    'observable_labels', 'file_stem', 'dataset_id', 'O', 'cfg', 'params'};
available = whos('-file', observable_file);
available_names = {available.name};
load_names = intersect(wanted_names, available_names, 'stable');
S = load(observable_file, load_names{:});

for i = 1:numel(load_names)
    obs_meta.(load_names{i}) = S.(load_names{i});
end

if isfield(S, 'obs_info') && ~isfield(obs_meta, 'observable_info')
    obs_meta.observable_info = S.obs_info;
end
if isfield(obs_meta, 'observable_info') && istable(obs_meta.observable_info) && ...
        ismember('observable_label', obs_meta.observable_info.Properties.VariableNames) && ...
        ~isfield(obs_meta, 'observable_labels')
    obs_meta.observable_labels = obs_meta.observable_info.observable_label;
end

if isfield(S, 'O') && isstruct(S.O)
    o_fields = {'session_ids', 'session_lengths', 'session_dx', ...
        'session_start_idx', 'session_end_idx', 'border_idx', ...
        'dx', 'fs', 'file_stem', 'dataset_id', ...
        'observable_info', 'observable_labels'};
    for i = 1:numel(o_fields)
        name = o_fields{i};
        if (~isfield(obs_meta, name) || isempty(obs_meta.(name))) && ...
                isfield(S.O, name) && ~isempty(S.O.(name))
            obs_meta.(name) = S.O.(name);
        end
    end
end

if isfield(S, 'cfg') && isstruct(S.cfg)
    if (~isfield(obs_meta, 'file_stem') || isempty(obs_meta.file_stem)) && ...
            isfield(S.cfg, 'file_stem') && ~isempty(S.cfg.file_stem)
        obs_meta.file_stem = S.cfg.file_stem;
    end
    if (~isfield(obs_meta, 'dataset_id') || isempty(obs_meta.dataset_id)) && ...
            isfield(S.cfg, 'dataset_id') && ~isempty(S.cfg.dataset_id)
        obs_meta.dataset_id = S.cfg.dataset_id;
    end
end
end


function session = local_resolve_session_metadata(EDMD_outputs, obs_meta, concat_info)
session = struct();
fields = {'session_ids', 'session_lengths', 'session_dx', ...
    'session_start_idx', 'session_end_idx', 'border_idx'};
for i = 1:numel(fields)
    name = fields{i};
    if isfield(EDMD_outputs, name) && ~isempty(EDMD_outputs.(name))
        session.(name) = EDMD_outputs.(name);
    elseif isfield(obs_meta, name) && ~isempty(obs_meta.(name))
        session.(name) = obs_meta.(name);
    else
        session.(name) = [];
    end
end

T = size(EDMD_outputs.efuns, 1);
if isempty(session.session_start_idx) || isempty(session.session_end_idx)
    session.session_start_idx = 1;
    session.session_end_idx = T;
    session.session_lengths = T;
    session.session_ids = 1;
    session.border_idx = [];
    session.source = 'single_trace_fallback';
else
    session.source = 'metadata';
end
session.session_start_idx = double(session.session_start_idx(:));
session.session_end_idx = double(session.session_end_idx(:));
if isempty(session.session_lengths)
    session.session_lengths = session.session_end_idx - session.session_start_idx + 1;
end
if isempty(session.session_ids)
    session.session_ids = (1:numel(session.session_start_idx)).';
end
if isempty(session.border_idx)
    session.border_idx = session.session_end_idx(1:end-1);
end
if ~isempty(concat_info)
    session.concat_total_length = local_get_field(concat_info, 'total_length', T);
end
end


function dt = local_resolve_dt(EDMD_outputs, obs_meta, default_dt)
dt = [];
candidates = {'dt', 'dx', 'sampling_period', 'sample_period'};
for i = 1:numel(candidates)
    name = candidates{i};
    if isfield(EDMD_outputs, name) && ~isempty(EDMD_outputs.(name))
        dt = double(EDMD_outputs.(name));
        break;
    end
    if isfield(obs_meta, name) && ~isempty(obs_meta.(name))
        dt = double(obs_meta.(name));
        break;
    end
end
if isempty(dt) && isfield(EDMD_outputs, 'fs') && ~isempty(EDMD_outputs.fs)
    dt = 1 / double(EDMD_outputs.fs);
end
if isempty(dt) && isfield(obs_meta, 'fs') && ~isempty(obs_meta.fs)
    dt = 1 / double(obs_meta.fs);
end
if isempty(dt) || ~isfinite(dt(1)) || dt(1) <= 0
    dt = default_dt;
else
    dt = dt(1);
end
end


function [time_vec, session_border_t] = local_time_axis_and_borders(T, dt, session)
time_vec = (0:T-1)' * dt;
border_idx = double(session.border_idx(:));
border_idx = border_idx(border_idx >= 1 & border_idx <= T);
session_border_t = time_vec(border_idx);
end


function path_out = local_save_figure(fig, fig_dir, name, params)
path_out = '';
if isempty(fig) || ~isvalid(fig)
    return;
end
if params.save_png
    path_out = fullfile(fig_dir, name);
    exportgraphics(fig, path_out, 'Resolution', params.resolution);
end
if params.save_fig
    [~, stem] = fileparts(name);
    savefig(fig, fullfile(fig_dir, [stem, '.fig']));
end
end


function S = local_set_default(S, name, value)
if ~isfield(S, name) || isempty(S.(name))
    S.(name) = value;
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function run_info = local_empty_run()
run_info = struct('dataset_stem', '', 'run_name', '', 'output_dir', '', ...
    'autodl_root', '', 'observable_mode', '', 'residual_form', '', ...
    'n_output_files', 0);
end


function row = local_empty_row()
row = struct('dataset_stem', '', 'run_name', '', 'observable_mode', '', ...
    'residual_form', '', 'status', '', 'message', '', 'post_file', '', ...
    'main_png', '', 'deconv_png', '', 'timescale_png', '', ...
    'runtime_sec', NaN);
end


function row = local_make_row(run_info, status, message, post_file, ...
    main_png, deconv_png, timescale_png, runtime_sec)
row = local_empty_row();
row.dataset_stem = run_info.dataset_stem;
row.run_name = run_info.run_name;
row.observable_mode = run_info.observable_mode;
row.residual_form = run_info.residual_form;
row.status = status;
row.message = message;
row.post_file = post_file;
row.main_png = main_png;
row.deconv_png = deconv_png;
row.timescale_png = timescale_png;
row.runtime_sec = runtime_sec;
end


function mode = local_parse_observable_mode(run_name)
known = {'global_svd100', 'HP_svd100', 'slow_band_power', ...
    'roi_mean', 'eleHP', 'HP', 'svd'};
mode = '';
for i = 1:numel(known)
    pat = ['_', known{i}];
    if contains(run_name, pat, 'IgnoreCase', true)
        mode = known{i};
        return;
    end
end
end


function form = local_parse_residual_form(run_name)
forms = {'projected_vlambda', 'projected_kv'};
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
        dataset_dir = 'F12.M01';
    otherwise
        dataset_dir = dataset_stem;
end
end
