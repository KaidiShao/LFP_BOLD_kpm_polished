function ctx = prepare_bold_reskoopnet_postprocessing_run(run_info, params)
%PREPARE_BOLD_RESKOOPNET_POSTPROCESSING_RUN Resolve shared inputs for one pipeline 7 run.

if nargin < 1 || ~isstruct(run_info)
    error('run_info must be a struct.');
end
if nargin < 2 || ~isstruct(params)
    error('params must be a struct.');
end

run_info = local_normalize_run_info(run_info);

ctx = struct();
ctx.run_info = run_info;
ctx.out_root = fullfile( ...
    io_project.get_pipeline_stage_dir(params.processed_root, ...
    run_info.dataset_stem, 7, 'bold_postprocessing'), ...
    run_info.run_name);
ctx.mat_dir = fullfile(ctx.out_root, 'mat');
ctx.fig_dir = fullfile(ctx.out_root, 'fig');
if exist(ctx.mat_dir, 'dir') ~= 7, mkdir(ctx.mat_dir); end
if exist(ctx.fig_dir, 'dir') ~= 7, mkdir(ctx.fig_dir); end

ctx.post_file = fullfile(ctx.mat_dir, [run_info.run_name, '_bold_post.mat']);
ctx.skip_existing = isfield(params, 'skip_existing') && logical(params.skip_existing) && ...
    local_is_valid_bold_post_file(ctx.post_file);
if ctx.skip_existing
    return;
end

source_cfg = struct();
source_cfg.mode = 'chunk_dir';
source_cfg.data_dir = run_info.output_dir;
source_cfg.concat = params.concat;
[ctx.EDMD_outputs, ctx.concat_info, ctx.source_info] = io_edmd.load_edmd_source(source_cfg);

ctx.observable_file = local_find_observable_file(params, run_info, ctx.EDMD_outputs);
ctx.obs_meta = local_load_observable_metadata(ctx.observable_file);
ctx.session = local_resolve_session_metadata( ...
    ctx.EDMD_outputs, ctx.obs_meta, ctx.concat_info);
ctx.dt = local_resolve_dt(ctx.EDMD_outputs, ctx.obs_meta, params.default_dt);
[ctx.time_vec, ctx.session_border_t] = local_time_axis_and_borders( ...
    size(ctx.EDMD_outputs.efuns, 1), ctx.dt, ctx.session);
end


function run_info = local_normalize_run_info(run_info)
required = {'dataset_stem', 'run_name', 'output_dir'};
for i = 1:numel(required)
    name = required{i};
    if ~isfield(run_info, name) || isempty(run_info.(name))
        error('run_info.%s is required.', name);
    end
    run_info.(name) = char(string(run_info.(name)));
end
if exist(run_info.output_dir, 'dir') ~= 7
    error('run_info.output_dir does not exist: %s', run_info.output_dir);
end
if ~isfield(run_info, 'dataset_id') || isempty(run_info.dataset_id)
    run_info.dataset_id = io_project.get_dataset_id_from_stem(run_info.dataset_stem);
else
    run_info.dataset_id = char(string(run_info.dataset_id));
end
if ~isfield(run_info, 'observable_mode') || isempty(run_info.observable_mode)
    run_info.observable_mode = local_parse_observable_mode(run_info.run_name);
end
if ~isfield(run_info, 'residual_form') || isempty(run_info.residual_form)
    run_info.residual_form = local_parse_residual_form(run_info.run_name);
end
end


function observable_file = local_find_observable_file(params, run_info, EDMD_outputs)
observable_file = '';

if isfield(run_info, 'observable_file') && ~isempty(run_info.observable_file)
    observable_file = char(string(run_info.observable_file));
    if exist(observable_file, 'file') == 2
        return;
    end
end

if isfield(EDMD_outputs, 'observable_file') && ~isempty(EDMD_outputs.observable_file)
    observable_file = char(string(EDMD_outputs.observable_file));
    if exist(observable_file, 'file') == 2
        return;
    end
end

mode = run_info.observable_mode;
if isempty(mode) && isfield(EDMD_outputs, 'observable_mode')
    mode = char(string(EDMD_outputs.observable_mode));
end
if isempty(mode)
    return;
end

obs_dir = io_project.get_pipeline_stage_dir( ...
    params.processed_root, run_info.dataset_stem, 3, 'bold_observables');
candidate = fullfile(obs_dir, ...
    sprintf('%s_bold_observables_%s.mat', run_info.dataset_id, mode));
if exist(candidate, 'file') == 2
    observable_file = candidate;
    return;
end

L = dir(fullfile(obs_dir, sprintf('*_bold_observables_%s.mat', mode)));
if isscalar(L)
    observable_file = fullfile(L(1).folder, L(1).name);
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
    'observable_labels', 'file_stem', 'dataset_id', 'O', 'cfg', 'cfg_saved', 'params'};
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
if (~isfield(obs_meta, 'cfg') || isempty(obs_meta.cfg)) && ...
        isfield(S, 'cfg_saved') && isstruct(S.cfg_saved)
    obs_meta.cfg = S.cfg_saved;
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
if ~isempty(concat_info) && isstruct(concat_info) && ...
        isfield(concat_info, 'total_length') && ~isempty(concat_info.total_length)
    session.concat_total_length = concat_info.total_length;
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


function tf = local_is_valid_bold_post_file(post_file)
tf = false;
if exist(post_file, 'file') ~= 2
    return;
end
try
    vars = whos('-file', post_file);
catch
    return;
end
tf = any(strcmp({vars.name}, 'BOLD_POST'));
end


function mode = local_parse_observable_mode(run_name)
known = {'global_slow_band_power_svd100', 'roi_mean_slow_band_power', ...
    'slow_band_power_svd', 'gsvd100_ds', 'global_svd100', 'HP_svd100', 'slow_band_power', ...
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
