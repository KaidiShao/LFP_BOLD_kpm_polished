function manifest_file = run_blp_eigenfunction_dimred_top30_plots( ...
        cfg, result_file, run_info, method_tag, params, repo_root)
%RUN_BLP_EIGENFUNCTION_DIMRED_TOP30_PLOTS Export pipeline 5 dimred top-window figures.

if nargin < 5 || isempty(params)
    params = struct();
end

if nargin < 6 || isempty(repo_root)
    repo_root = fileparts(fileparts(mfilename('fullpath')));
end

top30_script = fullfile(repo_root, 'scripts', ...
    'script_plot_e10gb1_efun_svd_top30_state_diversity_windows.m');
if exist(top30_script, 'file') ~= 2
    error('Top30 plotting script was not found: %s', top30_script);
end
if isempty(result_file) || exist(result_file, 'file') ~= 2
    error('Result MAT file for top30 plotting was not found: %s', result_file);
end

condition_tag = local_resolve_condition_tag(cfg, run_info);
method_tag = local_resolve_method_tag(cfg, method_tag);
dataset_name = local_resolve_dataset_name(cfg, run_info);
wrapper_params = params;
processed_root = local_resolve_processed_root(cfg, wrapper_params);
stage_root = io_project.get_pipeline_stage_dir( ...
    processed_root, dataset_name, 5, 'efun_dimred_top30');
save_dir = fullfile(stage_root, condition_tag, method_tag);

cfg = cfg; %#ok<NASGU>
result_file = result_file; %#ok<NASGU>
top_window_params = struct(); %#ok<NASGU>
top_window_params.n_top_windows = local_get_param(wrapper_params, 'top30_n_windows', 30);
top_window_params.make_overview = true;
top_window_params.make_state_space = false;
top_window_params.make_consensus_state_space = true;
top_window_params.close_figures = local_get_param(wrapper_params, 'close_figures', true);
top_window_params.force_recompute_norm_cache = ...
    local_get_param(wrapper_params, 'top30_force_recompute_norm_cache', false);
top_window_params.save_dir = save_dir;
top_window_params.norm_cache_file = fullfile(save_dir, 'norm.mat');

manifest_file = '';
run(top30_script);

if exist('manifest_file', 'var') ~= 1 || isempty(manifest_file) %#ok<EXIST>
    error('Top30 plotting script did not return manifest_file.');
end

sync_pipeline5_dataset_window_figures( ...
    dataset_name, manifest_file, condition_tag, method_tag, processed_root);
end


function value = local_get_param(params, field_name, default_value)
if isfield(params, field_name) && ~isempty(params.(field_name))
    value = params.(field_name);
else
    value = default_value;
end
end


function processed_root = local_resolve_processed_root(cfg, params)
processed_root = '';
if isstruct(params) && isfield(params, 'processed_root') && ~isempty(params.processed_root)
    processed_root = char(string(params.processed_root));
elseif isstruct(cfg) && isfield(cfg, 'dataset') && isstruct(cfg.dataset) && ...
        isfield(cfg.dataset, 'processed_root') && ~isempty(cfg.dataset.processed_root)
    processed_root = char(string(cfg.dataset.processed_root));
else
    processed_root = io_project.get_project_processed_root();
end
end


function condition_tag = local_resolve_condition_tag(cfg, run_info)
condition_tag = '';

if isstruct(cfg) && isfield(cfg, 'output') && isstruct(cfg.output)
    if isfield(cfg.output, 'condition_tag') && ~isempty(cfg.output.condition_tag)
        condition_tag = char(string(cfg.output.condition_tag));
    elseif isfield(cfg.output, 'output_run_name') && ~isempty(cfg.output.output_run_name)
        condition_tag = char(string(cfg.output.output_run_name));
    end
end

if isempty(condition_tag)
    path_hint = local_get_cfg_dir_hint(cfg);
    if ~isempty(path_hint)
        [method_dir, ~, ~] = fileparts(path_hint);
        [condition_dir, ~, ~] = fileparts(method_dir);
        [~, condition_tag] = fileparts(condition_dir);
    end
end

if isempty(condition_tag)
    if nargin >= 2 && isstruct(run_info) && ...
            isfield(run_info, 'observable_mode') && ~isempty(run_info.observable_mode) && ...
            isfield(run_info, 'residual_form') && ~isempty(run_info.residual_form)
        condition_tag = build_blp_eigenfunction_condition_tag(run_info);
    else
        condition_tag = 'unknown_condition';
    end
end
end


function method_tag = local_resolve_method_tag(cfg, method_tag_in)
method_tag = char(string(method_tag_in));
if ~isempty(method_tag)
    return;
end

path_hint = local_get_cfg_dir_hint(cfg);
if ~isempty(path_hint)
    [method_dir, ~, ~] = fileparts(path_hint);
    [~, method_tag] = fileparts(method_dir);
end

if isempty(method_tag)
    method_tag = 'unknown_method';
end
end


function path_hint = local_get_cfg_dir_hint(cfg)
path_hint = '';
if isfield(cfg, 'save') && isfield(cfg.save, 'dir') && ~isempty(cfg.save.dir)
    path_hint = char(string(cfg.save.dir));
elseif isfield(cfg, 'output') && isfield(cfg.output, 'root') && ~isempty(cfg.output.root)
    path_hint = char(string(cfg.output.root));
end
end


function dataset_name = local_resolve_dataset_name(cfg, run_info)
dataset_name = '';
if isstruct(cfg)
    if isfield(cfg, 'file_stem') && ~isempty(cfg.file_stem)
        dataset_name = char(string(cfg.file_stem));
    elseif isfield(cfg, 'dataset') && isstruct(cfg.dataset) && ...
            isfield(cfg.dataset, 'name') && ~isempty(cfg.dataset.name)
        dataset_name = char(string(cfg.dataset.name));
    end
end

if isempty(dataset_name) && nargin >= 2 && isstruct(run_info) && ...
        isfield(run_info, 'dataset') && ~isempty(run_info.dataset)
    dataset_name = char(string(run_info.dataset));
end

if isempty(dataset_name)
    error('Unable to resolve dataset name for pipeline 5 top-window plots.');
end
end
