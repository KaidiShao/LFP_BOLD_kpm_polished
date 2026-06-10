% Redraw current Pipeline 5 run-level figures from existing reduction MAT files.
%
% This helper does not recompute eigenfunction dimension reduction. It reloads
% canonical P5 reduction MAT files, redraws the current scatter renderers, and
% refreshes the per-dataset pipeline5_summary_figures cache.
%
% Optional variables before run(...):
%   dataset_stems       = {'e10gb1','e10fV1','e10gh1','f12m01','e10gw1'};
%   condition_tags      = {'abs_projected_vlambda','complex_split_projected_vlambda'};
%   method_filter       = {'svd','nmf','mds','umap'};
%   component_counts    = 3:8;
%   include_state_space = true;
%   include_consensus   = true;
%   include_spectrum    = true;
%   continue_on_error   = true;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
set(groot, 'defaultFigureVisible', 'off');

if ~exist('dataset_stems', 'var') || isempty(dataset_stems)
    dataset_stems = {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01', 'e10gw1'};
end
if ~exist('condition_tags', 'var') || isempty(condition_tags)
    condition_tags = {'abs_projected_vlambda', 'complex_split_projected_vlambda'};
end
if ~exist('method_filter', 'var') || isempty(method_filter)
    method_filter = {'svd', 'nmf', 'mds', 'umap'};
end
if ~exist('component_counts', 'var') || isempty(component_counts)
    component_counts = 3:8;
end
if ~exist('include_state_space', 'var') || isempty(include_state_space)
    include_state_space = true;
end
if ~exist('include_consensus', 'var') || isempty(include_consensus)
    include_consensus = true;
end
if ~exist('include_spectrum', 'var') || isempty(include_spectrum)
    include_spectrum = true;
end
if ~exist('continue_on_error', 'var') || isempty(continue_on_error)
    continue_on_error = true;
end

dataset_stems = cellstr(string(dataset_stems(:)).');
condition_tags = cellstr(string(condition_tags(:)).');
method_filter = cellstr(lower(string(method_filter(:))).');
component_counts = double(component_counts(:).');
include_state_space = logical(include_state_space);
include_consensus = logical(include_consensus);
include_spectrum = logical(include_spectrum);
continue_on_error = logical(continue_on_error);

processed_root = io_project.get_project_processed_root();
stage_name = io_project.get_pipeline_stage_name(5, 'eigenfunction_reduction');
manifest_dir = fullfile(processed_root, 'postprocessing_manifests', ...
    'pipeline5_redraw_current_scatter');
if exist(manifest_dir, 'dir') ~= 7
    mkdir(manifest_dir);
end
manifest_file = fullfile(manifest_dir, sprintf( ...
    'pipeline5_redraw_current_scatter_%s.csv', datestr(now, 'yyyymmdd_HHMMSS')));

fprintf('Redrawing current Pipeline 5 figures from MAT files.\n');
fprintf('Processed root: %s\n', processed_root);
fprintf('Datasets      : %s\n', strjoin(dataset_stems, ', '));
fprintf('Conditions    : %s\n', strjoin(condition_tags, ', '));
fprintf('Methods       : %s\n', strjoin(method_filter, ', '));
fprintf('k             : %s\n\n', strjoin(cellstr(string(component_counts)), ', '));

rows = struct([]);
row_count = 0;
consensus_cache = struct();

for i_ds = 1:numel(dataset_stems)
    dataset = dataset_stems{i_ds};
    for i_cond = 1:numel(condition_tags)
        condition_tag = condition_tags{i_cond};
        for i_k = 1:numel(component_counts)
            k = component_counts(i_k);
            for i_method = 1:numel(method_filter)
                method = method_filter{i_method};
                method_tag = sprintf('%s_k%02d', method, k);
                method_root = fullfile(processed_root, dataset, stage_name, ...
                    condition_tag, method_tag);
                result_file = local_latest_reduction_result_file(method_root);

                fprintf('[%s | %s | %s]\n', dataset, condition_tag, method_tag);
                if isempty(result_file)
                    fprintf('  missing source MAT\n\n');
                    row_count = row_count + 1;
                    rows(row_count).dataset = string(dataset); %#ok<SAGROW>
                    rows(row_count).condition = string(condition_tag);
                    rows(row_count).method_tag = string(method_tag);
                    rows(row_count).status = "missing_source";
                    rows(row_count).message = "";
                    rows(row_count).result_file = "";
                    rows(row_count).state_space_png = "";
                    rows(row_count).consensus_state_space_png = "";
                    rows(row_count).spectrum_diagnostics_png = "";
                    continue;
                end

                try
                    S = load_mat_file_with_short_path(result_file, 'result');
                    if ~isfield(S, 'result')
                        error('MAT file does not contain variable result.');
                    end
                    result = S.result;
                    if ~isfield(result, 'cfg') || isempty(result.cfg)
                        error('Result has no embedded cfg.');
                    end

                    cfg = local_refresh_plot_cfg(result.cfg, processed_root, ...
                        dataset, stage_name, condition_tag, method_tag, ...
                        include_state_space, include_consensus, include_spectrum);

                    [C_consensus, consensus_cache] = local_consensus_for_dataset( ...
                        dataset, include_consensus, processed_root, consensus_cache);

                    params = struct();
                    params.close_figures = true;
                    plot_paths = plot_blp_eigenfunction_reduction_outputs( ...
                        result, cfg, params, C_consensus);

                    run_info = struct();
                    run_info.dataset = dataset;
                    run_info.run_name = local_get_nested(cfg, ...
                        {'output', 'source_run_name'}, '');
                    sync_pipeline5_eigenfunction_reduction_ssc_summary( ...
                        cfg, run_info, method_tag, plot_paths, processed_root);

                    fprintf('  OK\n\n');
                    row_count = row_count + 1;
                    rows(row_count).dataset = string(dataset); %#ok<SAGROW>
                    rows(row_count).condition = string(condition_tag);
                    rows(row_count).method_tag = string(method_tag);
                    rows(row_count).status = "ok";
                    rows(row_count).message = "";
                    rows(row_count).result_file = string(result_file);
                    rows(row_count).state_space_png = local_path_field(plot_paths, 'state_space_png');
                    rows(row_count).consensus_state_space_png = local_path_field(plot_paths, 'consensus_state_space_png');
                    rows(row_count).spectrum_diagnostics_png = local_path_field(plot_paths, 'spectrum_diagnostics_png');
                catch ME
                    fprintf(2, '  FAILED: %s\n\n', ME.message);
                    row_count = row_count + 1;
                    rows(row_count).dataset = string(dataset); %#ok<SAGROW>
                    rows(row_count).condition = string(condition_tag);
                    rows(row_count).method_tag = string(method_tag);
                    rows(row_count).status = "error";
                    rows(row_count).message = string(ME.message);
                    rows(row_count).result_file = string(result_file);
                    rows(row_count).state_space_png = "";
                    rows(row_count).consensus_state_space_png = "";
                    rows(row_count).spectrum_diagnostics_png = "";
                    if ~continue_on_error
                        rethrow(ME);
                    end
                end
            end
        end
    end
end

if ~isempty(rows)
    writetable(struct2table(rows), manifest_file);
end

fprintf('Redraw manifest:\n  %s\n', manifest_file);


function result_file = local_latest_reduction_result_file(method_root)
result_file = '';
mat_dir = fullfile(method_root, 'mat');
if exist(mat_dir, 'dir') ~= 7
    return;
end

L = dir(fullfile(mat_dir, '*.mat'));
if isempty(L)
    return;
end
names = string({L.name});
keep = contains(names, "_efun__time__", 'IgnoreCase', true) | ...
    contains(names, "_efun__spectrum__", 'IgnoreCase', true);
L = L(keep);
if isempty(L)
    return;
end
[~, idx] = max([L.datenum]);
result_file = fullfile(L(idx).folder, L(idx).name);
end


function cfg = local_refresh_plot_cfg(cfg, processed_root, dataset, stage_name, ...
        condition_tag, method_tag, include_state_space, include_consensus, include_spectrum)
method_root = fullfile(processed_root, dataset, stage_name, condition_tag, method_tag);
figure_dir = fullfile(method_root, 'fig');
if exist(figure_dir, 'dir') ~= 7
    mkdir(figure_dir);
end

cfg.output.root = method_root;
cfg.output.figure_dir = figure_dir;
cfg.output.condition_tag = condition_tag;
cfg.output.output_run_name = condition_tag;
cfg.dataset.name = dataset;
cfg.dataset.processed_root = processed_root;

if ~isfield(cfg, 'viz') || isempty(cfg.viz)
    defaults = build_blp_eigenfunction_reduction_defaults();
    cfg.viz = defaults.viz;
end

defaults = build_blp_eigenfunction_reduction_defaults();
if ~isfield(cfg.viz, 'state_space') || isempty(cfg.viz.state_space)
    cfg.viz.state_space = defaults.viz.state_space;
end
if ~isfield(cfg.viz, 'state_space_consensus') || isempty(cfg.viz.state_space_consensus)
    cfg.viz.state_space_consensus = defaults.viz.state_space_consensus;
end
if ~isfield(cfg.viz, 'spectrum') || isempty(cfg.viz.spectrum)
    cfg.viz.spectrum = defaults.viz.spectrum;
end

cfg.viz.state_space.enable = include_state_space;
cfg.viz.state_space.save_dir = figure_dir;
cfg.viz.state_space.save_tag = [method_tag '_ss'];
cfg.viz.state_space.save_figure = true;
cfg.viz.state_space.render_mode = 'scatter';

cfg.viz.state_space_consensus.enable = include_consensus;
cfg.viz.state_space_consensus.save_dir = figure_dir;
cfg.viz.state_space_consensus.save_tag = [method_tag '_ssc'];
cfg.viz.state_space_consensus.save_figure = true;
cfg.viz.state_space_consensus.render_mode = 'scatter';

cfg.viz.spectrum.enable = include_spectrum;
cfg.viz.spectrum.save_dir = figure_dir;
cfg.viz.spectrum.save_tag = [method_tag '_spec'];
cfg.viz.spectrum.save_figure = true;
end


function [C_consensus, cache] = local_consensus_for_dataset( ...
        dataset, include_consensus, processed_root, cache)
C_consensus = [];
if ~include_consensus
    return;
end

field_name = matlab.lang.makeValidName(dataset);
if isfield(cache, field_name)
    C_consensus = cache.(field_name);
    return;
end

loader_cfg = struct();
loader_cfg.file_stem = dataset;
[C_consensus, source_consensus_file] = io_results.load_consensus_state_results( ...
    loader_cfg, processed_root, []);
fprintf('  consensus: %s\n', source_consensus_file);
cache.(field_name) = C_consensus;
end


function value = local_path_field(S, field_name)
value = "";
if isstruct(S) && isfield(S, field_name)
    value = string(S.(field_name));
end
end


function value = local_get_nested(S, names, default_value)
value = default_value;
current = S;
for i = 1:numel(names)
    name = names{i};
    if ~isstruct(current) || ~isfield(current, name)
        return;
    end
    current = current.(name);
end
value = current;
end
