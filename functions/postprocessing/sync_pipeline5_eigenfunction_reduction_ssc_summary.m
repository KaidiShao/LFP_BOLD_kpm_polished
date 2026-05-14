function sync_pipeline5_eigenfunction_reduction_ssc_summary( ...
        cfg, run_info, method_tag, plot_paths, output_root)
%SYNC_PIPELINE5_EIGENFUNCTION_REDUCTION_SSC_SUMMARY Copy run-level figures into shared summary folders.

if nargin < 5 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if nargin < 4 || isempty(plot_paths)
    return;
end

condition_tag = local_resolve_condition_tag(cfg, run_info);
dataset_stem = local_resolve_dataset_stem(cfg, run_info);
summary_root = fullfile(output_root, dataset_stem, ...
    'pipeline5_summary_figures', condition_tag, 'eigenfunction_reduction');

local_copy_plot(plot_paths, 'state_space_png', ...
    fullfile(summary_root, 'state_space'), method_tag);
local_copy_plot(plot_paths, 'consensus_state_space_png', ...
    fullfile(summary_root, 'consensus_state_space'), method_tag);
local_copy_plot(plot_paths, 'spectrum_diagnostics_png', ...
    fullfile(summary_root, 'spectrum_diagnostics'), method_tag);
end


function local_copy_plot(plot_paths, field_name, dst_dir, method_tag)
if ~isfield(plot_paths, field_name)
    return;
end

src = char(string(plot_paths.(field_name)));
if isempty(src) || exist(src, 'file') ~= 2
    return;
end

if exist(dst_dir, 'dir') ~= 7
    mkdir(dst_dir);
end

old_files = dir(fullfile(dst_dir, sprintf('%s__*.png', char(string(method_tag)))));
for i = 1:numel(old_files)
    delete(fullfile(old_files(i).folder, old_files(i).name));
end

[~, name, ext] = fileparts(src);
dst = fullfile(dst_dir, sprintf('%s__%s%s', char(string(method_tag)), name, ext));
copyfile(src, dst, 'f');
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
    condition_tag = build_blp_eigenfunction_condition_tag(run_info);
end
end


function dataset_stem = local_resolve_dataset_stem(cfg, run_info)
dataset_stem = '';
if isstruct(cfg)
    if isfield(cfg, 'file_stem') && ~isempty(cfg.file_stem)
        dataset_stem = char(string(cfg.file_stem));
    elseif isfield(cfg, 'dataset') && isstruct(cfg.dataset) && ...
            isfield(cfg.dataset, 'name') && ~isempty(cfg.dataset.name)
        dataset_stem = char(string(cfg.dataset.name));
    end
end

if isempty(dataset_stem) && isstruct(run_info) && ...
        isfield(run_info, 'dataset') && ~isempty(run_info.dataset)
    dataset_stem = char(string(run_info.dataset));
end

if isempty(dataset_stem)
    error('Unable to resolve dataset stem for pipeline 5 summary sync.');
end
end
