function sync_pipeline6_dataset_window_figures(cfg, window_cache, run_table, output_root)
%SYNC_PIPELINE6_DATASET_WINDOW_FIGURES Copy comparable window figures into one folder tree.

if nargin < 4 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if isempty(run_table) || height(run_table) == 0
    return;
end

raw_manifest = local_load_pipeline2_window_manifest(cfg, output_root);
dataset_dir = io_project.get_dataset_window_figures_dir(output_root, cfg);
if exist(dataset_dir, 'dir') ~= 7
    mkdir(dataset_dir);
end

for i = 1:height(run_table)
    row = run_table(i, :);
    top_dir = fullfile(dataset_dir, local_window_folder_name(row));
    if exist(top_dir, 'dir') ~= 7
        mkdir(top_dir);
    end

    local_copy_pipeline2_plot(raw_manifest, row, top_dir);
    local_copy_pipeline6_plot(row.postprocess_main_png, top_dir, ...
        local_pipeline6_file_name(row, 'postprocess_main'));
    local_copy_pipeline6_plot(row.deconv_png, top_dir, ...
        local_pipeline6_file_name(row, 'deconv'));
    local_copy_pipeline6_plot(row.deconv_localwin_png, top_dir, ...
        local_pipeline6_file_name(row, 'deconv_localwin_norm'));
end
end


function manifest = local_load_pipeline2_window_manifest(cfg, output_root)
manifest_file = fullfile( ...
    io_project.get_pipeline_stage_dir(output_root, cfg, 2, 'figures_consensus_state_top_window_plots'), ...
    'plot_manifest.csv');

if exist(manifest_file, 'file') ~= 2
    error('Pipeline 2 top-window plot manifest is missing: %s', manifest_file);
end

manifest = readtable(manifest_file, 'TextType', 'string');
manifest.Properties.UserData.manifest_file = manifest_file;
end


function folder_name = local_window_folder_name(row)
rank_val = double(local_table_value(row, 'window_rank'));
start_idx = double(local_table_value(row, 'global_start_idx'));
end_idx = double(local_table_value(row, 'global_end_idx'));

folder_name = sprintf('top_win_rank_%02d_%d_%d', rank_val, start_idx, end_idx);
end


function local_copy_pipeline2_plot(manifest, row, top_dir)
rank_val = double(local_table_value(row, 'window_rank'));
start_idx = double(local_table_value(row, 'global_start_idx'));
end_idx = double(local_table_value(row, 'global_end_idx'));

mask = double(manifest.state_diversity_rank) == rank_val & ...
    double(manifest.global_start_idx) == start_idx & ...
    double(manifest.global_end_idx) == end_idx;
if ~any(mask)
    return;
end

src = char(manifest.plot_file(find(mask, 1, 'first')));
src = local_resolve_existing_plot_path(src, manifest.Properties.UserData.manifest_file);
if isempty(src)
    return;
end

dst = fullfile(top_dir, 'pipeline2_consensus_state_window.png');
copyfile(src, dst, 'f');
end


function local_copy_pipeline6_plot(src_in, top_dir, dst_name)
src = char(string(src_in));
if isempty(src)
    return;
end

src = local_resolve_existing_plot_path(src, '');
if isempty(src)
    return;
end

dst = fullfile(top_dir, dst_name);
copyfile(src, dst, 'f');
end


function dst_name = local_pipeline6_file_name(row, suffix)
obs = char(local_table_value(row, 'observable_mode'));
resid = char(local_table_value(row, 'residual_form'));
run_name = char(local_table_value(row, 'run_name'));
obs = matlab.lang.makeValidName(obs);
resid = matlab.lang.makeValidName(resid);
run_name = local_filename_safe(run_name);
dst_name = sprintf('pipeline6_%s_%s_%s_%s.png', run_name, obs, resid, suffix);
end


function value = local_table_value(row, var_name)
if ismember(var_name, row.Properties.VariableNames)
    value = row.(var_name)(1);
else
    error('Window row is missing variable: %s', var_name);
end
end


function out = local_filename_safe(name_in)
out = regexprep(char(string(name_in)), '[^a-zA-Z0-9_\-]+', '_');
out = regexprep(out, '_+', '_');
out = regexprep(out, '^_+|_+$', '');
if isempty(out)
    out = 'untitled_run';
end
end


function resolved = local_resolve_existing_plot_path(path_in, manifest_file)
resolved = '';
path_in = char(string(path_in));
if isempty(path_in)
    return;
end

if exist(path_in, 'file') == 2
    resolved = path_in;
    return;
end

[~, name, ext] = fileparts(path_in);
if ~isempty(manifest_file)
    manifest_dir = fileparts(manifest_file);
    candidate = fullfile(manifest_dir, [name ext]);
    if exist(candidate, 'file') == 2
        resolved = candidate;
        return;
    end
end

path_parts = strsplit(path_in, filesep);
stage_candidates = {'01_postprocess_main', '03_deconv', '04_deconv_localwin_norm'};
for i = 1:numel(stage_candidates)
    idx = find(strcmpi(path_parts, stage_candidates{i}), 1, 'last');
    if isempty(idx) || idx < 2
        continue;
    end
    run_dir = fullfile(path_parts{1:idx-1});
    candidate = fullfile(run_dir, stage_candidates{i}, [name ext]);
    if exist(candidate, 'file') == 2
        resolved = candidate;
        return;
    end
end
end
