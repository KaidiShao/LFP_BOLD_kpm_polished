function sync_pipeline5_dataset_window_figures( ...
        cfg, manifest_file, condition_tag, method_tag, output_root)
%SYNC_PIPELINE5_DATASET_WINDOW_FIGURES Copy dimred top-window figures into shared window folders.

if nargin < 5 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

if isempty(manifest_file) || exist(manifest_file, 'file') ~= 2
    error('Top30 manifest file was not found: %s', char(string(manifest_file)));
end

manifest = readtable(manifest_file, 'TextType', 'string');
dataset_dir = io_project.get_dataset_window_figures_dir(output_root, cfg);
if exist(dataset_dir, 'dir') ~= 7
    mkdir(dataset_dir);
end

for i = 1:height(manifest)
    row = manifest(i, :);
    top_dir = fullfile(dataset_dir, local_window_folder_name(row));
    if exist(top_dir, 'dir') ~= 7
        mkdir(top_dir);
    end

    local_copy_plot(row, manifest_file, 'eigenfunction_overview_png', ...
        fullfile(top_dir, sprintf('pipeline5_%s_%s_overview.png', condition_tag, method_tag)));

    if ismember('state_space_png', row.Properties.VariableNames)
        local_copy_plot(row, manifest_file, 'state_space_png', ...
            fullfile(top_dir, sprintf('pipeline5_%s_%s_state_space.png', condition_tag, method_tag)));
    end

    local_copy_plot(row, manifest_file, 'consensus_state_space_png', ...
        fullfile(top_dir, sprintf('pipeline5_%s_%s_consensus_state_space.png', condition_tag, method_tag)));
end
end


function folder_name = local_window_folder_name(row)
rank_val = double(local_row_value(row, 'state_diversity_rank'));
start_idx = double(local_row_value(row, 'global_start_idx'));
end_idx = double(local_row_value(row, 'global_end_idx'));
folder_name = sprintf('top_win_rank_%02d_%d_%d', rank_val, start_idx, end_idx);
end


function value = local_row_value(row, var_name)
if ~ismember(var_name, row.Properties.VariableNames)
    error('Top30 manifest row is missing variable: %s', var_name);
end
value = row.(var_name)(1);
end


function local_copy_plot(row, manifest_file, var_name, dst_file)
if ~ismember(var_name, row.Properties.VariableNames)
    return;
end

src = char(string(row.(var_name)(1)));
src = local_resolve_existing_plot_path(src, manifest_file);
if isempty(src)
    return;
end

copyfile(src, dst_file, 'f');
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
manifest_dir = fileparts(manifest_file);
candidate = fullfile(manifest_dir, [name ext]);
if exist(candidate, 'file') == 2
    resolved = candidate;
end
end
