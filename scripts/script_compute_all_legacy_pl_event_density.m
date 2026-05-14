% Batch-compute legacy PL event density for all matching DataPons datasets.
%
% Scans E:\DataPons recursively for exact legacy detector files:
%   <stem>_<group>_mevt_pl.mat
%   <stem>_<group>_sevt_pl.mat
%
% Outputs density .mat files under each processed dataset's:
%   pipeline2_legacy_event_density
%
% Outputs comparison PNGs under both:
%   pipeline2_figures_legacy_event_density
%   summary_figures\pipeline2_legacy_event_density
%
% Optional controls, set before run(...) if needed:
%   data_root = 'E:\DataPons';
%   output_root = 'E:\DataPons_processed';
%   density_bin_sec = 2;
%   density_smooth_sigma_sec = 2;
%   close_figures = true;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('data_root', 'var') || isempty(data_root)
    data_root = 'E:\DataPons';
end
if ~exist('output_root', 'var') || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end
if ~exist('density_bin_sec', 'var') || isempty(density_bin_sec)
    density_bin_sec = 2;
end
if ~exist('density_smooth_sigma_sec', 'var') || isempty(density_smooth_sigma_sec)
    density_smooth_sigma_sec = density_bin_sec;
end
if ~exist('close_figures', 'var') || isempty(close_figures)
    close_figures = true;
end

records = discover_legacy_pl_event_files(data_root);
if isempty(records)
    fprintf('No legacy mevt_pl/sevt_pl files found under:\n  %s\n', data_root);
    return;
end

summary_dir = fullfile(io_project.get_summary_figures_root(output_root), ...
    'pipeline2_legacy_event_density');
if exist(summary_dir, 'dir') ~= 7
    mkdir(summary_dir);
end

fprintf('Found %d dataset/group legacy event set(s).\n', numel(records));

for i = 1:numel(records)
    rec = records(i);
    cfg = build_legacy_cfg(rec);

    fprintf('\n[%d/%d] %s %s\n', i, numel(records), cfg.dataset_id, rec.group_name);

    params = struct();
    params.bin_sec = double(density_bin_sec);
    params.smooth_sigma_sec = double(density_smooth_sigma_sec);
    params.group_name = rec.group_name;

    E_mevt = [];
    E_sevt = [];

    if ~isempty(rec.mevt_file)
        params.legacy_variable = 'mevt_pl';
        E_mevt = compute_legacy_blp_event_density( ...
            cfg, output_root, rec.mevt_file, 'mevt_pl', params);
        fprintf('  mevt_pl -> %s\n', E_mevt.save_file);
    else
        fprintf('  mevt_pl missing, skipped.\n');
    end

    if ~isempty(rec.sevt_file)
        params.legacy_variable = 'sevt_pl';
        E_sevt = compute_legacy_blp_event_density( ...
            cfg, output_root, rec.sevt_file, 'sevt_pl', params);
        fprintf('  sevt_pl -> %s\n', E_sevt.save_file);
    else
        fprintf('  sevt_pl missing, skipped.\n');
    end

    fig_dir = io_project.get_pipeline_stage_dir(output_root, cfg, 2, ...
        'figures_legacy_event_density');
    if exist(fig_dir, 'dir') ~= 7
        mkdir(fig_dir);
    end

    fig_name = sprintf('%s_%s_legacy_pl_event_density.png', cfg.file_stem, rec.group_name);
    fig_file = fullfile(fig_dir, fig_name);
    summary_file = fullfile(summary_dir, fig_name);

    plot_title = sprintf('%s %s legacy PL event density', cfg.dataset_id, rec.group_name);
    fig = plot_legacy_blp_event_density_pair(E_mevt, E_sevt, plot_title, fig_file);
    copyfile(fig_file, summary_file);

    if close_figures && ~isempty(fig) && isvalid(fig)
        close(fig);
    end

    fprintf('  figure -> %s\n', fig_file);
    fprintf('  summary -> %s\n', summary_file);
end


function records = discover_legacy_pl_event_files(data_root)
files = [ ...
    dir(fullfile(data_root, '**', '*_mevt_pl.mat')); ...
    dir(fullfile(data_root, '**', '*_sevt_pl.mat'))];

records = struct( ...
    'dataset_dir', {}, ...
    'dataset_id', {}, ...
    'file_stem', {}, ...
    'group_name', {}, ...
    'mevt_file', {}, ...
    'sevt_file', {});

for i = 1:numel(files)
    if files(i).isdir
        continue;
    end

    name = files(i).name;
    tok = regexp(name, ...
        '^(?<stem>[A-Za-z0-9]+)_(?<group>spont\d*|vspont)_(?<detector>mevt_pl|sevt_pl)\.mat$', ...
        'names');
    if isempty(tok)
        continue;
    end

    group_dir = files(i).folder;
    [dataset_dir, group_name_from_dir] = fileparts(group_dir);
    if ~strcmpi(group_name_from_dir, tok.group)
        continue;
    end

    [~, dataset_id] = fileparts(dataset_dir);
    output_stem = dataset_id_to_stem(dataset_id);
    key_idx = find_record(records, dataset_dir, output_stem, tok.group);
    if isempty(key_idx)
        key_idx = numel(records) + 1;
        records(key_idx).dataset_dir = dataset_dir;
        records(key_idx).dataset_id = dataset_id;
        records(key_idx).file_stem = output_stem;
        records(key_idx).group_name = tok.group;
        records(key_idx).mevt_file = '';
        records(key_idx).sevt_file = '';
    end

    full_path = fullfile(files(i).folder, files(i).name);
    switch lower(tok.detector)
        case 'mevt_pl'
            records(key_idx).mevt_file = full_path;
        case 'sevt_pl'
            records(key_idx).sevt_file = full_path;
    end
end
end


function idx = find_record(records, dataset_dir, file_stem, group_name)
idx = [];
for k = 1:numel(records)
    if strcmpi(records(k).dataset_dir, dataset_dir) && ...
            strcmpi(records(k).file_stem, file_stem) && ...
            strcmpi(records(k).group_name, group_name)
        idx = k;
        return;
    end
end
end


function cfg = build_legacy_cfg(rec)
cfg = struct();
cfg.dataset_id = rec.dataset_id;
cfg.raw_data_root = [rec.dataset_dir, filesep];
cfg.file_stem = rec.file_stem;
end


function file_stem = dataset_id_to_stem(dataset_id)
file_stem = lower(regexprep(char(string(dataset_id)), '[^A-Za-z0-9]', ''));
end
