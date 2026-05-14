% Plot three-row event-density comparisons for datasets with pipeline2 output.
%
% Rows:
%   1. current pipeline2 event density
%   2. legacy mevt_pl event density
%   3. legacy sevt_pl event density
%
% Reads:
%   E:\DataPons_processed\<dataset>\pipeline2_event_density
%   E:\DataPons_processed\<dataset>\pipeline2_legacy_event_density
%
% Writes:
%   E:\DataPons_processed\<dataset>\pipeline2_figures_legacy_event_density
%   E:\DataPons_processed\summary_figures\pipeline2_pipeline2_plus_legacy_event_density
%
% Optional controls:
%   output_root = 'E:\DataPons_processed';
%   close_figures = true;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

if ~exist('output_root', 'var') || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end
if ~exist('close_figures', 'var') || isempty(close_figures)
    close_figures = true;
end

summary_dir = fullfile(io_project.get_summary_figures_root(output_root), ...
    'pipeline2_pipeline2_plus_legacy_event_density');
if exist(summary_dir, 'dir') ~= 7
    mkdir(summary_dir);
end

pipeline_files = dir(fullfile(output_root, '*', 'pipeline2_event_density', '*_event_density_2s.mat'));
for i = 1:numel(pipeline_files)
    pipeline_file = fullfile(pipeline_files(i).folder, pipeline_files(i).name);
    dataset_dir = fileparts(fileparts(pipeline_file));
    [~, dataset_dir_name] = fileparts(dataset_dir);
    file_stem = lower(dataset_dir_name);

    S = load(pipeline_file, 'E');
    E_pipeline2 = S.E;

    legacy_dir = fullfile(dataset_dir, 'pipeline2_legacy_event_density');
    if exist(legacy_dir, 'dir') ~= 7
        continue;
    end

    groups = discover_legacy_groups(legacy_dir, file_stem);
    for g = 1:numel(groups)
        group_name = groups{g};
        E_mevt = load_optional_legacy(legacy_dir, file_stem, group_name, 'mevt_pl');
        E_sevt = load_optional_legacy(legacy_dir, file_stem, group_name, 'sevt_pl');
        if isempty(E_mevt) && isempty(E_sevt)
            continue;
        end

        fig_dir = fullfile(dataset_dir, 'pipeline2_figures_legacy_event_density');
        if exist(fig_dir, 'dir') ~= 7
            mkdir(fig_dir);
        end

        fig_name = sprintf('%s_%s_pipeline2_plus_legacy_pl_event_density.png', file_stem, group_name);
        fig_file = fullfile(fig_dir, fig_name);
        summary_file = fullfile(summary_dir, fig_name);
        plot_title = sprintf('%s %s pipeline2 + legacy PL event density', file_stem, group_name);

        fig = plot_pipeline2_legacy_event_density_triplet( ...
            E_pipeline2, E_mevt, E_sevt, plot_title, fig_file);
        copyfile(fig_file, summary_file);

        if close_figures && ~isempty(fig) && isvalid(fig)
            close(fig);
        end

        fprintf('Saved comparison figure:\n  %s\n', fig_file);
    end
end


function groups = discover_legacy_groups(legacy_dir, file_stem)
files = dir(fullfile(legacy_dir, sprintf('%s_*_legacy_event_density_2s.mat', file_stem)));
if isempty(files)
    files = dir(fullfile(legacy_dir, sprintf('%s_*_legacy_event_density_2s.mat', lower(file_stem))));
end
groups = cell(0, 1);
expr = sprintf('^%s_(?<group>.+?)_(?:mevt_pl|sevt_pl)_legacy_event_density_2s\\.mat$', regexptranslate('escape', lower(file_stem)));

for i = 1:numel(files)
    tok = regexp(files(i).name, expr, 'names');
    if isempty(tok)
        continue;
    end
    groups{end+1, 1} = tok.group; %#ok<AGROW>
end

groups = unique(groups, 'stable');
end


function E = load_optional_legacy(legacy_dir, file_stem, group_name, detector_name)
file_name = sprintf('%s_%s_%s_legacy_event_density_2s.mat', lower(file_stem), group_name, detector_name);
file_path = fullfile(legacy_dir, file_name);
if exist(file_path, 'file') ~= 2
    E = [];
    return;
end

S = load(file_path, 'E');
E = S.E;
end
