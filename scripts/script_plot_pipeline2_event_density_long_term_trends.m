% Batch plot long-term pipeline2 event-density summaries from saved outputs.
%
% Usage from MATLAB:
%   run('scripts/script_plot_pipeline2_event_density_long_term_trends.m');
%
% Optional controls, set before run(...) if needed:
%   event_density_files = {'E:\DataPons_processed\e10gb1\pipeline2_event_density\e10gb1_event_density_2s.mat'};
%   plot_params.visible = 'on';
%   plot_params.time_unit = 'hours';
%   plot_params.trend_window_sec = 300;
%   plot_params.show_raw_mean = false;

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;
clc;

output_root = io_project.get_project_processed_root();
summary_dir = io_project.get_pipeline_summary_dir(output_root, 2, 'event_density');
if exist(summary_dir, 'dir') ~= 7
    mkdir(summary_dir);
end

if ~exist('event_density_files', 'var') || isempty(event_density_files)
    event_density_files = local_discover_event_density_files(output_root);
else
    event_density_files = cellstr(string(event_density_files(:)));
end

if isempty(event_density_files)
    error('No saved pipeline2 event-density files were found under: %s', output_root);
end

if ~exist('plot_params', 'var') || isempty(plot_params)
    plot_params = struct();
end

if ~isfield(plot_params, 'save') || isempty(plot_params.save)
    plot_params.save = true;
end
if ~isfield(plot_params, 'output_dir') || isempty(plot_params.output_dir)
    plot_params.output_dir = summary_dir;
end
if ~isfield(plot_params, 'visible') || isempty(plot_params.visible)
    plot_params.visible = 'off';
end
if ~isfield(plot_params, 'close_after_save') || isempty(plot_params.close_after_save)
    plot_params.close_after_save = true;
end
if ~isfield(plot_params, 'time_unit') || isempty(plot_params.time_unit)
    plot_params.time_unit = 'hours';
end
if ~isfield(plot_params, 'trend_window_sec') || isempty(plot_params.trend_window_sec)
    plot_params.trend_window_sec = 300;
end
if ~isfield(plot_params, 'show_session_boundaries') || isempty(plot_params.show_session_boundaries)
    plot_params.show_session_boundaries = true;
end
if ~isfield(plot_params, 'show_raw_mean') || isempty(plot_params.show_raw_mean)
    plot_params.show_raw_mean = false;
end

n_files = numel(event_density_files);
batch_rows = repmat(local_empty_manifest_row(), n_files, 1);

for i = 1:n_files
    input_file = char(event_density_files{i});
    fprintf('\n============================================================\n');
    fprintf('Plotting pipeline2 event-density long-term trend (%d/%d)\n', i, n_files);
    fprintf('  %s\n', input_file);

    batch_rows(i).input_file = string(input_file);

    try
        [~, info] = plot_blp_event_density_long_term_trend(input_file, plot_params);
        batch_rows(i) = local_populate_manifest_row(batch_rows(i), info);
        batch_rows(i).status = "ok";
        fprintf('Saved figure(s):\n');
        disp(info.output_files);
    catch ME
        batch_rows(i).status = "error";
        batch_rows(i).error_message = string(ME.message);
        fprintf(2, 'Failed: %s\n', ME.message);
    end
end

plot_manifest = struct2table(batch_rows);
manifest_file = fullfile(plot_params.output_dir, 'pipeline2_event_density_long_term_trend_manifest.csv');
writetable(plot_manifest, manifest_file);

fprintf('\nSaved manifest to:\n  %s\n', manifest_file);


function files = local_discover_event_density_files(output_root)
stage_name = io_project.get_pipeline_stage_name(2, 'event_density');
listing = dir(fullfile(output_root, '*', stage_name, '*_event_density_*.mat'));
files = arrayfun(@(x) string(fullfile(x.folder, x.name)), listing(:), 'UniformOutput', true);
files = sort(cellstr(files));
end


function row = local_empty_manifest_row()
row = struct( ...
    'dataset_id', "", ...
    'file_stem', "", ...
    'input_file', "", ...
    'output_png', "", ...
    'output_pdf', "", ...
    'output_fig', "", ...
    'band_labels', "", ...
    'time_unit', "", ...
    'total_duration_hr', NaN, ...
    'n_sessions', NaN, ...
    'status', "", ...
    'error_message', "");
end


function row = local_populate_manifest_row(row, info)
if isfield(info, 'dataset_id') && ~isempty(info.dataset_id)
    row.dataset_id = string(info.dataset_id);
end
if isfield(info, 'file_stem') && ~isempty(info.file_stem)
    row.file_stem = string(info.file_stem);
end
if isfield(info, 'band_labels') && ~isempty(info.band_labels)
    row.band_labels = strjoin(cellstr(info.band_labels(:).'), ',');
end
if isfield(info, 'time_unit') && ~isempty(info.time_unit)
    row.time_unit = string(info.time_unit);
end
if isfield(info, 'total_duration_sec') && ~isempty(info.total_duration_sec)
    row.total_duration_hr = double(info.total_duration_sec) / 3600;
end
if isfield(info, 'n_sessions') && ~isempty(info.n_sessions)
    row.n_sessions = double(info.n_sessions);
end
if isfield(info, 'output_files') && ~isempty(info.output_files)
    output_files = cellstr(info.output_files(:));
    for i = 1:numel(output_files)
        [~, ~, ext] = fileparts(output_files{i});
        switch lower(ext)
            case '.png'
                row.output_png = string(output_files{i});
            case '.pdf'
                row.output_pdf = string(output_files{i});
            case '.fig'
                row.output_fig = string(output_files{i});
        end
    end
end
end
