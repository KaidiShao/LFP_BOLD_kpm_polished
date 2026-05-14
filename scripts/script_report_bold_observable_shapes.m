% Report true BOLD observable matrix sizes for one or more cfg_*.m definitions.

clear;
clc;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

processed_root = io_project.get_project_processed_root();

if ~exist('cfg_names', 'var') || isempty(cfg_names)
    cfg_names = {'E10gb1', 'E10gH1', 'E10fV1', 'F12m01'};
end
cfg_names = cellstr(string(cfg_names(:)).');

if ~exist('observable_modes', 'var') || isempty(observable_modes)
    observable_modes = {'eleHP', 'HP', 'roi_mean', 'roi_mean_slow_band_power', ...
        'slow_band_power', 'slow_band_power_svd', 'svd', 'HP_svd100', ...
        'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100'};
end
observable_modes = cellstr(string(observable_modes(:)).');

rows = cell(0, 6);
dataset_report_files = strings(0, 1);

fprintf('BOLD observable shape report\n');
fprintf('Processed root:\n  %s\n\n', processed_root);
fprintf('%-10s %-16s %12s %14s %12s %s\n', ...
    'dataset', 'mode', 'samples', 'observables', 'MB', 'file');

for i_dataset = 1:numel(cfg_names)
    cfg_fun = ['cfg_' cfg_names{i_dataset}];
    if exist(cfg_fun, 'file') ~= 2
        error('Config function not found on path: %s', cfg_fun);
    end

    cfg = feval(cfg_fun);
    dataset_id = cfg.dataset_id;
    dataset_dir = io_project.get_pipeline_stage_dir(processed_root, cfg, 3, 'bold_observables');
    dataset_rows = cell(0, 6);

    for i_mode = 1:numel(observable_modes)
        mode_name = observable_modes{i_mode};
        input_file = fullfile(dataset_dir, ...
            sprintf('%s_bold_observables_%s.mat', dataset_id, mode_name));

        if exist(input_file, 'file') ~= 2
            warning('Observable file not found, skipping: %s', input_file);
            continue;
        end

        [n_samples, n_observables] = read_obs_size(input_file);
        file_info = dir(input_file);
        file_mb = file_info.bytes / 1024 / 1024;

        fprintf('%-10s %-16s %12d %14d %12.1f %s\n', ...
            dataset_id, mode_name, n_samples, n_observables, file_mb, input_file);

        dataset_rows(end+1, :) = {dataset_id, mode_name, n_samples, ...
            n_observables, file_mb, input_file}; %#ok<SAGROW>
        rows(end+1, :) = {dataset_id, mode_name, n_samples, ...
            n_observables, file_mb, input_file}; %#ok<SAGROW>
    end

    if ~isempty(dataset_rows)
        T_dataset = cell2table(dataset_rows, 'VariableNames', ...
            {'dataset_id', 'observable_mode', 'n_samples', ...
            'n_observables', 'file_mb', 'file'});
        dataset_output_file = fullfile(dataset_dir, ...
            sprintf('%s_bold_observable_shape_report.csv', cfg.file_stem));
        writetable(T_dataset, dataset_output_file);
        dataset_report_files(end+1, 1) = string(dataset_output_file); %#ok<SAGROW>
    end
end

if exist('output_file', 'var') && ~isempty(output_file)
    T = cell2table(rows, 'VariableNames', ...
        {'dataset_id', 'observable_mode', 'n_samples', ...
        'n_observables', 'file_mb', 'file'});
    writetable(T, output_file);
    fprintf('\nSaved combined report:\n  %s\n', output_file);
end

if ~isempty(dataset_report_files)
    fprintf('\nSaved dataset reports:\n');
    fprintf('  %s\n', strjoin(cellstr(dataset_report_files), sprintf('\n  ')));
end

function [n_samples, n_observables] = read_obs_size(input_file)
try
    M = matfile(input_file);
    sz = size(M, 'obs');
catch
    S = load(input_file, 'obs');
    sz = size(S.obs);
end

if numel(sz) < 2
    sz(2) = 1;
end
n_samples = sz(1);
n_observables = sz(2);
end
