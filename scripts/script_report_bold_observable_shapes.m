% Report true BOLD observable matrix sizes.
%
% Run from MATLAB:
%   cd('D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished')
%   run('scripts\script_report_bold_observable_shapes.m')

clear;
clc;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

input_root = fullfile(get_project_processed_root(), 'bold_observables');
output_file = fullfile(input_root, 'bold_observable_shape_report.csv');

dataset_ids = {'E10.gb1', 'E10.gH1', 'E10.fV1', 'F12.m01'};
observable_modes = {'eleHP', 'HP', 'roi_mean', 'slow_band_power', ...
    'svd', 'HP_svd100', 'global_svd100'};

rows = cell(0, 6);

fprintf('BOLD observable shape report\n');
fprintf('Input root:\n  %s\n\n', input_root);
fprintf('%-10s %-16s %12s %14s %12s %s\n', ...
    'dataset', 'mode', 'samples', 'observables', 'MB', 'file');

for i_dataset = 1:numel(dataset_ids)
    dataset_id = dataset_ids{i_dataset};
    dataset_dir = fullfile(input_root, dataset_id);

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

        rows(end+1, :) = {dataset_id, mode_name, n_samples, ...
            n_observables, file_mb, input_file}; %#ok<SAGROW>
    end
end

T = cell2table(rows, 'VariableNames', ...
    {'dataset_id', 'observable_mode', 'n_samples', ...
    'n_observables', 'file_mb', 'file'});
writetable(T, output_file);

fprintf('\nSaved report:\n  %s\n', output_file);

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
