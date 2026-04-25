% Plot BOLD observable QC figures before ResKoopNet.
%
% Run from MATLAB:
%   cd('D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished')
%   run('scripts\script_plot_bold_pre_reskoopnet_qc.m')

clear;
clc;

project_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(project_root));

input_root = fullfile(get_project_processed_root(), 'bold_observables');
output_root = fullfile(get_project_processed_root(), 'bold_observables_qc');

dataset_ids = {'E10.gb1', 'E10.gH1', 'E10.fV1', 'F12.m01'};
processed_dataset_dirs = {'e10gb1', 'e10gh1', 'e10fV1', 'f12m01'};
observable_modes = {'eleHP', 'HP', 'roi_mean', 'slow_band_power', ...
    'svd', 'HP_svd100', 'global_svd100'};

params = struct();
params.save = true;
params.visible = 'off';
params.close_after_save = true;
params.max_variables = 12;
params.max_heatmap_variables = Inf;
params.max_heatmap_samples = 4000;
params.segment_session = 'longest';
params.segment_start_sec = 0;
params.segment_duration_sec = 300;
params.make_session_gallery = false;
params.save_formats = {'png'};

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

        params.output_dir = fullfile(output_root, dataset_id, mode_name);
        params.figure_prefix = sprintf('%s_%s_pre_reskoopnet_qc', dataset_id, mode_name);
        params_mode = params;
        if strcmpi(mode_name, 'slow_band_power')
            params_mode.max_heatmap_variables = 200;
        end

        fprintf('Plotting BOLD pre-ResKoopNet QC: %s | %s\n', dataset_id, mode_name);
        [~, info] = plot_bold_pre_reskoopnet_qc(input_file, params_mode);
        fprintf('  Observable matrix: [%d samples x %d observables]\n', ...
            info.n_samples, info.n_variables);
        fprintf('  Saved: %s\n', strjoin(cellstr(info.output_files), sprintf('\n         ')));

        dataset_qc_dir = fullfile(get_project_processed_root(), ...
            processed_dataset_dirs{i_dataset}, 'bold_observables_qc', mode_name);
        copied_files = copy_qc_files_to_dataset_dir(info.output_files, dataset_qc_dir);
        fprintf('  Copied to dataset folder: %s\n', strjoin(cellstr(copied_files), sprintf('\n                            ')));

        flat_files = copy_qc_files_to_dataset_dir(info.output_files, project_root);
        fprintf('  Copied to project root: %s\n', strjoin(cellstr(flat_files), sprintf('\n                          ')));
    end
end

fprintf('\nBOLD QC figures saved under:\n  %s\n', output_root);

function copied_files = copy_qc_files_to_dataset_dir(source_files, target_dir)
if exist(target_dir, 'dir') ~= 7
    mkdir(target_dir);
end

source_files = string(source_files(:));
copied_files = strings(0, 1);
for i_file = 1:numel(source_files)
    src = char(source_files(i_file));
    if exist(src, 'file') ~= 2
        warning('QC figure file not found, cannot copy: %s', src);
        continue;
    end

    [~, name, ext] = fileparts(src);
    dst = fullfile(target_dir, [name, ext]);
    copyfile(src, dst);
    copied_files(end+1, 1) = string(dst); %#ok<SAGROW>
end
end
