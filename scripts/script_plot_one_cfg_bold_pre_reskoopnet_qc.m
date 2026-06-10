% Plot BOLD observable QC figures before ResKoopNet for one cfg_*.m entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   observable_modes = {'global_svd100', 'gsvd100_ds', 'HP_svd100', 'roi_mean'};
%   run('scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m');

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

close all force;

cfg_name = char(string(cfg_name));

if ~exist('observable_modes', 'var') || isempty(observable_modes)
    observable_modes = {'global_svd100', 'gsvd100_ds', 'HP_svd100', 'roi_mean'};
end
observable_modes = cellstr(string(observable_modes(:)).');

if ~exist('save_qc', 'var') || isempty(save_qc)
    save_qc = true;
end
if ~exist('qc_visible', 'var') || isempty(qc_visible)
    qc_visible = 'off';
end

cfg_fun = ['cfg_' cfg_name];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
processed_root = io_project.get_project_processed_root();
input_dir = io_project.get_pipeline_stage_dir(processed_root, cfg.file_stem, 3, 'bold_observables');
dataset_qc_root = io_project.get_pipeline_stage_dir( ...
    processed_root, cfg.file_stem, 3, 'figures_bold_pre_reskoopnet_qc');
central_qc_root = fullfile(processed_root, 'bold_observables_qc', cfg.dataset_id);

params = struct();
params.save = logical(save_qc);
params.visible = qc_visible;
params.close_after_save = true;
params.max_variables = 12;
params.max_heatmap_variables = Inf;
params.max_heatmap_samples = 4000;
params.segment_session = 'longest';
params.segment_start_sec = 0;
params.segment_duration_sec = 300;
params.make_session_gallery = false;
params.save_formats = {'png'};

for i_mode = 1:numel(observable_modes)
    mode_name = observable_modes{i_mode};
    input_file = fullfile(input_dir, ...
        sprintf('%s_bold_observables_%s.mat', cfg.dataset_id, mode_name));

    if exist(input_file, 'file') ~= 2
        warning('Observable file not found, skipping QC: %s', input_file);
        continue;
    end

    params_mode = params;
    params_mode.output_dir = fullfile(central_qc_root, mode_name);
    params_mode.figure_prefix = sprintf('%s_%s_pre_reskoopnet_qc', cfg.dataset_id, mode_name);
    if strcmpi(mode_name, 'slow_band_power')
        params_mode.max_heatmap_variables = 200;
    end

    fprintf('Plotting BOLD pre-ResKoopNet QC: %s | %s\n', cfg.dataset_id, mode_name);
    [~, info] = plot_bold_pre_reskoopnet_qc(input_file, params_mode);
    fprintf('  Observable matrix: [%d samples x %d observables]\n', ...
        info.n_samples, info.n_variables);
    fprintf('  Saved: %s\n', strjoin(cellstr(info.output_files), sprintf('\n         ')));

    dataset_qc_dir = fullfile(dataset_qc_root, mode_name);
    copied_files = local_copy_qc_files(info.output_files, dataset_qc_dir);
    fprintf('  Copied to dataset folder: %s\n', ...
        strjoin(cellstr(copied_files), sprintf('\n                            ')));
end


function copied_files = local_copy_qc_files(source_files, target_dir)
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
    copied_files(end + 1, 1) = string(dst); %#ok<SAGROW>
end
end
