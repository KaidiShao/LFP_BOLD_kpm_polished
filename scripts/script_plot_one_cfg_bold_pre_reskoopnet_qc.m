% Canonical interactive BOLD pre-ResKoopNet QC entry.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m');
%
% Optional controls, set before run(...) if needed:
%   observable_modes = {'eleHP', 'HP', 'roi_mean'};
%   save_qc = true;
%   qc_visible = 'off';
%
% For copying saved QC figures to the summary folder or project root, use:
%   scripts/script_distribute_one_cfg_bold_pre_reskoopnet_qc_figures.m

if ~exist('cfg_name', 'var') || isempty(cfg_name)
    error(['Set cfg_name before running this script, for example: ' ...
        'cfg_name = ''E10gH1'';']);
end

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

clear info
clc;

cfg_fun = ['cfg_' char(string(cfg_name))];
if exist(cfg_fun, 'file') ~= 2
    error('Config function not found on path: %s', cfg_fun);
end

cfg = feval(cfg_fun);
processed_root = io_project.get_project_processed_root();

if ~exist('observable_modes', 'var') || isempty(observable_modes)
    observable_modes = {'eleHP', 'HP', 'roi_mean', 'roi_mean_slow_band_power', ...
        'slow_band_power', 'slow_band_power_svd', 'svd', 'HP_svd100', ...
        'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100'};
end
observable_modes = cellstr(string(observable_modes(:)).');

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

if exist('save_qc', 'var') && ~isempty(save_qc)
    params.save = logical(save_qc);
end
if exist('qc_visible', 'var') && ~isempty(qc_visible)
    params.visible = char(string(qc_visible));
end

qc_outputs = struct([]);
dataset_dir = io_project.get_pipeline_stage_dir(processed_root, cfg, 3, 'bold_observables');
qc_stage_dir = io_project.get_pipeline_stage_dir(processed_root, cfg, 3, 'figures_bold_pre_reskoopnet_qc');

for i_mode = 1:numel(observable_modes)
    mode_name = observable_modes{i_mode};
    input_file = fullfile(dataset_dir, ...
        sprintf('%s_bold_observables_%s.mat', cfg.dataset_id, mode_name));

    if exist(input_file, 'file') ~= 2
        warning('Observable file not found, skipping: %s', input_file);
        continue;
    end

    params_mode = params;
    params_mode.output_dir = fullfile(qc_stage_dir, mode_name);
    params_mode.figure_prefix = sprintf('%s_%s_pre_reskoopnet_qc', cfg.dataset_id, mode_name);
    if strcmpi(mode_name, 'slow_band_power')
        params_mode.max_heatmap_variables = 200;
    end

    fprintf('Plotting BOLD pre-ResKoopNet QC: %s | %s\n', cfg.dataset_id, mode_name);
    [~, info] = plot_bold_pre_reskoopnet_qc(input_file, params_mode);
    fprintf('  Observable matrix: [%d samples x %d observables]\n', ...
        info.n_samples, info.n_variables);
    fprintf('  Saved: %s\n', strjoin(cellstr(info.output_files), sprintf('\n         ')));

    qc_outputs(i_mode).mode_name = mode_name; %#ok<SAGROW>
    qc_outputs(i_mode).input_file = input_file; %#ok<SAGROW>
    qc_outputs(i_mode).output_files = string(info.output_files); %#ok<SAGROW>
    qc_outputs(i_mode).n_samples = info.n_samples; %#ok<SAGROW>
    qc_outputs(i_mode).n_variables = info.n_variables; %#ok<SAGROW>
end

fprintf('\nBOLD QC finished for %s.\n', cfg.dataset_id);
