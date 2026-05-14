% Distribute saved BOLD pre-ResKoopNet QC figures after plotting.
%
% Usage from MATLAB:
%   cfg_name = 'E10gH1';
%   run('scripts/script_distribute_one_cfg_bold_pre_reskoopnet_qc_figures.m');
%
% Optional controls, set before run(...) if needed:
%   observable_modes = {'eleHP', 'HP', 'roi_mean'};
%   qc_copy_to_summary_dir = true;
%   qc_copy_to_project_root = false;

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
project_root = fileparts(fileparts(mfilename('fullpath')));

if ~exist('observable_modes', 'var') || isempty(observable_modes)
    observable_modes = {'eleHP', 'HP', 'roi_mean', 'roi_mean_slow_band_power', ...
        'slow_band_power', 'slow_band_power_svd', 'svd', 'HP_svd100', ...
        'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100'};
end
observable_modes = cellstr(string(observable_modes(:)).');

copy_to_summary_dir = true;
copy_to_project_root = false;
if exist('qc_copy_to_summary_dir', 'var') && ~isempty(qc_copy_to_summary_dir)
    copy_to_summary_dir = logical(qc_copy_to_summary_dir);
end
if exist('qc_copy_to_project_root', 'var') && ~isempty(qc_copy_to_project_root)
    copy_to_project_root = logical(qc_copy_to_project_root);
end
if ~copy_to_summary_dir && ~copy_to_project_root
    error(['Nothing to do: enable qc_copy_to_summary_dir and/or ' ...
        'qc_copy_to_project_root before running this script.']);
end

qc_stage_dir = io_project.get_pipeline_stage_dir(processed_root, cfg, 3, 'figures_bold_pre_reskoopnet_qc');
summary_qc_dir = io_project.get_pipeline_summary_dir(processed_root, 3, 'figures_bold_pre_reskoopnet_qc_summary');
qc_distribution_outputs = struct([]);

for i_mode = 1:numel(observable_modes)
    mode_name = observable_modes{i_mode};
    mode_dir = fullfile(qc_stage_dir, mode_name);
    if exist(mode_dir, 'dir') ~= 7
        warning('QC figure folder not found, skipping: %s', mode_dir);
        continue;
    end

    source_files = local_collect_qc_files(mode_dir);
    if isempty(source_files)
        warning('No QC figure files found in: %s', mode_dir);
        continue;
    end

    fprintf('Distributing BOLD QC figures: %s | %s\n', cfg.dataset_id, mode_name);

    summary_copies = strings(0, 1);
    if copy_to_summary_dir
        summary_copies = local_copy_qc_files(source_files, summary_qc_dir);
        fprintf('  Copied to summary folder: %s\n', ...
            strjoin(cellstr(summary_copies), sprintf('\n                           ')));
    end

    root_copies = strings(0, 1);
    if copy_to_project_root
        root_copies = local_copy_qc_files(source_files, project_root);
        fprintf('  Copied to project root: %s\n', ...
            strjoin(cellstr(root_copies), sprintf('\n                          ')));
    end

    qc_distribution_outputs(i_mode).mode_name = mode_name; %#ok<SAGROW>
    qc_distribution_outputs(i_mode).source_files = source_files; %#ok<SAGROW>
    qc_distribution_outputs(i_mode).summary_copies = summary_copies; %#ok<SAGROW>
    qc_distribution_outputs(i_mode).project_root_copies = root_copies; %#ok<SAGROW>
end

fprintf('\nBOLD QC figure distribution finished for %s.\n', cfg.dataset_id);


function source_files = local_collect_qc_files(mode_dir)
listing = dir(fullfile(mode_dir, '*'));
listing = listing(~[listing.isdir]);
source_files = strings(0, 1);
for i_file = 1:numel(listing)
    source_files(end+1, 1) = string(fullfile(listing(i_file).folder, listing(i_file).name)); %#ok<SAGROW>
end
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
    copied_files(end+1, 1) = string(dst); %#ok<SAGROW>
end
end
