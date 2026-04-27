% Batch-generate BOLD observables and session-aware snapshots.
%
% This script reads session files from each dataset's roits folder. It does
% not modify cfg files.

this_script_dir = fileparts(mfilename('fullpath'));
repo_root = fileparts(this_script_dir);
addpath(genpath(repo_root));

output_root = fullfile(get_project_processed_root(), 'bold_observables');
if exist(output_root, 'dir') ~= 7
    mkdir(output_root);
end

% -------------------------------------------------------------------------
% Manual BOLD job list.
% Exact semantic ROI names come from cfg.bold.role_map and are validated
% against the current roiTs{1,r}.name values during loading.
% -------------------------------------------------------------------------
bold_jobs = struct([]);

bold_jobs(1).cfg_name = 'E10gb1';

bold_jobs(2).cfg_name = 'E10gH1';

bold_jobs(3).cfg_name = 'E10fV1';

bold_jobs(4).cfg_name = 'F12m01';

% Observable branches to generate.
% eleHP and HP are separate voxel-level observable branches.
% roi_mean uses one mean trace per ROI across all available ROIs.
% slow_band_power uses eleHP+HP voxels.
% svd uses all available ROI voxels, not ROI-mean traces.
% HP_svd100 uses only HP voxels. global_svd100 uses all available ROI voxels.
observable_modes = {'eleHP', 'HP', 'roi_mean', 'slow_band_power', ...
    'svd', 'HP_svd100', 'global_svd100'};
n_svd_components = 50;
n_svd100_components = 100;
snapshot_lag = 1;
save_precision = 'single';
force_recompute = false;

% Minimal preprocessing before observable construction.
pre = struct();
pre.demean = false;
pre.zscore = false;
pre.detrend = false;
pre.notch_hz = [];
pre.bandpass_hz = [];
pre.save_filtered = false;

for i_job = 1:numel(bold_jobs)
    job = bold_jobs(i_job);
    validate_bold_job(job);

    cfg_fun = ['cfg_' job.cfg_name];
    if exist(cfg_fun, 'file') ~= 2
        error('Config function not found on path: %s', cfg_fun);
    end

    cfg = feval(cfg_fun);
    cfg = local_apply_bold_defaults(cfg);
    role_map = local_require_bold_role_map(cfg);

    fprintf('\n============================================================\n');
    fprintf('Building BOLD observables for %s\n', cfg.dataset_id);
    fprintf('BOLD roits folder:\n  %s\n', fullfile(cfg.raw_data_root, cfg.bold.data_subfolder));
    fprintf('Configured BOLD ROI role map:\n');
    fprintf('  eleHP -> %s\n', role_map.elehp);
    fprintf('  HP    -> %s\n', role_map.hp);

    dataset_out_dir = fullfile(output_root, cfg.dataset_id);
    if exist(dataset_out_dir, 'dir') ~= 7
        mkdir(dataset_out_dir);
    end

    for i_mode = 1:numel(observable_modes)
        mode_name = observable_modes{i_mode};
        out_file = fullfile(dataset_out_dir, ...
            sprintf('%s_bold_observables_%s.mat', cfg.dataset_id, mode_name));

        if exist(out_file, 'file') == 2 && ~force_recompute
            fprintf('Skipping existing %s\n  %s\n', mode_name, out_file);
            continue;
        end

        cfg_mode = cfg;
        cfg_mode.bold.variable_mode = 'voxels';
        cfg_mode.bold.selected_region_names = select_region_names_for_observable_mode(role_map, mode_name);

        fprintf('Observable mode %s uses variable_mode=%s, selected_region_names=%s\n', ...
            mode_name, cfg_mode.bold.variable_mode, strjoin(cfg_mode.bold.selected_region_names, ', '));

        D = load_bold_dataset(cfg_mode);
        P = preprocess_bold_sessions(D, pre);

        params = struct();
        params.source = 'data';
        params.observable_branch = mode_name;
        if ismember(lower(mode_name), {'elehp', 'hp'})
            params.mode = 'identity';
        elseif ismember(lower(mode_name), {'hp_svd100', 'global_svd100'})
            params.mode = 'svd';
        else
            params.mode = mode_name;
        end
        params.precision = save_precision;
        if ismember(lower(mode_name), {'hp_svd100', 'global_svd100'})
            params.n_components = n_svd100_components;
        else
            params.n_components = n_svd_components;
        end

        if strcmpi(mode_name, 'slow_band_power')
            params.include_raw = true;
            params.band_power_transform = 'power';
        end

        O = build_bold_observables(P, params);
        snap = make_session_aware_snapshot_pairs(O, snapshot_lag);

        obs = O.data;
        obs_info = O.observable_info;
        dx = O.dx;
        dt = O.dx;
        fs = [];
        if isfield(O, 'fs')
            fs = O.fs;
        end
        sampling_period = dx;
        sample_period = dx;
        sampling_frequency = fs;
        session_ids = O.session_ids;
        session_lengths = O.session_lengths;
        session_dx = O.session_dx;
        session_start_idx = O.session_start_idx;
        session_end_idx = O.session_end_idx;
        border_idx = O.border_idx;
        snapshot_lag_saved = snapshot_lag;
        snapshot_valid_idx_x = snap.valid_idx_x;
        snapshot_valid_idx_y = snap.valid_idx_y;
        snapshot_session_idx = snap.session_idx;
        snapshot_session_id = snap.session_id;
        cfg = cfg_mode;

        save(out_file, ...
            'obs', 'obs_info', ...
            'dx', 'dt', 'fs', 'sampling_period', 'sample_period', 'sampling_frequency', ...
            'session_ids', 'session_lengths', 'session_dx', ...
            'session_start_idx', 'session_end_idx', 'border_idx', ...
            'snapshot_lag_saved', 'snapshot_valid_idx_x', 'snapshot_valid_idx_y', ...
            'snapshot_session_idx', 'snapshot_session_id', ...
            'O', 'snap', 'params', 'pre', 'cfg', '-v7.3');
        fprintf('Saved %s observables: [%d x %d]\n  %s\n', ...
            mode_name, size(O.data, 1), size(O.data, 2), out_file);
        fprintf('Snapshot pairs: X [%d x %d], Y [%d x %d]\n', ...
            size(snap.X, 1), size(snap.X, 2), size(snap.Y, 1), size(snap.Y, 2));
    end
end

fprintf('\nFinished BOLD observable batch: %d dataset(s), %d mode(s).\n', ...
    numel(bold_jobs), numel(observable_modes));

function validate_bold_job(job)
if ~isfield(job, 'cfg_name') || isempty(job.cfg_name)
    error('bold_jobs cfg_name is empty.');
end
end

function cfg = local_apply_bold_defaults(cfg)
if ~isfield(cfg, 'bold') || ~isstruct(cfg.bold)
    cfg.bold = struct();
end
if ~isfield(cfg.bold, 'data_subfolder') || isempty(cfg.bold.data_subfolder)
    cfg.bold.data_subfolder = 'roits';
end
if ~isfield(cfg.bold, 'input_var') || isempty(cfg.bold.input_var)
    cfg.bold.input_var = 'roiTs';
end
end

function role_map = local_require_bold_role_map(cfg)
if ~isfield(cfg, 'bold') || ~isstruct(cfg.bold) || ...
        ~isfield(cfg.bold, 'role_map') || ~isstruct(cfg.bold.role_map)
    error(['cfg.bold.role_map is required for BOLD observable batching.\n' ...
        'Expected fields:\n' ...
        '  cfg.bold.role_map.elehp\n' ...
        '  cfg.bold.role_map.hp']);
end

required_fields = {'elehp', 'hp'};
role_map = struct();
for i = 1:numel(required_fields)
    field_name = required_fields{i};
    if ~isfield(cfg.bold.role_map, field_name) || isempty(cfg.bold.role_map.(field_name))
        error('cfg.bold.role_map.%s is required for dataset %s.', ...
            field_name, cfg.dataset_id);
    end
    role_map.(field_name) = char(string(cfg.bold.role_map.(field_name)));
end
end

function selected_region_names = select_region_names_for_observable_mode(role_map, mode_name)
switch lower(mode_name)
    case 'elehp'
        selected_region_names = {role_map.elehp};
    case 'hp'
        selected_region_names = {role_map.hp};
    case 'hp_svd100'
        selected_region_names = {role_map.hp};
    case 'roi_mean'
        selected_region_names = {};
    case 'slow_band_power'
        selected_region_names = {role_map.elehp, role_map.hp};
    case {'svd', 'global_svd100'}
        selected_region_names = {};
    otherwise
        error('Unsupported observable mode: %s', mode_name);
end
selected_region_names = cellstr(string(selected_region_names(:)).');
end
