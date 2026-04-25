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
% -------------------------------------------------------------------------
bold_jobs = struct([]);

bold_jobs(1).cfg_name = 'E10gb1';
bold_jobs(1).input_dir = 'E:\DataPons\E10.gb1\roits';
bold_jobs(1).file_stem = 'e10gb1';
bold_jobs(1).input_var = 'roiTs';
bold_jobs(1).dx = 2;
bold_jobs(1).eleHP_region_names = {'eleHP', 'ele_HP'};
bold_jobs(1).HP_region_names = {'HP'};
bold_jobs(1).slow_band_power_region_names = {'eleHP', 'ele_HP', 'HP'};

bold_jobs(2).cfg_name = 'E10gH1';
bold_jobs(2).input_dir = 'E:\DataPons\E10.gH1\roits';
bold_jobs(2).file_stem = 'e10gh1';
bold_jobs(2).input_var = 'roiTs';
bold_jobs(2).dx = 2;
bold_jobs(2).eleHP_region_names = {'eleHP', 'ele_HP'};
bold_jobs(2).HP_region_names = {'HP'};
bold_jobs(2).slow_band_power_region_names = {'eleHP', 'ele_HP', 'HP'};

bold_jobs(3).cfg_name = 'E10fV1';
bold_jobs(3).input_dir = 'E:\DataPons\E10.fV1\roits';
bold_jobs(3).file_stem = 'e10fv1';
bold_jobs(3).input_var = 'roiTs';
bold_jobs(3).dx = 2;
bold_jobs(3).eleHP_region_names = {'eleHP', 'ele_HP'};
bold_jobs(3).HP_region_names = {'HP'};
bold_jobs(3).slow_band_power_region_names = {'eleHP', 'ele_HP', 'HP'};

bold_jobs(4).cfg_name = 'F12m01';
bold_jobs(4).input_dir = 'E:\DataPons\F12.m01\roits';
bold_jobs(4).file_stem = 'f12m01';
bold_jobs(4).input_var = 'roiTs';
bold_jobs(4).dx = 2;
bold_jobs(4).eleHP_region_names = {'eleHP', 'ele_HP'};
bold_jobs(4).HP_region_names = {'HP'};
bold_jobs(4).slow_band_power_region_names = {'eleHP', 'ele_HP', 'HP'};

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
    cfg.bold = struct();
    cfg.bold.input_dir = job.input_dir;
    cfg.bold.file_stem = job.file_stem;
    cfg.bold.input_var = job.input_var;
    cfg.bold.dx = job.dx;

    fprintf('\n============================================================\n');
    fprintf('Building BOLD observables for %s\n', cfg.dataset_id);
    fprintf('BOLD roits folder:\n  %s\n', cfg.bold.input_dir);

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
        cfg_mode.bold.selected_region_names = select_region_names_for_observable_mode(job, mode_name);

        fprintf('Observable mode %s uses variable_mode=%s, selected_region_names=%s\n', ...
            mode_name, cfg_mode.bold.variable_mode, strjoin(cfg_mode.bold.selected_region_names, ', '));

        D = load_bold_roits_folder_dataset(cfg_mode);
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
if ~isfield(job, 'input_dir') || isempty(job.input_dir)
    error('bold_jobs input_dir is empty for %s.', job.cfg_name);
end
if exist(job.input_dir, 'dir') ~= 7
    error('BOLD input_dir does not exist for %s:\n  %s', job.cfg_name, job.input_dir);
end
if ~isfield(job, 'file_stem') || isempty(job.file_stem)
    error('bold_jobs file_stem is empty for %s.', job.cfg_name);
end
required_region_fields = {'eleHP_region_names', 'HP_region_names', 'slow_band_power_region_names'};
for i_field = 1:numel(required_region_fields)
    field_name = required_region_fields{i_field};
    if ~isfield(job, field_name)
        error('bold_jobs.%s is missing for %s.', field_name, job.cfg_name);
    end
end
end

function selected_region_names = select_region_names_for_observable_mode(job, mode_name)
switch lower(mode_name)
    case 'elehp'
        selected_region_names = job.eleHP_region_names;
    case 'hp'
        selected_region_names = job.HP_region_names;
    case 'hp_svd100'
        selected_region_names = job.HP_region_names;
    case 'roi_mean'
        selected_region_names = {};
    case 'slow_band_power'
        selected_region_names = job.slow_band_power_region_names;
    case {'svd', 'global_svd100'}
        selected_region_names = {};
    otherwise
        error('Unsupported observable mode: %s', mode_name);
end
selected_region_names = cellstr(string(selected_region_names(:)).');
end

function D = load_bold_roits_folder_dataset(cfg)
session_ids = collect_included_session_ids(cfg);
if isempty(session_ids)
    error('No included sessions were found in cfg.sessions.');
end

bold = cfg.bold;
data_cells = cell(numel(session_ids), 1);
session_lengths = zeros(numel(session_ids), 1);
session_dx = zeros(numel(session_ids), 1);
source_files = cell(numel(session_ids), 1);

variable_info_ref = table();
region_labels_ref = {};
n_var_region_ref = [];

for i_sess = 1:numel(session_ids)
    session_id = session_ids(i_sess);
    file_path = fullfile(bold.input_dir, sprintf('%s_%04d_roits.mat', ...
        bold.file_stem, session_id));
    if exist(file_path, 'file') ~= 2
        error('Missing BOLD roits file for session %g:\n  %s', session_id, file_path);
    end

    fprintf('Loading BOLD session %g:\n  %s\n', session_id, file_path);
    S = load(file_path, bold.input_var);
    if ~isfield(S, bold.input_var)
        error('Variable "%s" was not found in %s.', bold.input_var, file_path);
    end

    [x, variable_info, region_labels, n_var_region, dx] = read_single_session_roiTs( ...
        S.(bold.input_var), bold);

    if i_sess == 1
        variable_info_ref = variable_info;
        region_labels_ref = region_labels;
        n_var_region_ref = n_var_region;
    elseif size(x, 2) ~= height(variable_info_ref)
        error('Variable count mismatch in %s.', file_path);
    end

    data_cells{i_sess} = x;
    session_lengths(i_sess) = size(x, 1);
    session_dx(i_sess) = dx;
    source_files{i_sess} = file_path;
end

D = struct();
D.data = cat(1, data_cells{:});
D.data_storage = 'memory';
D.n_time = size(D.data, 1);
D.session_ids = session_ids(:);
D.session_lengths = session_lengths;
D.session_dx = session_dx;
D.session_end_idx = cumsum(session_lengths);
D.session_start_idx = [1; D.session_end_idx(1:end-1) + 1];
D.border_idx = D.session_end_idx(1:end-1);
D.dx = unique_dx_or_empty(session_dx);
if ~isempty(D.dx)
    D.fs = 1 / D.dx;
else
    D.fs = [];
end
D.selected_variables = (1:height(variable_info_ref)).';
D.variable_labels = variable_info_ref.variable_label;
D.variable_info = variable_info_ref;
D.region_labels = region_labels_ref;
D.n_var_region = n_var_region_ref;
D.coords = [];
D.input_dir = bold.input_dir;
D.input_var = bold.input_var;
D.source_files = source_files;
D.variable_mode = bold.variable_mode;
D.selected_region_names = get_optional_field(bold, 'selected_region_names', {});
D.selected_regions = unique(variable_info_ref.region_idx, 'stable');
D.dataset_id = cfg.dataset_id;
D.file_stem = cfg.file_stem;
end

function session_ids = collect_included_session_ids(cfg)
session_ids = [];
for i_group = 1:numel(cfg.sessions)
    if cfg.sessions(i_group).include
        session_ids = [session_ids; cfg.sessions(i_group).session_id(:)]; %#ok<AGROW>
    end
end
end

function [x, variable_info, region_labels, n_var_region, dx] = read_single_session_roiTs(roiTs, bold)
if isstruct(roiTs) && isfield(roiTs, 'dat')
    x = double(roiTs.dat);
    dx = get_optional_field(roiTs, 'dx', bold.dx);
    n_var = size(x, 2);
    variable_info = table((1:n_var).', (1:n_var).', repmat({''}, n_var, 1), ...
        ones(n_var, 1), cellstr(compose("var%04d", (1:n_var).')), ...
        'VariableNames', {'variable_idx', 'region_idx', 'region_label', ...
        'within_region_idx', 'variable_label'});
    region_labels = {};
    n_var_region = ones(n_var, 1);
    return;
end

if ~iscell(roiTs)
    error('Expected roiTs to be a cell array of ROI structs or a struct with .dat.');
end

n_regions_total = size(roiTs, 2);
region_labels_all = read_roi_labels_from_roiTs(roiTs);
selected_regions = resolve_selected_regions_from_names(region_labels_all, bold, n_regions_total);
region_labels = region_labels_all(selected_regions);

x_cells = cell(numel(selected_regions), 1);
n_var_region = zeros(numel(selected_regions), 1);
var_region_idx = [];
var_region_label = {};
var_within_idx = [];
var_labels = {};
var_idx = 0;
dx = [];

for i_region = 1:numel(selected_regions)
    region_id = selected_regions(i_region);
    R = roiTs{1, region_id};
    if ~isstruct(R) || ~isfield(R, 'dat')
        error('roiTs{1,%d} does not contain field .dat.', region_id);
    end

    dat = double(R.dat);
    if isempty(dx)
        dx = get_optional_field(R, 'dx', bold.dx);
    end

    if strcmpi(bold.variable_mode, 'roi_mean')
        xr = mean(dat, 2, 'omitnan');
    else
        xr = dat;
    end

    x_cells{i_region} = xr;
    n_var_region(i_region) = size(xr, 2);

    for j_var = 1:size(xr, 2)
        var_idx = var_idx + 1;
        var_region_idx(end+1, 1) = region_id; %#ok<AGROW>
        var_region_label{end+1, 1} = region_labels_all{region_id}; %#ok<AGROW>
        var_within_idx(end+1, 1) = j_var; %#ok<AGROW>
        if strcmpi(bold.variable_mode, 'roi_mean')
            var_labels{end+1, 1} = sprintf('%s_mean', region_labels_all{region_id}); %#ok<AGROW>
        else
            var_labels{end+1, 1} = sprintf('%s_v%04d', region_labels_all{region_id}, j_var); %#ok<AGROW>
        end
    end
end

if isempty(dx)
    dx = bold.dx;
end

x = cat(2, x_cells{:});
variable_info = table((1:var_idx).', var_region_idx, var_region_label, ...
    var_within_idx, var_labels, ...
    'VariableNames', {'variable_idx', 'region_idx', 'region_label', ...
    'within_region_idx', 'variable_label'});
end

function labels = read_roi_labels_from_roiTs(roiTs)
n_regions = size(roiTs, 2);
labels = cell(n_regions, 1);
for i_region = 1:n_regions
    R = roiTs{1, i_region};
    labels{i_region} = read_single_roi_label(R, i_region);
end
end

function label = read_single_roi_label(R, i_region)
label = '';
if isstruct(R)
    candidate_fields = {'name', 'Name', 'label', 'Label', 'roi', 'ROI', ...
        'roi_name', 'roiName', 'region', 'Region', 'region_name', ...
        'regionName', 'atlas_name', 'atlasName', 'atlas_label', 'atlasLabel'};
    for i_field = 1:numel(candidate_fields)
        field_name = candidate_fields{i_field};
        if isfield(R, field_name) && ~isempty(R.(field_name))
            value = R.(field_name);
            if ischar(value) || isstring(value)
                label = char(string(value));
                return;
            end
            if iscell(value) && numel(value) == 1 && (ischar(value{1}) || isstring(value{1}))
                label = char(string(value{1}));
                return;
            end
        end
    end
end
label = sprintf('roi%02d', i_region);
end

function selected_regions = resolve_selected_regions_from_names(region_labels_all, bold, n_regions_total)
if ~isfield(bold, 'selected_region_names') || isempty(bold.selected_region_names)
    selected_regions = 1:n_regions_total;
    return;
end

target_names = cellstr(string(bold.selected_region_names(:)));
selected_regions = zeros(1, numel(target_names));
normalized_labels = cellfun(@normalize_roi_label, region_labels_all, 'UniformOutput', false);

for i_target = 1:numel(target_names)
    target = normalize_roi_label(target_names{i_target});
    match_idx = find(strcmp(normalized_labels, target), 1, 'first');
    if isempty(match_idx)
        error(['Could not resolve BOLD ROI name "%s" from roiTs metadata.\n' ...
            'Available ROI labels read from this file are:\n  %s\n' ...
            'Please make sure roiTs ROI structs contain label/name metadata, or update the reader.'], ...
            target_names{i_target}, strjoin(region_labels_all(:).', ', '));
    end
    selected_regions(i_target) = match_idx;
end
end

function label = normalize_roi_label(label)
label = lower(char(string(label)));
label = regexprep(label, '[^a-z0-9]', '');
end

function val = get_optional_field(S, name, fallback)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    val = S.(name);
else
    val = fallback;
end
end

function dx = unique_dx_or_empty(session_dx)
dx0 = session_dx(1);
if all(abs(session_dx - dx0) < 1e-12)
    dx = dx0;
else
    dx = [];
    warning('BOLD sampling interval dx is inconsistent across sessions.');
end
end
