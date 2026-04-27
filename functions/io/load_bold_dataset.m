function D = load_bold_dataset(cfg, opts)
%LOAD_BOLD_DATASET Load selected raw BOLD roiTs sessions and concatenate them in time.
%
% Expected cfg fields:
%   cfg.raw_data_root
%   cfg.file_stem
%   cfg.sessions(i).session_id
%   cfg.sessions(i).include
%
% Optional cfg.bold fields:
%   .data_subfolder       default 'roits'
%   .input_var            default 'roiTs'
%   .variable_mode        'voxels' (default) or 'roi_mean'
%   .selected_region_names exact ROI names from roiTs{1,r}.name
%
% Formal input contract:
%   - each session file is <file_stem>_<session_id>_roits.mat
%   - files live under fullfile(cfg.raw_data_root, cfg.bold.data_subfolder)
%   - each file contains variable cfg.bold.input_var
%   - roiTs must be a cell array of ROI structs
%   - each selected ROI struct must contain fields .dat, .dx, and .name
%
% Output D.data is time x variables. Session metadata is always preserved.

if nargin < 2
    opts = struct();
end

if ~isfield(opts, 'metadata_only') || isempty(opts.metadata_only)
    opts.metadata_only = false;
end

bold = local_resolve_bold_cfg(cfg);
session_ids = collect_included_session_ids(cfg);

if isempty(session_ids)
    error('No included sessions were found in cfg.sessions.');
end

data_dir = fullfile(cfg.raw_data_root, bold.data_subfolder);
if exist(data_dir, 'dir') ~= 7
    error('BOLD data_subfolder does not exist:\n  %s', data_dir);
end

if ~opts.metadata_only
    data_cells = cell(numel(session_ids), 1);
end

session_lengths = zeros(numel(session_ids), 1);
session_dx = zeros(numel(session_ids), 1);
source_files = cell(numel(session_ids), 1);

variable_info_ref = table();
region_labels_ref = {};
n_var_region_ref = [];
coords_ref = [];
selected_regions_ref = [];
all_region_names_ref = {};

for i_sess = 1:numel(session_ids)
    session_id = session_ids(i_sess);
    file_path = fullfile(data_dir, sprintf('%s_%04d_roits.mat', ...
        cfg.file_stem, session_id));
    if exist(file_path, 'file') ~= 2
        error('Missing BOLD roits file for session %g:\n  %s', session_id, file_path);
    end

    fprintf('Loading BOLD session %g:\n  %s\n', session_id, file_path);
    S = load(file_path, bold.input_var);
    if ~isfield(S, bold.input_var)
        error('Variable "%s" was not found in %s.', bold.input_var, file_path);
    end

    [x, variable_info, region_labels, n_var_region, coords, dx, ...
        selected_regions, all_region_names] = ...
        read_single_session_roiTs(S.(bold.input_var), bold, file_path);

    if i_sess == 1
        variable_info_ref = variable_info;
        region_labels_ref = region_labels;
        n_var_region_ref = n_var_region;
        coords_ref = coords;
        selected_regions_ref = selected_regions;
        all_region_names_ref = all_region_names;
    else
        if ~isequal(variable_info_ref, variable_info)
            error(['BOLD ROI variable layout changed across sessions.\n' ...
                'Reference session file:\n  %s\n' ...
                'Mismatched session file:\n  %s'], source_files{1}, file_path);
        end
        if ~isequal(region_labels_ref, region_labels) || ...
                ~isequal(n_var_region_ref, n_var_region) || ...
                ~isequal(selected_regions_ref, selected_regions)
            error(['Selected BOLD ROI structure changed across sessions.\n' ...
                'Reference session file:\n  %s\n' ...
                'Mismatched session file:\n  %s'], source_files{1}, file_path);
        end
        if ~isequal(all_region_names_ref, all_region_names)
            error(['Full BOLD ROI name list changed across sessions.\n' ...
                'Reference session file:\n  %s\n' ...
                'Mismatched session file:\n  %s'], source_files{1}, file_path);
        end
    end

    if ~opts.metadata_only
        data_cells{i_sess} = x;
    end

    session_lengths(i_sess) = size(x, 1);
    session_dx(i_sess) = dx;
    source_files{i_sess} = file_path;
end

D = struct();
if opts.metadata_only
    D.data = [];
    D.data_storage = 'metadata_only';
    D.n_time = sum(session_lengths);
else
    D.data = cat(1, data_cells{:});
    D.data_storage = 'memory';
    D.n_time = size(D.data, 1);
end

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
D.coords = coords_ref;
D.input_dir = data_dir;
D.input_var = bold.input_var;
D.source_files = source_files;
D.variable_mode = bold.variable_mode;
D.selected_region_names = bold.selected_region_names;
D.selected_regions = selected_regions_ref(:);
D.all_region_names = all_region_names_ref(:);
D.roi_label_field = 'name';
if isfield(cfg, 'dataset_id')
    D.dataset_id = cfg.dataset_id;
end
if isfield(cfg, 'file_stem')
    D.file_stem = cfg.file_stem;
end
end

function bold = local_resolve_bold_cfg(cfg)
if isfield(cfg, 'bold') && isstruct(cfg.bold)
    bold = cfg.bold;
else
    bold = struct();
end

if ~isfield(bold, 'data_subfolder') || isempty(bold.data_subfolder)
    bold.data_subfolder = 'roits';
end
if ~isfield(bold, 'input_var') || isempty(bold.input_var)
    bold.input_var = 'roiTs';
end
if ~isfield(bold, 'variable_mode') || isempty(bold.variable_mode)
    bold.variable_mode = 'voxels';
end
if ~isfield(bold, 'selected_region_names') || isempty(bold.selected_region_names)
    bold.selected_region_names = {};
else
    bold.selected_region_names = cellstr(string(bold.selected_region_names(:)).');
end

switch lower(bold.variable_mode)
    case {'voxels', 'roi_mean'}
        % supported
    otherwise
        error('cfg.bold.variable_mode must be ''voxels'' or ''roi_mean''.');
end
end

function session_ids = collect_included_session_ids(cfg)
session_ids = [];
for i = 1:numel(cfg.sessions)
    if cfg.sessions(i).include
        session_ids = [session_ids; cfg.sessions(i).session_id(:)]; %#ok<AGROW>
    end
end
end

function [x, variable_info, region_labels, n_var_region, coords_all, dx, ...
    selected_regions, all_region_names] = read_single_session_roiTs(roiTs, bold, file_path)
if ~iscell(roiTs) || isempty(roiTs)
    error(['BOLD roiTs must be a nonempty cell array of ROI structs.\n' ...
        'File:\n  %s'], file_path);
end

n_regions_total = size(roiTs, 2);
all_region_names = read_roi_names_from_roiTs(roiTs, file_path);
selected_regions = resolve_selected_regions_from_names(all_region_names, bold, file_path);
region_labels = all_region_names(selected_regions);

x_cells = cell(numel(selected_regions), 1);
n_var_region = zeros(numel(selected_regions), 1);
coords_cells = cell(numel(selected_regions), 1);
var_region_idx = [];
var_region_label = {};
var_within_idx = [];
var_labels = {};
var_idx = 0;
dx = [];

for i_region = 1:numel(selected_regions)
    region_id = selected_regions(i_region);
    R = roiTs{1, region_id};
    if ~isstruct(R)
        error('roiTs{1,%d} is not a struct in file:\n  %s', region_id, file_path);
    end
    if ~isfield(R, 'dat') || isempty(R.dat)
        error('roiTs{1,%d} does not contain nonempty field .dat in file:\n  %s', ...
            region_id, file_path);
    end
    if ~isfield(R, 'dx') || isempty(R.dx)
        error('roiTs{1,%d} does not contain nonempty field .dx in file:\n  %s', ...
            region_id, file_path);
    end

    dat = double(R.dat);
    region_dx = double(R.dx);
    if ~isscalar(region_dx) || ~isfinite(region_dx) || region_dx <= 0
        error('roiTs{1,%d}.dx must be a positive finite scalar in file:\n  %s', ...
            region_id, file_path);
    end
    if isempty(dx)
        dx = region_dx;
    elseif abs(region_dx - dx) > 1e-12
        error(['Selected BOLD ROIs disagree on .dx within one session file.\n' ...
            'File:\n  %s\nFirst dx = %.12g, roiTs{1,%d}.dx = %.12g'], ...
            file_path, dx, region_id, region_dx);
    end

    if strcmpi(bold.variable_mode, 'roi_mean')
        xr = mean(dat, 2, 'omitnan');
    else
        xr = dat;
    end

    x_cells{i_region} = xr;
    n_var_region(i_region) = size(xr, 2);

    if isfield(R, 'coords') && ~isempty(R.coords)
        coords_cells{i_region} = double(R.coords);
    else
        coords_cells{i_region} = [];
    end

    for j_var = 1:size(xr, 2)
        var_idx = var_idx + 1;
        var_region_idx(end+1, 1) = region_id; %#ok<AGROW>
        var_region_label{end+1, 1} = all_region_names{region_id}; %#ok<AGROW>
        var_within_idx(end+1, 1) = j_var; %#ok<AGROW>
        if strcmpi(bold.variable_mode, 'roi_mean')
            var_labels{end+1, 1} = sprintf('%s_mean', all_region_names{region_id}); %#ok<AGROW>
        else
            var_labels{end+1, 1} = sprintf('%s_v%04d', all_region_names{region_id}, j_var); %#ok<AGROW>
        end
    end
end

x = cat(2, x_cells{:});
coords_all = cat_coords(coords_cells, n_var_region, bold.variable_mode, file_path, selected_regions);
variable_info = table((1:var_idx).', var_region_idx, var_region_label, ...
    var_within_idx, var_labels, ...
    'VariableNames', {'variable_idx', 'region_idx', 'region_label', ...
    'within_region_idx', 'variable_label'});
end

function labels = read_roi_names_from_roiTs(roiTs, file_path)
n_regions = size(roiTs, 2);
labels = cell(n_regions, 1);
for i_region = 1:n_regions
    R = roiTs{1, i_region};
    if ~isstruct(R)
        error('roiTs{1,%d} is not a struct in file:\n  %s', i_region, file_path);
    end
    labels{i_region} = read_single_roi_name(R, i_region, file_path);
end
end

function label = read_single_roi_name(R, i_region, file_path)
if ~isfield(R, 'name') || isempty(R.name)
    error(['roiTs{1,%d} is missing required ROI label field .name.\n' ...
        'File:\n  %s'], i_region, file_path);
end

value = R.name;
if ischar(value) || isstring(value)
    label = char(string(value));
    return;
end
if iscell(value) && numel(value) == 1 && (ischar(value{1}) || isstring(value{1}))
    label = char(string(value{1}));
    return;
end

error(['roiTs{1,%d}.name must be a string-like scalar.\n' ...
    'File:\n  %s'], i_region, file_path);
end

function selected_regions = resolve_selected_regions_from_names(region_labels_all, bold, file_path)
if isempty(bold.selected_region_names)
    selected_regions = 1:numel(region_labels_all);
    return;
end

target_names = bold.selected_region_names;
selected_regions = zeros(1, numel(target_names));

for i_target = 1:numel(target_names)
    match_idx = find(strcmp(region_labels_all, target_names{i_target}), 1, 'first');
    if isempty(match_idx)
        error(['Could not resolve BOLD ROI name "%s" from roiTs{1,r}.name.\n' ...
            'File:\n  %s\nAvailable ROI names:\n  %s'], ...
            target_names{i_target}, file_path, strjoin(region_labels_all(:).', ', '));
    end
    selected_regions(i_target) = match_idx;
end
end

function coords_all = cat_coords(coords_cells, n_var_region, variable_mode, file_path, selected_regions)
coords_all = [];
if strcmpi(variable_mode, 'roi_mean')
    return;
end

for i = 1:numel(coords_cells)
    c = coords_cells{i};
    if isempty(c)
        coords_all = [];
        return;
    end
    if size(c, 1) ~= n_var_region(i)
        error(['roiTs{1,%d}.coords row count does not match the number of voxel variables.\n' ...
            'File:\n  %s'], selected_regions(i), file_path);
    end
    coords_all = [coords_all; c]; %#ok<AGROW>
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
