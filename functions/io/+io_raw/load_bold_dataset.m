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
%   .data_subfolder        default 'roits'
%   .input_var             default 'roiTs'
%   .variable_mode         'voxels' (default) or 'roi_mean'
%   .selected_region_names exact roiTs{1,r}.name values to keep
%
% Each session file is expected to be:
%   <file_stem>_<session_id>_roits.mat
%
% Output D.data is time x variables. Session metadata is always preserved.

if nargin < 2
    opts = struct();
end

if ~isfield(opts, 'metadata_only') || isempty(opts.metadata_only)
    opts.metadata_only = false;
end

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

data_dir = fullfile(cfg.raw_data_root, bold.data_subfolder);
if exist(data_dir, 'dir') ~= 7
    error('BOLD data_subfolder does not exist:\n  %s', data_dir);
end

session_ids = io_utils.collect_included_session_entries(cfg);

if isempty(session_ids)
    error('No included sessions were found in cfg.sessions.');
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

    fprintf('[%d/%d] Loading BOLD session %04d ...\n', ...
        i_sess, numel(session_ids), session_id);

    S = load(file_path, bold.input_var);
    if ~isfield(S, bold.input_var)
        error('Variable "%s" was not found in %s.', bold.input_var, file_path);
    end

    [x, variable_info, region_labels, n_var_region, coords, dx, ...
        selected_regions, all_region_names] = local_read_roiTs_session( ...
        S.(bold.input_var), bold.variable_mode, bold.selected_region_names, file_path);

    if isempty(variable_info_ref)
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

        if strcmpi(bold.variable_mode, 'voxels') && ...
                ~local_coords_match(coords_ref, coords)
            error(['BOLD voxel coordinates changed across sessions.\n' ...
                'Reference session file:\n  %s\n' ...
                'Mismatched session file:\n  %s\n' ...
                'For variable_mode = ''voxels'', ROI voxel coordinates must be identical ' ...
                'across sessions.'], source_files{1}, file_path);
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

    fprintf('         done, size = [%d, %d], dx = %.12g\n', ...
        size(x, 1), size(x, 2), dx);
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
[D.session_start_idx, D.session_end_idx, D.border_idx] = ...
    io_utils.build_session_index_metadata(session_lengths);
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

[D.dx, D.fs] = io_utils.resolve_uniform_dx(session_dx, ...
    'BOLD sampling interval dx is inconsistent across sessions.');
end

function [x, variable_info, region_labels, n_var_region, coords_all, dx, ...
    selected_regions, all_region_names] = local_read_roiTs_session( ...
    roiTs, variable_mode, selected_region_names, file_path)
if ~iscell(roiTs) || isempty(roiTs)
    error(['BOLD roiTs must be a nonempty cell array of ROI structs.\n' ...
        'File:\n  %s'], file_path);
end

all_region_names = local_read_all_roi_names(roiTs, file_path);
selected_regions = local_resolve_selected_regions( ...
    all_region_names, selected_region_names, file_path);
region_labels = all_region_names(selected_regions);
[x, variable_info, n_var_region, coords_all, dx] = ...
    local_extract_selected_region_data( ...
        roiTs, selected_regions, all_region_names, variable_mode, file_path);
end

function all_region_names = local_read_all_roi_names(roiTs, file_path)
n_regions = size(roiTs, 2);
all_region_names = cell(n_regions, 1);

for i_region = 1:n_regions
    R = roiTs{1, i_region};
    if ~isstruct(R)
        error('roiTs{1,%d} is not a struct in file:\n  %s', i_region, file_path);
    end
    all_region_names{i_region} = local_read_single_roi_name(R, i_region, file_path);
end
end

function selected_regions = local_resolve_selected_regions( ...
    all_region_names, selected_region_names, file_path)
if isempty(selected_region_names)
    selected_regions = 1:numel(all_region_names);
    return;
end

selected_regions = zeros(1, numel(selected_region_names));
for i_target = 1:numel(selected_region_names)
    match_idx = find(strcmp(all_region_names, selected_region_names{i_target}), 1, 'first');
    if isempty(match_idx)
        error(['Could not resolve BOLD ROI name "%s" from roiTs{1,r}.name.\n' ...
            'File:\n  %s\nAvailable ROI names:\n  %s'], ...
            selected_region_names{i_target}, file_path, strjoin(all_region_names(:).', ', '));
    end
    selected_regions(i_target) = match_idx;
end
end

function [x, variable_info, n_var_region, coords_all, dx] = ...
    local_extract_selected_region_data( ...
        roiTs, selected_regions, all_region_names, variable_mode, file_path)
x_cells = cell(numel(selected_regions), 1);
coords_cells = cell(numel(selected_regions), 1);
n_var_region = zeros(numel(selected_regions), 1);
var_region_idx_cells = cell(numel(selected_regions), 1);
var_region_label_cells = cell(numel(selected_regions), 1);
var_within_idx_cells = cell(numel(selected_regions), 1);
var_label_cells = cell(numel(selected_regions), 1);
dx = [];

for i_region = 1:numel(selected_regions)
    region_id = selected_regions(i_region);
    region_name = all_region_names{region_id};
    R = roiTs{1, region_id};

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

    if strcmpi(variable_mode, 'roi_mean')
        xr = mean(dat, 2, 'omitnan');
        coords_cells{i_region} = [];
        var_label_cells{i_region} = {sprintf('%s_mean', region_name)};
    else
        xr = dat;
        if isfield(R, 'coords') && ~isempty(R.coords)
            coords_cells{i_region} = double(R.coords);
            if size(coords_cells{i_region}, 1) ~= size(xr, 2)
                error(['roiTs{1,%d}.coords row count does not match the number of voxel variables.\n' ...
                    'File:\n  %s'], region_id, file_path);
            end
        else
            coords_cells{i_region} = [];
        end
        var_label_cells{i_region} = arrayfun( ...
            @(j_var) sprintf('%s_v%04d', region_name, j_var), ...
            (1:size(xr, 2)).', 'UniformOutput', false);
    end

    x_cells{i_region} = xr;
    n_var_region(i_region) = size(xr, 2);
    var_region_idx_cells{i_region} = repmat(region_id, size(xr, 2), 1);
    var_region_label_cells{i_region} = repmat({region_name}, size(xr, 2), 1);
    var_within_idx_cells{i_region} = (1:size(xr, 2)).';
end

x = cat(2, x_cells{:});
coords_all = local_cat_region_coords(coords_cells, variable_mode);

variable_info = table( ...
    (1:sum(n_var_region)).', ...
    cat(1, var_region_idx_cells{:}), ...
    cat(1, var_region_label_cells{:}), ...
    cat(1, var_within_idx_cells{:}), ...
    cat(1, var_label_cells{:}), ...
    'VariableNames', {'variable_idx', 'region_idx', 'region_label', ...
    'within_region_idx', 'variable_label'});
end

function coords_all = local_cat_region_coords(coords_cells, variable_mode)
if strcmpi(variable_mode, 'roi_mean')
    coords_all = [];
    return;
end

for i_region = 1:numel(coords_cells)
    if isempty(coords_cells{i_region})
        coords_all = [];
        return;
    end
end

coords_all = cat(1, coords_cells{:});
end

function label = local_read_single_roi_name(R, i_region, file_path)
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

function tf = local_coords_match(coords_ref, coords)
if isempty(coords_ref) && isempty(coords)
    tf = true;
    return;
end

if xor(isempty(coords_ref), isempty(coords))
    tf = false;
    return;
end

if ~isequal(size(coords_ref), size(coords))
    tf = false;
    return;
end

same_nan = isnan(coords_ref) & isnan(coords);
same_val = abs(coords_ref - coords) < 1e-12;
tf = all(same_nan(:) | same_val(:));
end
