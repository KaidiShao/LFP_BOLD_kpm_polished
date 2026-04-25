function D = load_bold_dataset(cfg, opts)
% Load selected BOLD roiTs sessions and concatenate them in time.
%
% Expected cfg fields to add manually:
%   cfg.bold.input_file        MAT file containing one concatenated roiTs
%   cfg.bold.input_var         variable name, default 'roiTs'
%   cfg.bold.dx                BOLD sampling interval / TR in seconds
%                              or cfg.bold.session_dx for per-session TRs
%
% Optional cfg.bold fields:
%   .selected_regions          region indices to load
%   .region_labels             cell array of region names
%   .variable_mode             'voxels' (default) or 'roi_mean'
%   .session_lengths           fallback if roiTs has no session metadata
%   .session_length_field      optional roiTs field name for session lengths
%   .demean_on_load            false (default)
%
% Output D.data is time x variables. Session metadata is always preserved.

if nargin < 2
    opts = struct();
end

if ~isfield(opts, 'metadata_only') || isempty(opts.metadata_only)
    opts.metadata_only = false;
end

if ~isfield(cfg, 'bold') || ~isstruct(cfg.bold)
    error(['cfg.bold is missing. Add BOLD fields manually, for example: ' ...
        'cfg.bold.input_file and cfg.bold.dx.']);
end

bold = cfg.bold;

if ~isfield(bold, 'input_file') || isempty(bold.input_file)
    error('cfg.bold.input_file must point to the MAT file containing BOLD roiTs data.');
end
if ~isfield(bold, 'input_var') || isempty(bold.input_var)
    bold.input_var = 'roiTs';
end
if ~isfield(bold, 'dx') || isempty(bold.dx)
    if isfield(bold, 'session_dx') && ~isempty(bold.session_dx)
        bold.dx = [];
    else
        error('cfg.bold.dx or cfg.bold.session_dx must be set to the BOLD sampling interval / TR in seconds.');
    end
end
if ~isfield(bold, 'variable_mode') || isempty(bold.variable_mode)
    bold.variable_mode = 'voxels';
end
if ~isfield(bold, 'demean_on_load') || isempty(bold.demean_on_load)
    bold.demean_on_load = false;
end

session_ids = collect_included_session_ids(cfg);

if isempty(session_ids)
    error('No included sessions were found in cfg.sessions.');
end

fprintf('Loading BOLD source:\n  %s\n', bold.input_file);
S = load(bold.input_file, bold.input_var);
if ~isfield(S, bold.input_var)
    error('Variable "%s" was not found in %s.', bold.input_var, bold.input_file);
end

src = S.(bold.input_var);
[session_lengths_all, session_ids_all] = read_roiTs_session_metadata(src, bold);
session_lengths = select_session_lengths(session_lengths_all, session_ids_all, session_ids, bold);

data_cells = {};
session_dx = resolve_session_dx(bold, numel(session_ids));

fprintf('Loading single concatenated BOLD roiTs object ...\n');

[x_all, variable_info_ref, region_labels_ref, n_var_region_ref, coords_ref] = ...
    read_roiTs_object(src, bold);

if bold.demean_on_load
    idx_end = cumsum(session_lengths);
    idx_start = [1; idx_end(1:end-1) + 1];
    for k = 1:numel(session_ids)
        idx = idx_start(k):idx_end(k);
        x_all(idx, :) = x_all(idx, :) - mean(x_all(idx, :), 1, 'omitnan');
    end
end

if size(x_all, 1) ~= sum(session_lengths)
    error('Concatenated roiTs length (%d) does not match sum(cfg.bold.session_lengths) (%d).', ...
        size(x_all, 1), sum(session_lengths));
end

if ~opts.metadata_only
    data_cells{1, 1} = x_all;
end

fprintf('         done, size = [%d, %d]\n', size(x_all, 1), size(x_all, 2));

D = struct();
if opts.metadata_only
    D.data = [];
    D.data_storage = 'source';
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
D.input_file = bold.input_file;
D.input_var = bold.input_var;
D.variable_mode = bold.variable_mode;
if isfield(cfg, 'dataset_id')
    D.dataset_id = cfg.dataset_id;
end
if isfield(cfg, 'file_stem')
    D.file_stem = cfg.file_stem;
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

function session_dx = resolve_session_dx(bold, n_sessions)
if isfield(bold, 'session_dx') && ~isempty(bold.session_dx)
    session_dx = double(bold.session_dx(:));
    if isscalar(session_dx)
        session_dx = repmat(session_dx, n_sessions, 1);
    elseif numel(session_dx) ~= n_sessions
        error('cfg.bold.session_dx must be scalar or match the number of included sessions.');
    end
else
    session_dx = repmat(double(bold.dx), n_sessions, 1);
end
end

function [session_lengths, session_ids] = read_roiTs_session_metadata(roiTs, bold)
session_lengths = [];
session_ids = [];

candidate_length_fields = {'session_lengths', 'session_length', ...
    'T_sess', 'sess_lengths', 'run_lengths'};
candidate_id_fields = {'session_ids', 'session_id', 'sess_ids', 'sess_id'};

if isfield(bold, 'session_length_field') && ~isempty(bold.session_length_field)
    candidate_length_fields = [{bold.session_length_field}, candidate_length_fields];
end
if isfield(bold, 'session_id_field') && ~isempty(bold.session_id_field)
    candidate_id_fields = [{bold.session_id_field}, candidate_id_fields];
end

if isstruct(roiTs)
    for i = 1:numel(candidate_length_fields)
        name = candidate_length_fields{i};
        if isfield(roiTs, name) && ~isempty(roiTs.(name))
            session_lengths = double(roiTs.(name)(:));
            break;
        end
    end

    for i = 1:numel(candidate_id_fields)
        name = candidate_id_fields{i};
        if isfield(roiTs, name) && ~isempty(roiTs.(name))
            session_ids = double(roiTs.(name)(:));
            break;
        end
    end
end

if isempty(session_lengths) && iscell(roiTs) && ~isempty(roiTs)
    R = roiTs{1};
    if isstruct(R)
        for i = 1:numel(candidate_length_fields)
            name = candidate_length_fields{i};
            if isfield(R, name) && ~isempty(R.(name))
                session_lengths = double(R.(name)(:));
                break;
            end
        end

        for i = 1:numel(candidate_id_fields)
            name = candidate_id_fields{i};
            if isfield(R, name) && ~isempty(R.(name))
                session_ids = double(R.(name)(:));
                break;
            end
        end
    end
end
end

function session_lengths = select_session_lengths(session_lengths_all, session_ids_all, session_ids, bold)
if isempty(session_lengths_all)
    if isfield(bold, 'session_lengths') && ~isempty(bold.session_lengths)
        session_lengths_all = double(bold.session_lengths(:));
    else
        error(['Could not find session lengths in roiTs. Add a session-length field to roiTs ' ...
            'or set cfg.bold.session_lengths as a fallback.']);
    end
end

if ~isempty(session_ids_all)
    session_lengths = zeros(numel(session_ids), 1);
    for k = 1:numel(session_ids)
        idx = find(session_ids_all == session_ids(k), 1, 'first');
        if isempty(idx)
            error('Selected session id %g was not found in roiTs session metadata.', session_ids(k));
        end
        session_lengths(k) = session_lengths_all(idx);
    end
else
    if numel(session_lengths_all) ~= numel(session_ids)
        error(['roiTs session lengths do not include session IDs, so their count must match ' ...
            'the number of included cfg sessions.']);
    end
    session_lengths = session_lengths_all(:);
end
end

function [x, variable_info, region_labels, n_var_region, coords_all] = read_roiTs_object(roiTs, bold)
if isstruct(roiTs) && isfield(roiTs, 'observables') && ...
        strcmpi(bold.variable_mode, 'observables')
    x = double(roiTs.observables);
    n_var = size(x, 2);
    variable_info = table((1:n_var).', repmat({''}, n_var, 1), ...
        (1:n_var).', cellstr(compose("obs%04d", (1:n_var).')), ...
        'VariableNames', {'variable_idx', 'region_label', 'within_region_idx', 'variable_label'});
    region_labels = {};
    n_var_region = n_var;
    coords_all = [];
    return;
end

if isstruct(roiTs) && isfield(roiTs, 'dat') && ~iscell(roiTs)
    x = double(roiTs.dat);
    n_var = size(x, 2);
    variable_info = table((1:n_var).', repmat({''}, n_var, 1), ...
        (1:n_var).', cellstr(compose("var%04d", (1:n_var).')), ...
        'VariableNames', {'variable_idx', 'region_label', 'within_region_idx', 'variable_label'});
    region_labels = {};
    n_var_region = n_var;
    coords_all = get_optional_field(roiTs, 'coords', []);
    return;
end

if ~iscell(roiTs)
    error('roiTs must be a struct with .dat/.observables or a cell array of region structs.');
end

n_regions_total = size(roiTs, 2);
if isfield(bold, 'selected_regions') && ~isempty(bold.selected_regions)
    selected_regions = bold.selected_regions(:).';
else
    selected_regions = 1:n_regions_total;
end

region_labels_all = resolve_region_labels(bold, n_regions_total);
region_labels = region_labels_all(selected_regions);

x_cells = cell(numel(selected_regions), 1);
n_var_region = zeros(numel(selected_regions), 1);
coords_cells = cell(numel(selected_regions), 1);
var_region_idx = [];
var_region_label = {};
var_within_idx = [];
var_labels = {};
var_idx = 0;

for i = 1:numel(selected_regions)
    r = selected_regions(i);
    R = roiTs{1, r};
    if ~isfield(R, 'dat')
        error('roiTs{1,%d} does not contain field .dat.', r);
    end

    dat = double(R.dat);
    if strcmpi(bold.variable_mode, 'roi_mean')
        xr = mean(dat, 2, 'omitnan');
    else
        xr = dat;
    end

    x_cells{i} = xr;
    n_var_region(i) = size(xr, 2);
    coords_cells{i} = get_optional_field(R, 'coords', []);

    for j = 1:size(xr, 2)
        var_idx = var_idx + 1;
        var_region_idx(end+1, 1) = r; %#ok<AGROW>
        var_region_label{end+1, 1} = region_labels_all{r}; %#ok<AGROW>
        var_within_idx(end+1, 1) = j; %#ok<AGROW>
        if strcmpi(bold.variable_mode, 'roi_mean')
            var_labels{end+1, 1} = sprintf('%s_mean', region_labels_all{r}); %#ok<AGROW>
        else
            var_labels{end+1, 1} = sprintf('%s_v%04d', region_labels_all{r}, j); %#ok<AGROW>
        end
    end
end

x = cat(2, x_cells{:});
coords_all = cat_coords(coords_cells, n_var_region, bold.variable_mode);
variable_info = table((1:var_idx).', var_region_idx, var_region_label, ...
    var_within_idx, var_labels, ...
    'VariableNames', {'variable_idx', 'region_idx', 'region_label', ...
    'within_region_idx', 'variable_label'});
end

function labels = resolve_region_labels(bold, n_regions)
if isfield(bold, 'region_labels') && ~isempty(bold.region_labels)
    labels = cellstr(string(bold.region_labels(:)));
    if numel(labels) < n_regions
        error('cfg.bold.region_labels has fewer entries than roiTs regions.');
    end
else
    labels = cell(n_regions, 1);
    for i = 1:n_regions
        labels{i} = sprintf('roi%02d', i);
    end
end
end

function val = get_optional_field(S, name, fallback)
if isstruct(S) && isfield(S, name)
    val = S.(name);
else
    val = fallback;
end
end

function coords_all = cat_coords(coords_cells, n_var_region, variable_mode)
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
        coords_all = [];
        return;
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
