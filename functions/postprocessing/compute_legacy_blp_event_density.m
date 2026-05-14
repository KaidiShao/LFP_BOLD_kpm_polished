function E = compute_legacy_blp_event_density(cfg, output_root, source_event_file, detector_name, params)
%COMPUTE_LEGACY_BLP_EVENT_DENSITY Compute density for legacy PL event files.
%
% This adapter treats legacy event files such as mevt_pl and sevt_pl as
% parallel event detectors. It bins their final event onsets by the legacy
% split labels instead of trying to coerce them into pipeline2 DetectResults.

if nargin < 2 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end
if nargin < 3 || isempty(source_event_file)
    error('source_event_file is required.');
end
if nargin < 4 || isempty(detector_name)
    detector_name = infer_detector_name(source_event_file);
end
if nargin < 5
    params = struct();
end

params = apply_default_params(params);
detector_name = char(string(detector_name));

if exist(source_event_file, 'file') ~= 2
    error('Legacy event file was not found: %s', source_event_file);
end

S_file = load(source_event_file);
legacy_var = resolve_legacy_variable(S_file, params.legacy_variable, detector_name);
L = S_file.(legacy_var);

[onset_sec, split_idx] = extract_onsets_and_splits(L);
[event_labels, event_ranges] = extract_event_labels_and_ranges(L, split_idx);
n_types = numel(event_labels);

t_end = resolve_duration_sec(L, onset_sec, params.bin_sec);
edges = 0:params.bin_sec:t_end;
if edges(end) < t_end
    edges = [edges, t_end];
end
t_centers = edges(1:end-1) + diff(edges) / 2;
n_bins = numel(t_centers);

counts_by_type = zeros(n_bins, n_types);
density_by_type = zeros(n_bins, n_types);
smoothed_density_by_type = zeros(n_bins, n_types);
total_event_count = zeros(n_types, 1);

g = build_gaussian_kernel(params.smooth_sigma_sec, params.bin_sec);

for k = 1:n_types
    this_onset = onset_sec(split_idx == k);
    this_onset = this_onset(isfinite(this_onset) & this_onset >= 0 & this_onset <= t_end);
    counts = histcounts(this_onset, edges);
    dens = counts(:) ./ diff(edges(:));

    counts_by_type(:, k) = counts(:);
    density_by_type(:, k) = dens;
    smoothed_density_by_type(:, k) = conv(dens, g, 'same');
    total_event_count(k) = numel(this_onset);
end

save_dir = io_project.get_pipeline_stage_dir(output_root, cfg, 2, 'legacy_event_density');
if exist(save_dir, 'dir') ~= 7
    mkdir(save_dir);
end

group_name = char(string(params.group_name));
if isempty(group_name)
    group_name = infer_group_name(source_event_file);
end

save_tag = build_save_tag(cfg.file_stem, group_name, detector_name, params.bin_sec);
save_file = fullfile(save_dir, [save_tag, '.mat']);
source_event_file_signature = io_utils.build_file_signature(source_event_file);

E = struct();
E.save_file = save_file;
E.source_event_file = source_event_file;
E.source_event_file_signature = source_event_file_signature;
E.dataset_id = cfg.dataset_id;
E.file_stem = cfg.file_stem;
E.pipeline_stage = 'pipeline2_legacy_event_density';
E.detector_name = detector_name;
E.legacy_variable = legacy_var;
E.group_name = group_name;

E.bin_sec = params.bin_sec;
E.bin_hz = 1 / params.bin_sec;
E.smooth_sigma_sec = params.smooth_sigma_sec;
if params.smooth_sigma_sec > 0
    E.smooth_sigma_hz = 1 / params.smooth_sigma_sec;
else
    E.smooth_sigma_hz = [];
end

E.t_edges = edges(:);
E.t_centers = t_centers(:);
E.t_end = t_end;
E.event_labels = event_labels(:);
E.event_ranges = event_ranges(:);
E.onset_sec = onset_sec(:);
E.split_idx = split_idx(:);
E.counts_by_type = counts_by_type;
E.density_by_type = density_by_type;
E.smoothed_density_by_type = smoothed_density_by_type;
E.total_event_count = total_event_count;
E.legacy_meta = extract_legacy_metadata(L);
E.params = params;
E.types = build_type_summaries(event_labels, event_ranges, edges, t_centers, ...
    counts_by_type, density_by_type, smoothed_density_by_type, total_event_count);

save(save_file, 'E', '-v7.3');
end


function params = apply_default_params(params)
if ~isfield(params, 'bin_sec')
    params.bin_sec = 2;
end
if ~isfield(params, 'smooth_sigma_sec')
    params.smooth_sigma_sec = params.bin_sec;
end
if ~isfield(params, 'legacy_variable')
    params.legacy_variable = '';
end
if ~isfield(params, 'group_name')
    params.group_name = '';
end

if ~isscalar(params.bin_sec) || params.bin_sec <= 0
    error('params.bin_sec must be a positive scalar.');
end
if ~isscalar(params.smooth_sigma_sec) || params.smooth_sigma_sec < 0
    error('params.smooth_sigma_sec must be a nonnegative scalar.');
end
end


function detector_name = infer_detector_name(source_event_file)
[~, detector_name] = fileparts(source_event_file);
detector_name = regexprep(detector_name, '^e10gb1_spont1?_', '');
detector_name = regexprep(detector_name, '^e10gb1_', '');
end


function group_name = infer_group_name(source_event_file)
parent_dir = fileparts(source_event_file);
[~, group_name] = fileparts(parent_dir);
end


function legacy_var = resolve_legacy_variable(S_file, requested_var, detector_name)
if ~isempty(requested_var)
    legacy_var = char(string(requested_var));
    if ~isfield(S_file, legacy_var)
        error('Requested legacy variable "%s" was not found in the source file.', legacy_var);
    end
    return;
end

names = fieldnames(S_file);
names = names(~strncmp(names, '__', 2));
if ismember(detector_name, names)
    legacy_var = detector_name;
elseif numel(names) == 1
    legacy_var = names{1};
else
    error('Could not infer the legacy event variable. Candidates: %s', strjoin(names, ', '));
end
end


function [onset_sec, split_idx] = extract_onsets_and_splits(L)
if ~isfield(L, 'onset') || isempty(L.onset)
    error('Legacy event struct must contain a nonempty onset field.');
end
if ~isfield(L, 'split') || isempty(L.split)
    error('Legacy event struct must contain a nonempty split field.');
end

onset_sec = double(L.onset(:));
split_idx = double(L.split(:));

if numel(onset_sec) ~= numel(split_idx)
    error('Legacy onset and split must have the same number of elements.');
end

valid = isfinite(onset_sec) & isfinite(split_idx) & split_idx == round(split_idx) & split_idx >= 1;
onset_sec = onset_sec(valid);
split_idx = split_idx(valid);
end


function [event_labels, event_ranges] = extract_event_labels_and_ranges(L, split_idx)
n_types = max(split_idx);

if isfield(L, 'bname') && ~isempty(L.bname)
    event_labels = normalize_cellstr(L.bname);
else
    event_labels = cell(n_types, 1);
    for k = 1:n_types
        event_labels{k} = sprintf('event_type_%02d', k);
    end
end

if numel(event_labels) < n_types
    for k = numel(event_labels)+1:n_types
        event_labels{k, 1} = sprintf('event_type_%02d', k);
    end
else
    event_labels = event_labels(1:n_types);
end

if isfield(L, 'brange') && ~isempty(L.brange)
    event_ranges = normalize_range_list(L.brange, n_types);
else
    event_ranges = cell(n_types, 1);
    event_ranges(:) = {[]};
end
end


function t_end = resolve_duration_sec(L, onset_sec, bin_sec)
if isfield(L, 'Duration') && ~isempty(L.Duration) && isfinite(double(L.Duration))
    t_end = double(L.Duration);
else
    t_end = ceil(max(onset_sec) / bin_sec) * bin_sec;
end
if isempty(t_end) || ~isfinite(t_end) || t_end <= 0
    error('Could not resolve a positive legacy recording duration.');
end
end


function meta = extract_legacy_metadata(L)
meta = struct();
fields = {'EleSite', 'EvtField', 'grpname', 'session', 'Duration', 'evtnum', 'evtdensity'};
for i = 1:numel(fields)
    f = fields{i};
    if isfield(L, f)
        meta.(f) = L.(f);
    end
end
end


function out = normalize_cellstr(value)
if iscell(value)
    out = cellfun(@char, value(:), 'UniformOutput', false);
elseif isstring(value)
    out = cellstr(value(:));
elseif ischar(value)
    out = {value};
else
    out = cellstr(string(value(:)));
end
end


function ranges = normalize_range_list(value, n_types)
ranges = cell(n_types, 1);
ranges(:) = {[]};

if iscell(value)
    n = min(numel(value), n_types);
    for k = 1:n
        ranges{k} = double(value{k});
    end
elseif isnumeric(value)
    if size(value, 1) == n_types
        for k = 1:n_types
            ranges{k} = double(value(k, :));
        end
    elseif size(value, 2) == n_types
        for k = 1:n_types
            ranges{k} = double(value(:, k)).';
        end
    else
        ranges{1} = double(value);
    end
else
    n = min(numel(value), n_types);
    for k = 1:n
        ranges{k} = char(string(value(k)));
    end
end
end


function g = build_gaussian_kernel(sigma_sec, bin_sec)
if isempty(sigma_sec) || sigma_sec <= 0
    g = 1;
    return;
end

sig_bins = sigma_sec / bin_sec;
hw = max(1, ceil(4 * sig_bins));
x = (-hw:hw);
g = exp(-(x .^ 2) / (2 * sig_bins ^ 2));
g = g / sum(g);
end


function types = build_type_summaries(event_labels, event_ranges, edges, t_centers, ...
    counts_by_type, density_by_type, smoothed_density_by_type, total_event_count)
n_types = numel(event_labels);
types = struct([]);

for k = 1:n_types
    types(k).name = event_labels{k};
    types(k).type_index = k;
    types(k).event_range = event_ranges{k};
    types(k).t_edges = edges(:);
    types(k).t_centers = t_centers(:);
    types(k).counts = counts_by_type(:, k);
    types(k).density = density_by_type(:, k);
    types(k).smoothed_density = smoothed_density_by_type(:, k);
    types(k).total_event_count = total_event_count(k);
end
end


function tag = build_save_tag(file_stem, group_name, detector_name, bin_sec)
sec_tag = strrep(sprintf('%gs', bin_sec), '.', 'p');
if isempty(group_name)
    tag = sprintf('%s_%s_legacy_event_density_%s', file_stem, detector_name, sec_tag);
else
    tag = sprintf('%s_%s_%s_legacy_event_density_%s', file_stem, group_name, detector_name, sec_tag);
end
end
