function result = export_p8_p10_roi_profile_consistency_sources(varargin)
%EXPORT_P8_P10_ROI_PROFILE_CONSISTENCY_SOURCES Export current P8/P10 ROI profiles.
%
% The exporter intentionally writes numeric ROI vectors only.  It does not
% create activation maps or ROI bar figures, so it can be used as a light
% P11 consistency source after P8/P10 numeric xcorr runs finish.

opts = local_parse_inputs(varargin{:});

repo_root = fileparts(fileparts(mfilename('fullpath')));
addpath(genpath(repo_root));
close all force;
set(groot, 'defaultFigureVisible', 'off');

if exist(opts.output_dir, 'dir') ~= 7
    mkdir(opts.output_dir);
end

if opts.group_by_strict_subprocess
    p8_hits_file = opts.strict_band_hits_file;
    p10_hits_file = opts.strict_band_hits_file;
    p8_out = fullfile(opts.output_dir, 'p8_roi_profiles_by_subprocess_long.csv');
    p10_out = fullfile(opts.output_dir, 'p10_roi_profiles_by_subprocess_long.csv');
else
    p8_hits_file = fullfile(opts.repo_root, 'results', ...
        'pipeline8_cross_session_consistency_current', ...
        'p8_top_xcorr_hits_readable.csv');
    p10_hits_file = fullfile(opts.repo_root, 'results', ...
        'pipeline10_cross_session_consistency_current', ...
        'p10_top_xcorr_hits_readable.csv');
    p8_out = fullfile(opts.output_dir, 'p8_roi_profiles_long.csv');
    p10_out = fullfile(opts.output_dir, 'p10_roi_profiles_long.csv');
end

fprintf('[ROI consistency] Reading current P7 candidates...\n');
p8_candidates = local_current_p7_candidates(opts);
fprintf('[ROI consistency] Current P7 candidates: %d\n', numel(p8_candidates));

fprintf('[ROI consistency] Exporting P8 ROI profiles...\n');
if ~opts.export_p8
    fprintf('[ROI consistency] Skipping P8 export by request.\n');
    p8_count = 0;
elseif opts.group_by_strict_subprocess
    p8_count = local_export_p8_profiles_by_subprocess(p8_hits_file, p8_out, p8_candidates, opts);
else
    p8_count = local_export_p8_profiles(p8_hits_file, p8_out, p8_candidates, opts);
end

fprintf('[ROI consistency] Exporting P10 ROI profiles...\n');
if ~opts.export_p10
    fprintf('[ROI consistency] Skipping P10 export by request.\n');
    p10_count = 0;
elseif opts.group_by_strict_subprocess
    p10_count = local_export_p10_profiles_by_subprocess(p10_hits_file, p10_out, p8_candidates, opts);
else
    p10_count = local_export_p10_profiles(p10_hits_file, p10_out, p8_candidates, opts);
end

result = struct();
result.output_dir = opts.output_dir;
result.p8_output_file = p8_out;
result.p10_output_file = p10_out;
result.p8_rows = p8_count;
result.p10_rows = p10_count;

fprintf('[ROI consistency] P8 ROI rows: %d\n', p8_count);
fprintf('[ROI consistency] P10 ROI rows: %d\n', p10_count);
fprintf('[ROI consistency] Output dir: %s\n', opts.output_dir);
end


function opts = local_parse_inputs(varargin)
p = inputParser;
p.addParameter('processed_root', 'E:\DataPons_processed', @(x) ischar(x) || isstring(x));
p.addParameter('datapons_root', 'E:\DataPons', @(x) ischar(x) || isstring(x));
p.addParameter('datasets', {'e10gb1', 'e10fV1', 'e10gh1', 'f12m01'});
p.addParameter('run_tags', {'pv_gsvd100', 'pv_gsvd100_ds', 'pv_hp100', 'pv_roi'});
p.addParameter('output_dir', fullfile('results', 'pipeline_roi_profile_consistency_current'), @(x) ischar(x) || isstring(x));
p.addParameter('feature_reduce', 'mean', @(x) ischar(x) || isstring(x));
p.addParameter('roi_value_mode', 'mean_abs', @(x) ischar(x) || isstring(x));
p.addParameter('current_p7_run_names', {});
p.addParameter('strict_band_hits_file', fullfile('results', ...
    'pipeline8_10_strict_band_coupling_current', ...
    'p8_p10_strict_band_hits_long.csv'), @(x) ischar(x) || isstring(x));
p.addParameter('group_by_strict_subprocess', false);
p.addParameter('export_p8', true);
p.addParameter('export_p10', true);
p.parse(varargin{:});

opts = p.Results;
opts.repo_root = fileparts(fileparts(mfilename('fullpath')));
opts.processed_root = char(string(opts.processed_root));
opts.datapons_root = char(string(opts.datapons_root));
opts.output_dir = char(string(opts.output_dir));
opts.feature_reduce = char(string(opts.feature_reduce));
opts.roi_value_mode = char(string(opts.roi_value_mode));
opts.strict_band_hits_file = char(string(opts.strict_band_hits_file));
opts.datasets = cellstr(string(opts.datasets(:)).');
opts.run_tags = cellstr(string(opts.run_tags(:)).');
opts.current_p7_run_names = cellstr(string(opts.current_p7_run_names(:)).');
opts.group_by_strict_subprocess = local_to_logical(opts.group_by_strict_subprocess);
opts.export_p8 = local_to_logical(opts.export_p8);
opts.export_p10 = local_to_logical(opts.export_p10);
if ~isfolder(opts.output_dir)
    opts.output_dir = fullfile(opts.repo_root, opts.output_dir);
end
if exist(opts.strict_band_hits_file, 'file') ~= 2
    opts.strict_band_hits_file = fullfile(opts.repo_root, opts.strict_band_hits_file);
end
end


function candidates = local_current_p7_candidates(opts)
params = build_bold_cross_modal_coupling_params();
params.processed_root = opts.processed_root;
params.datapons_root = opts.datapons_root;
params.dataset_stems = opts.datasets;
params.current_best_p7_only = true;
if ~isempty(opts.current_p7_run_names)
    params.current_p7_run_names = opts.current_p7_run_names;
end
params.load_existing_xcorr = true;
params.headless = true;
params.close_figures = true;
candidates = discover_completed_bold_cross_modal_coupling_runs(params);
if isempty(candidates)
    return;
end
keep = ismember(lower(string({candidates.run_tag})), lower(string(opts.run_tags)));
candidates = candidates(keep);
end


function n_written = local_export_p8_profiles(hits_file, out_file, candidates, opts)
header = local_profile_header();
fid = fopen(out_file, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(header, ','));
n_written = 0;

if exist(hits_file, 'file') ~= 2
    warning('P8 hit table not found: %s', hits_file);
    return;
end

read_opts = detectImportOptions(hits_file, 'FileType', 'text', ...
    'TextType', 'string', 'PreserveVariableNames', true);
read_opts = setvartype(read_opts, read_opts.VariableNames, 'string');
T = readtable(hits_file, read_opts);
if isempty(T)
    return;
end
mask = strcmpi(T.source_level, "per_density_feature") & ...
    strcmpi(T.density_source_kind, "dimred") & ...
    (strcmpi(T.feature_family, "efun") | strcmpi(T.feature_family, "deconv_efun"));
T = T(mask, :);
if strcmpi(string(pipeline_name), "P10") && ...
        ismember('p9_feature', T.Properties.VariableNames) && ...
        ismember('p9_method_k', T.Properties.VariableNames)
    keep_context = ~ismissing(T.p9_feature) & strlength(T.p9_feature) > 0 & ...
        ~ismissing(T.p9_method_k) & strlength(T.p9_method_k) > 0;
    T = T(keep_context, :);
end
T = local_keep_datasets_and_tags(T, opts);
if isempty(T)
    return;
end

dataset_tag_values = T.dataset + "|" + T.run_tag;
dataset_tag_values = dataset_tag_values(~ismissing(dataset_tag_values) & strlength(dataset_tag_values) > 0);
dataset_tags = unique(dataset_tag_values, 'stable');
for i_key = 1:numel(dataset_tags)
    parts = split(dataset_tags(i_key), "|");
    if numel(parts) < 2 || any(ismissing(parts(1:2)))
        continue;
    end
    dataset = char(parts(1));
    run_tag = char(parts(2));
    idx_run = T.dataset == string(dataset) & T.run_tag == string(run_tag);
    Tr = T(idx_run, :);
    cand = local_find_candidate(candidates, dataset, run_tag);
    if isempty(cand)
        warning('No current P7 candidate for %s | %s', dataset, run_tag);
        continue;
    end

    fprintf('[P8 ROI] %s | %s | groups from %d hit rows\n', dataset, run_tag, height(Tr));
    plot_ctx = build_bold_activation_plot_context(cand.bold_post_file, struct( ...
        'datapons_root', opts.datapons_root, ...
        'feature_reduce', opts.feature_reduce));
    needed = unique(str2double(string(Tr.bold_mode_index)));
    needed = needed(isfinite(needed) & needed >= 1);
    cache = local_compute_roi_vectors(plot_ctx, needed, opts);

    group_keys = unique(Tr.feature_family + "|" + Tr.bold_feature + "|" + Tr.density_name, 'stable');
    for i_group = 1:numel(group_keys)
        gparts = split(group_keys(i_group), "|");
        idx_group = Tr.feature_family == gparts(1) & ...
            Tr.bold_feature == gparts(2) & Tr.density_name == gparts(3);
        Tg = Tr(idx_group, :);
        n_written = n_written + local_write_group_profiles( ...
            fid, 'P8', Tg, cache, 'bold_mode_index', "", "", NaN, "", opts);
    end
    clear plot_ctx cache
end
end


function n_written = local_export_p10_profiles(hits_file, out_file, candidates, opts)
header = local_profile_header();
fid = fopen(out_file, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(header, ','));
n_written = 0;

if exist(hits_file, 'file') ~= 2
    warning('P10 hit table not found: %s', hits_file);
    return;
end

read_opts = detectImportOptions(hits_file, 'FileType', 'text', ...
    'TextType', 'string', 'PreserveVariableNames', true);
read_opts = setvartype(read_opts, read_opts.VariableNames, 'string');
T = readtable(hits_file, read_opts);
if isempty(T)
    return;
end
mask = strcmpi(T.source_level, "per_density_feature") & ...
    strcmpi(T.density_source_kind, "dimred");
T = T(mask, :);
T = local_keep_datasets_and_tags(T, opts);
if isempty(T)
    return;
end

context_values = T.dataset + "|" + T.run_tag + "|" + T.p9_feature + "|" + T.p9_method_k;
context_values = context_values(~ismissing(context_values) & strlength(context_values) > 0);
context_keys = unique(context_values, 'stable');
for i_key = 1:numel(context_keys)
    parts = split(context_keys(i_key), "|");
    if numel(parts) < 4 || any(ismissing(parts(1:4)))
        continue;
    end
    dataset = char(parts(1));
    run_tag = char(parts(2));
    p9_feature = char(parts(3));
    p9_method_k = char(parts(4));
    idx_context = T.dataset == string(dataset) & T.run_tag == string(run_tag) & ...
        T.p9_feature == string(p9_feature) & T.p9_method_k == string(p9_method_k);
    Tc = T(idx_context, :);
    p9_file = local_find_p9_result_file(opts.processed_root, dataset, run_tag, p9_feature, p9_method_k);
    if isempty(p9_file)
        warning('No P9 result file for %s | %s | %s | %s', dataset, run_tag, p9_feature, p9_method_k);
        continue;
    end

    cand = local_find_candidate(candidates, dataset, run_tag);
    bold_post_file = '';
    if ~isempty(cand)
        bold_post_file = cand.bold_post_file;
    end

    fprintf('[P10 ROI] %s | %s | %s | groups from %d hit rows\n', ...
        dataset, run_tag, p9_method_k, height(Tc));
    component_ctx = build_bold_dimred_component_plot_context(p9_file, bold_post_file, struct( ...
        'datapons_root', opts.datapons_root, ...
        'feature_reduce', opts.feature_reduce));
    needed = unique(str2double(string(Tc.bold_component_index)));
    needed = needed(isfinite(needed) & needed >= 1);
    cache = local_compute_roi_vectors(component_ctx, needed, opts);

    group_keys = unique(Tc.feature_family + "|" + Tc.bold_feature + "|" + Tc.density_name, 'stable');
    for i_group = 1:numel(group_keys)
        gparts = split(group_keys(i_group), "|");
        idx_group = Tc.feature_family == gparts(1) & ...
            Tc.bold_feature == gparts(2) & Tc.density_name == gparts(3);
        Tg = Tc(idx_group, :);
        n_written = n_written + local_write_group_profiles( ...
            fid, 'P10', Tg, cache, 'bold_component_index', ...
            Tg.p9_feature(1), Tg.p9_method(1), local_to_double(Tg.p9_k(1)), Tg.p9_method_k(1), opts);
    end
    clear component_ctx cache
end
end


function n_written = local_export_p8_profiles_by_subprocess(hits_file, out_file, candidates, opts)
header = local_profile_header(true);
fid = fopen(out_file, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(header, ','));
n_written = 0;

T = local_read_strict_subprocess_hits(hits_file, 'P8', opts);
if isempty(T)
    return;
end

dataset_tags = unique(T.dataset + "|" + T.run_tag, 'stable');
for i_key = 1:numel(dataset_tags)
    parts = split(dataset_tags(i_key), "|");
    dataset = char(parts(1));
    run_tag = char(parts(2));
    idx_run = T.dataset == string(dataset) & T.run_tag == string(run_tag);
    Tr = T(idx_run, :);
    cand = local_find_candidate(candidates, dataset, run_tag);
    if isempty(cand)
        warning('No current P7 candidate for %s | %s', dataset, run_tag);
        continue;
    end

    fprintf('[P8 ROI subprocess] %s | %s | groups from %d strict hit rows\n', ...
        dataset, run_tag, height(Tr));
    plot_ctx = build_bold_activation_plot_context(cand.bold_post_file, struct( ...
        'datapons_root', opts.datapons_root, ...
        'feature_reduce', opts.feature_reduce));
    needed = unique(str2double(string(Tr.bold_mode_index)));
    needed = needed(isfinite(needed) & needed >= 1);
    cache = local_compute_roi_vectors(plot_ctx, needed, opts);

    group_keys = unique(Tr.feature_family + "|" + Tr.bold_feature + "|" + Tr.density_name, 'stable');
    for i_group = 1:numel(group_keys)
        gparts = split(group_keys(i_group), "|");
        idx_group = Tr.feature_family == gparts(1) & ...
            Tr.bold_feature == gparts(2) & Tr.density_name == gparts(3);
        Tg_parent = Tr(idx_group, :);
        n_written = n_written + local_write_subprocess_group_profiles( ...
            fid, 'P8', Tg_parent, cache, 'bold_mode_index', "", "", NaN, "", opts);
    end
    clear plot_ctx cache
end
end


function n_written = local_export_p10_profiles_by_subprocess(hits_file, out_file, candidates, opts)
header = local_profile_header(true);
fid = fopen(out_file, 'w');
cleanup = onCleanup(@() fclose(fid));
fprintf(fid, '%s\n', strjoin(header, ','));
n_written = 0;

T = local_read_strict_subprocess_hits(hits_file, 'P10', opts);
if isempty(T)
    return;
end

context_keys = unique(T.dataset + "|" + T.run_tag + "|" + T.p9_feature + "|" + T.p9_method_k, 'stable');
for i_key = 1:numel(context_keys)
    parts = split(context_keys(i_key), "|");
    dataset = char(parts(1));
    run_tag = char(parts(2));
    p9_feature = char(parts(3));
    p9_method_k = char(parts(4));
    idx_context = T.dataset == string(dataset) & T.run_tag == string(run_tag) & ...
        T.p9_feature == string(p9_feature) & T.p9_method_k == string(p9_method_k);
    Tc = T(idx_context, :);
    p9_file = local_find_p9_result_file(opts.processed_root, dataset, run_tag, p9_feature, p9_method_k);
    if isempty(p9_file)
        warning('No P9 result file for %s | %s | %s | %s', dataset, run_tag, p9_feature, p9_method_k);
        continue;
    end

    cand = local_find_candidate(candidates, dataset, run_tag);
    bold_post_file = '';
    if ~isempty(cand)
        bold_post_file = cand.bold_post_file;
    end

    fprintf('[P10 ROI subprocess] %s | %s | %s | groups from %d strict hit rows\n', ...
        dataset, run_tag, p9_method_k, height(Tc));
    component_ctx = build_bold_dimred_component_plot_context(p9_file, bold_post_file, struct( ...
        'datapons_root', opts.datapons_root, ...
        'feature_reduce', opts.feature_reduce));
    needed = unique(str2double(string(Tc.bold_component_index)));
    needed = needed(isfinite(needed) & needed >= 1);
    cache = local_compute_roi_vectors(component_ctx, needed, opts);

    group_keys = unique(Tc.feature_family + "|" + Tc.bold_feature + "|" + Tc.density_name, 'stable');
    for i_group = 1:numel(group_keys)
        gparts = split(group_keys(i_group), "|");
        idx_group = Tc.feature_family == gparts(1) & ...
            Tc.bold_feature == gparts(2) & Tc.density_name == gparts(3);
        Tg_parent = Tc(idx_group, :);
        n_written = n_written + local_write_subprocess_group_profiles( ...
            fid, 'P10', Tg_parent, cache, 'bold_component_index', ...
            Tg_parent.p9_feature(1), Tg_parent.p9_method(1), ...
            local_to_double(Tg_parent.p9_k(1)), Tg_parent.p9_method_k(1), opts);
    end
    clear component_ctx cache
end
end


function T = local_read_strict_subprocess_hits(hits_file, pipeline_name, opts)
T = table();
if exist(hits_file, 'file') ~= 2
    warning('Strict band hit table not found: %s', hits_file);
    return;
end

read_opts = detectImportOptions(hits_file, 'FileType', 'text', ...
    'TextType', 'string', 'PreserveVariableNames', true);
read_opts = setvartype(read_opts, read_opts.VariableNames, 'string');
T = readtable(hits_file, read_opts);
if isempty(T)
    return;
end
T = local_normalize_strict_hit_table(T);
mask = strcmpi(T.pipeline, string(pipeline_name)) & ...
    strcmpi(T.source_level, "per_feature") & ...
    strcmpi(T.density_source_kind, "dimred") & ...
    (strcmpi(T.feature_family, "efun") | strcmpi(T.feature_family, "deconv_efun"));
T = T(mask, :);
T = local_keep_datasets_and_tags(T, opts);
if isempty(T)
    return;
end
end


function T = local_normalize_strict_hit_table(T)
if ismember('bold_observable', T.Properties.VariableNames)
    T.observable = T.bold_observable;
end
if ismember('bold_feature_family', T.Properties.VariableNames)
    T.feature_family = T.bold_feature_family;
end
if ismember('density_class', T.Properties.VariableNames)
    source_kind = strings(height(T), 1);
    is_dimred = strcmpi(T.density_class, "dimred_efun_density");
    source_kind(is_dimred) = "dimred";
    source_kind(~is_dimred) = T.density_class(~is_dimred);
    T.density_source_kind = source_kind;
end
if ~ismember('density_display', T.Properties.VariableNames)
    T.density_display = T.density_name;
end
if ~ismember('strict_label', T.Properties.VariableNames)
    T.strict_label = repmat("label_missing", height(T), 1);
end

key_fields = {'pipeline', 'dataset', 'run_tag', 'observable', ...
    'p9_feature', 'p9_method', 'p9_method_k', 'source_level', ...
    'feature_family', 'bold_feature', 'density_name', ...
    'density_source_kind', 'density_condition', 'density_method_k', ...
    'density_index', 'strict_label'};
for i_field = 1:numel(key_fields)
    field_name = key_fields{i_field};
    if ismember(field_name, T.Properties.VariableNames)
        values = string(T.(field_name));
        values(ismissing(values)) = "";
        T.(field_name) = values;
    end
end

groups = strings(height(T), 1);
for i = 1:height(T)
    groups(i) = local_strict_label_group(T.strict_label(i));
end
T.strict_label_group = groups;
end


function n_written = local_write_subprocess_group_profiles(fid, pipeline_name, Tg_parent, cache, ...
        index_field, p9_feature, p9_method, p9_k, p9_method_k, opts)
n_written = 0;
if isempty(Tg_parent)
    return;
end

peaks_parent = str2double(string(Tg_parent.peak_abs_corr));
valid_parent = isfinite(peaks_parent);
parent_peak_sum = sum(peaks_parent(valid_parent), 'omitnan');
if parent_peak_sum <= 0
    parent_peak_sum = NaN;
end

label_groups = unique(Tg_parent.strict_label_group, 'stable');
group_peak_sums = nan(numel(label_groups), 1);
for i_group = 1:numel(label_groups)
    idx = Tg_parent.strict_label_group == label_groups(i_group);
    peaks_i = str2double(string(Tg_parent.peak_abs_corr(idx)));
    group_peak_sums(i_group) = sum(peaks_i(isfinite(peaks_i)), 'omitnan');
end
[~, order] = sort(group_peak_sums, 'descend', 'MissingPlacement', 'last');

for rank_i = 1:numel(order)
    i_group = order(rank_i);
    group_name = label_groups(i_group);
    idx = Tg_parent.strict_label_group == group_name;
    Tg = Tg_parent(idx, :);
    if isempty(Tg)
        continue;
    end
    label_values = unique(Tg.strict_label, 'stable');
    density_indices = unique(Tg.density_index, 'stable');
    peak_sum = group_peak_sums(i_group);
    if isfinite(parent_peak_sum) && parent_peak_sum > 0
        peak_fraction = peak_sum / parent_peak_sum;
    else
        peak_fraction = NaN;
    end
    extra = struct();
    extra.include_subprocess = true;
    extra.strict_label_group = group_name;
    extra.strict_label = strjoin(cellstr(string(label_values(:).')), ';');
    extra.subprocess_rank = rank_i;
    extra.subprocess_peak_abs_corr_sum = peak_sum;
    extra.subprocess_peak_abs_corr_fraction = peak_fraction;
    extra.subprocess_n_hit_rows = height(Tg);
    extra.subprocess_density_indices = strjoin(cellstr(string(density_indices(:).')), ';');
    n_written = n_written + local_write_group_profiles( ...
        fid, pipeline_name, Tg, cache, index_field, p9_feature, p9_method, ...
        p9_k, p9_method_k, opts, extra);
end
end


function T = local_keep_datasets_and_tags(T, opts)
dataset_values = string(T.dataset);
run_tag_values = string(T.run_tag);
keep = ~ismissing(dataset_values) & strlength(dataset_values) > 0 & ...
    ~ismissing(run_tag_values) & strlength(run_tag_values) > 0 & ...
    ismember(lower(dataset_values), lower(string(opts.datasets))) & ...
    ismember(lower(run_tag_values), lower(string(opts.run_tags)));
T = T(keep, :);
end


function cand = local_find_candidate(candidates, dataset, run_tag)
cand = [];
if isempty(candidates)
    return;
end
keep = strcmpi(string({candidates.dataset_stem}), string(dataset)) & ...
    strcmpi(string({candidates.run_tag}), string(run_tag));
idx = find(keep, 1, 'first');
if ~isempty(idx)
    cand = candidates(idx);
end
end


function p9_file = local_find_p9_result_file(processed_root, dataset, run_tag, p9_feature, p9_method_k)
p9_file = '';
mat_dir = fullfile(processed_root, dataset, ...
    'pipeline9_bold_eigenfunction_reduction', run_tag, p9_feature, p9_method_k, 'mat');
L = dir(fullfile(mat_dir, '*.mat'));
if isempty(L)
    return;
end
[~, order] = sort([L.datenum], 'descend');
L = L(order);
p9_file = fullfile(L(1).folder, L(1).name);
end


function cache = local_compute_roi_vectors(plot_ctx, raw_indices, opts)
raw_indices = raw_indices(:).';
n_mode = numel(raw_indices);
n_roi = numel(plot_ctx.region_labels);
values = nan(n_mode, n_roi);

for i_mode = 1:n_mode
    raw_idx = raw_indices(i_mode);
    if raw_idx < 1 || raw_idx > size(plot_ctx.kpm_modes, 1)
        continue;
    end
    mode_obs = double(plot_ctx.kpm_modes(raw_idx, :));
    if ~plot_ctx.direct_voxel_mode
        source_values = mode_obs * plot_ctx.coeff_t;
    else
        source_values = mode_obs;
    end
    voxel_values = map_bold_source_values_to_voxels( ...
        source_values(:), plot_ctx, struct('feature_reduce', opts.feature_reduce));
    values(i_mode, :) = local_reduce_voxels_to_regions( ...
        voxel_values(:), plot_ctx.voxel_region_idx(:), n_roi, opts.roi_value_mode);
end

cache = struct();
cache.raw_indices = raw_indices;
cache.roi_labels = cellstr(string(plot_ctx.region_labels(:)));
cache.values = values;
cache.index_map = containers.Map('KeyType', 'char', 'ValueType', 'double');
for i_mode = 1:n_mode
    cache.index_map(sprintf('%d', raw_indices(i_mode))) = i_mode;
end
end


function roi_values = local_reduce_voxels_to_regions(voxel_values, voxel_region_idx, n_roi, roi_value_mode)
roi_values = nan(1, n_roi);
for i_roi = 1:n_roi
    vals = voxel_values(voxel_region_idx == i_roi);
    vals = vals(isfinite(vals));
    if isempty(vals)
        continue;
    end
    switch lower(char(string(roi_value_mode)))
        case 'mean_abs'
            roi_values(i_roi) = mean(abs(vals), 'omitnan');
        case 'abs_mean'
            roi_values(i_roi) = abs(mean(vals, 'omitnan'));
        case 'real_mean'
            roi_values(i_roi) = mean(real(vals), 'omitnan');
        case 'imag_mean'
            roi_values(i_roi) = mean(imag(vals), 'omitnan');
        otherwise
            error('Unsupported roi_value_mode: %s', roi_value_mode);
    end
end
end


function n_written = local_write_group_profiles(fid, pipeline_name, Tg, cache, index_field, ...
        p9_feature, p9_method, p9_k, p9_method_k, opts, varargin)
n_written = 0;
if isempty(Tg)
    return;
end
extra = struct();
if ~isempty(varargin)
    extra = varargin{1};
end
include_subprocess = isfield(extra, 'include_subprocess') && extra.include_subprocess;

indices = str2double(string(Tg.(index_field)));
peaks = str2double(string(Tg.peak_abs_corr));
valid = isfinite(indices) & indices >= 1 & isfinite(peaks);
indices = indices(valid);
peaks = peaks(valid);
if isempty(indices)
    return;
end

row_ids = nan(size(indices));
for i = 1:numel(indices)
    key = sprintf('%d', indices(i));
    if isKey(cache.index_map, key)
        row_ids(i) = cache.index_map(key);
    end
end
valid = isfinite(row_ids);
indices = indices(valid);
peaks = peaks(valid);
row_ids = row_ids(valid);
if isempty(row_ids)
    return;
end

V = cache.values(row_ids, :);
roi_mean = mean(V, 1, 'omitnan');
if sum(peaks) > 0
    W = peaks(:) ./ sum(peaks);
    roi_weighted = W.' * V;
else
    roi_weighted = roi_mean;
end

selected_text = strjoin(cellstr(string(indices(:).')), ';');
mean_peak = mean(peaks, 'omitnan');
max_peak = max(peaks, [], 'omitnan');
n_hits = height(Tg);
n_selected = numel(indices);

for i_roi = 1:numel(cache.roi_labels)
    row = cell(1, 28 + 7 * double(include_subprocess));
    row{1} = pipeline_name;
    row{2} = Tg.dataset(1);
    row{3} = Tg.run_tag(1);
    row{4} = Tg.observable(1);
    row{5} = p9_feature;
    row{6} = p9_method;
    row{7} = local_num_to_text(p9_k);
    row{8} = p9_method_k;
    row{9} = Tg.feature_family(1);
    row{10} = Tg.bold_feature(1);
    row{11} = Tg.density_name(1);
    row{12} = local_get_table_value(Tg, 'density_display');
    row{13} = local_get_table_value(Tg, 'density_source_kind');
    row{14} = local_get_table_value(Tg, 'density_condition');
    row{15} = local_get_table_value(Tg, 'density_method');
    row{16} = local_get_table_value(Tg, 'density_k');
    row{17} = local_get_table_value(Tg, 'density_method_k');
    row{18} = local_get_table_value(Tg, 'density_threshold');
    row{19} = i_roi;
    row{20} = cache.roi_labels{i_roi};
    row{21} = roi_mean(i_roi);
    row{22} = roi_weighted(i_roi);
    row{23} = n_selected;
    row{24} = selected_text;
    row{25} = mean_peak;
    row{26} = max_peak;
    row{27} = n_hits;
    row{28} = opts.roi_value_mode;
    if include_subprocess
        row{29} = extra.strict_label_group;
        row{30} = extra.strict_label;
        row{31} = extra.subprocess_rank;
        row{32} = extra.subprocess_peak_abs_corr_sum;
        row{33} = extra.subprocess_peak_abs_corr_fraction;
        row{34} = extra.subprocess_n_hit_rows;
        row{35} = extra.subprocess_density_indices;
    end
    local_write_csv_row(fid, row);
    n_written = n_written + 1;
end
end


function header = local_profile_header(varargin)
include_subprocess = false;
if ~isempty(varargin)
    include_subprocess = varargin{1};
end
header = { ...
    'pipeline', 'dataset', 'run_tag', 'observable', ...
    'p9_feature', 'p9_method', 'p9_k', 'p9_method_k', ...
    'feature_family', 'bold_feature', ...
    'density_name', 'density_display', 'density_source_kind', ...
    'density_condition', 'density_method', 'density_k', 'density_method_k', ...
    'density_threshold', 'roi_index', 'roi_label', ...
    'roi_value_mean', 'roi_value_weighted', 'n_selected', ...
    'selected_indices', 'mean_peak_abs_corr', 'max_peak_abs_corr', ...
    'n_hit_rows', 'roi_value_mode'};
if include_subprocess
    header = [header, { ...
        'strict_label_group', 'strict_label', 'subprocess_rank', ...
        'subprocess_peak_abs_corr_sum', 'subprocess_peak_abs_corr_fraction', ...
        'subprocess_n_hit_rows', 'subprocess_density_indices'}];
end
end


function value = local_get_table_value(T, field_name)
if ismember(field_name, T.Properties.VariableNames)
    value = T.(field_name)(1);
else
    value = "";
end
end


function text = local_num_to_text(value)
if isstring(value) || ischar(value)
    text = char(string(value));
elseif isempty(value) || ~isfinite(value)
    text = '';
else
    text = sprintf('%g', value);
end
end


function value = local_to_double(value)
if isnumeric(value)
    value = double(value);
elseif isstring(value) || ischar(value)
    value = str2double(value);
else
    value = NaN;
end
if isempty(value)
    value = NaN;
else
    value = value(1);
end
end


function value = local_to_logical(value)
if islogical(value)
    value = value(1);
elseif isnumeric(value)
    value = value(1) ~= 0;
else
    text = lower(strtrim(char(string(value))));
    value = any(strcmp(text, {'true', '1', 'yes', 'y', 'on'}));
end
end


function group = local_strict_label_group(label)
label = string(label);
if ismissing(label) || strlength(label) == 0
    group = "unlabelled";
    return;
end
label = lower(strtrim(char(label)));
switch label
    case 'theta_selective'
        group = "theta";
    case {'ripple_gamma_no_theta', 'gamma_selective', 'ripple_selective'}
        group = "ripple_gamma";
    case 'mixed_or_partial'
        group = "mixed";
    case 'inactive'
        group = "inactive";
    otherwise
        group = "unlabelled";
end
end


function local_write_csv_row(fid, row)
for i = 1:numel(row)
    if i > 1
        fprintf(fid, ',');
    end
    fprintf(fid, '%s', local_csv_escape(row{i}));
end
fprintf(fid, '\n');
end


function text = local_csv_escape(value)
if isnumeric(value)
    if isempty(value) || ~isfinite(value)
        text = '';
    else
        text = sprintf('%.10g', value);
    end
    return;
end
text = char(string(value));
text = strrep(text, '"', '""');
if contains(text, ',') || contains(text, '"') || contains(text, newline) || contains(text, sprintf('\r'))
    text = ['"', text, '"'];
end
end
