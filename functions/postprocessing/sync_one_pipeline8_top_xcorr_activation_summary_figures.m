function summary_result = sync_one_pipeline8_top_xcorr_activation_summary_figures(run_info, source, output_root)
%SYNC_ONE_PIPELINE8_TOP_XCORR_ACTIVATION_SUMMARY_FIGURES Copy activation-map PNGs into shared summaries.

if nargin < 1 || ~isstruct(run_info)
    error('run_info struct is required.');
end
if nargin < 2
    source = struct();
end
if nargin < 3 || isempty(output_root)
    output_root = io_project.get_project_processed_root();
end

dataset_stem = char(string(local_get_field(run_info, 'dataset_stem', '')));
run_name = char(string(local_get_field(run_info, 'run_name', '')));
observable_mode = char(string(local_get_field(run_info, 'observable_mode', '')));
residual_form = char(string(local_get_field(run_info, 'residual_form', '')));
run_tag = char(string(local_get_field(run_info, 'run_tag', '')));
if isempty(run_tag)
    run_tag = make_bold_cross_modal_run_tag(run_info);
end
activation_root = char(string(local_resolve_activation_root(run_info, source)));

summary_result = struct();
summary_result.status = 'no_source_dir';
summary_result.scope = local_classify_observable_scope(observable_mode);
summary_result.dataset_stem = dataset_stem;
summary_result.run_name = run_name;
summary_result.run_tag = run_tag;
summary_result.observable_mode = observable_mode;
summary_result.residual_form = residual_form;
summary_result.summary_dir = '';
summary_result.copied_files = {};
summary_result.n_copied_files = 0;
summary_result.activation_root = activation_root;

if isempty(dataset_stem) || isempty(run_name) || isempty(observable_mode)
    summary_result.status = 'missing_run_info';
    return;
end
if isempty(activation_root) || exist(activation_root, 'dir') ~= 7
    return;
end

summary_root = io_project.get_pipeline_summary_dir(output_root, 8, 'figures_top_xcorr_activation_map');
summary_dir = fullfile(summary_root, sprintf('%s_observables', summary_result.scope));
if exist(summary_dir, 'dir') ~= 7
    mkdir(summary_dir);
end

group_dirs = local_collect_activation_group_dirs(activation_root);
if isempty(group_dirs)
    summary_result.status = 'no_activation_groups';
    summary_result.summary_dir = summary_dir;
    return;
end

copied_files = {};
for i_group = 1:numel(group_dirs)
    src_group = group_dirs{i_group};
    [~, group_name] = fileparts(src_group);
    L = dir(fullfile(src_group, '*.png'));
    for i_file = 1:numel(L)
        src_file = fullfile(L(i_file).folder, L(i_file).name);
        dst_name = local_build_flat_filename(dataset_stem, run_tag, group_name, L(i_file).name);
        dst_file = fullfile(summary_dir, dst_name);
        copyfile(src_file, dst_file, 'f');
        copied_files{end + 1, 1} = dst_file; %#ok<AGROW>
    end
end

if isempty(copied_files)
    summary_result.status = 'no_pngs';
else
    summary_result.status = 'ok';
end
summary_result.summary_dir = summary_dir;
summary_result.copied_files = copied_files;
summary_result.n_copied_files = numel(copied_files);
end


function activation_root = local_resolve_activation_root(run_info, source)
activation_root = local_get_field(source, 'activation_dir', '');
if isempty(activation_root)
    activation_root = local_get_field(source, 'primary_activation_dir', '');
end
if isempty(activation_root)
    activation_root = local_get_field(run_info, 'activation_root', '');
end
activation_root = char(string(activation_root));
end


function scope = local_classify_observable_scope(observable_mode)
global_modes = {'svd', 'global_svd100', 'gsvd100_ds', 'global_slow_band_power_svd100', ...
    'roi_mean', 'roi_mean_slow_band_power'};
if any(strcmpi(observable_mode, global_modes))
    scope = 'global';
else
    scope = 'local';
end
end


function group_dirs = local_collect_activation_group_dirs(activation_root)
patterns = { ...
    fullfile(activation_root, 'activation_maps_top*'), ...
    fullfile(activation_root, 'act'), ...
    fullfile(activation_root, 'act', 'efun'), ...
    fullfile(activation_root, 'act', 'deconv_efun')};
group_dirs = {};
for i_pat = 1:numel(patterns)
    L = dir(patterns{i_pat});
    L = L([L.isdir]);
    for i = 1:numel(L)
        group_dirs{end + 1, 1} = fullfile(L(i).folder, L(i).name); %#ok<AGROW>
    end
end
group_dirs = unique(group_dirs, 'stable');
end


function dst_name = local_build_flat_filename(dataset_stem, run_name, group_name, original_name)
[~, base, ext] = fileparts(original_name);
base = strrep(base, '_activation_reference', '');
group_slug = local_group_slug(group_name);
dst_name = sprintf('%s__%s__%s__%s%s', ...
    dataset_stem, run_name, group_slug, base, ext);
end


function slug = local_group_slug(group_name)
name = lower(char(string(group_name)));
switch name
    case 'activation_maps_top5_blp_event_density'
        slug = 'evt';
    case 'activation_maps_top5_blp_raw_eigenfunction_density'
        slug = 'raw';
    case 'activation_maps_top5_combined'
        slug = 'cmb';
    otherwise
        slug = regexprep(name, '^activation_maps_top\\d+_', '');
        slug = regexprep(slug, '[^a-z0-9]+', '_');
        slug = regexprep(slug, '^_+|_+$', '');
        if isempty(slug)
            slug = 'grp';
        end
    end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
