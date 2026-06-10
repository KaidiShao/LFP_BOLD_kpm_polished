function out = publish_bold_efun_density_cross_correlation_output(out, params)
%PUBLISH_BOLD_EFUN_DENSITY_CROSS_CORRELATION_OUTPUT Save and plot pipeline 8 xcorr outputs.

if params.save_results
    if exist(params.save_dir, 'dir') ~= 7
        mkdir(params.save_dir);
    end
    save_file = fullfile(params.save_dir, [params.save_tag, '.mat']);
    csv_file = fullfile(params.save_dir, [params.save_tag, '_peaks.csv']);
    top_csv_file = fullfile(params.save_dir, [params.save_tag, '_top.csv']);
    if local_get_field(params, 'save_mat', true)
        out_to_save = local_prepare_saved_output(out, params);
        save_mat_variable_atomic(save_file, 'out', out_to_save);
    else
        save_file = '';
    end

    if params.export_combined && ~isempty(out.peak_table)
        writetable(out.peak_table, csv_file);
    end
    if params.export_combined && ~isempty(out.top_table)
        writetable(out.top_table, top_csv_file);
    end

    density_group_files = repmat(struct( ...
        'density_name', '', ...
        'field_name', '', ...
        'peak_csv', '', ...
        'top_csv', ''), numel(out.density_group_names), 1);
    by_density_dir = fullfile(params.save_dir, 'density');
    if params.export_by_density && ~isempty(out.density_group_names) && exist(by_density_dir, 'dir') ~= 7
        mkdir(by_density_dir);
    end
    if params.export_by_density
        for i_group = 1:numel(out.density_group_names)
            field_name = out.density_group_fields{i_group};
            peak_csv_i = fullfile(by_density_dir, sprintf('%s_peaks__%s.csv', ...
                params.save_tag, field_name));
            top_csv_i = fullfile(by_density_dir, sprintf('%s_top__%s.csv', ...
                params.save_tag, field_name));
            peak_table_i = out.peak_table_by_density.(field_name);
            top_table_i = out.top_table_by_density.(field_name);
            if ~isempty(peak_table_i)
                writetable(peak_table_i, peak_csv_i);
            end
            if ~isempty(top_table_i)
                writetable(top_table_i, top_csv_i);
            end
            density_group_files(i_group).density_name = out.density_group_names{i_group};
            density_group_files(i_group).field_name = field_name;
            density_group_files(i_group).peak_csv = peak_csv_i;
            density_group_files(i_group).top_csv = top_csv_i;
        end
    end

    density_feature_group_files = local_export_density_feature_tables(out, params);

    out.save_paths = struct( ...
        'mat_file', save_file, ...
        'peak_csv', csv_file, ...
        'top_csv', top_csv_file, ...
        'by_density_dir', by_density_dir, ...
        'density_group_files', density_group_files, ...
        'by_density_feature_dir', fullfile(params.save_dir, 'feature'), ...
        'density_feature_group_files', density_feature_group_files);
end

if params.make_figures && params.export_combined
    out.figure_paths = plot_bold_efun_density_cross_correlation_summary(out, params.plot);
end
if params.make_figures && params.export_by_density
    out.figure_paths_by_density = local_plot_density_group_figures( ...
        out, out.density_group_names, out.density_group_fields, params);
end
if params.make_figures && isfield(params, 'export_by_density_feature') && ...
        params.export_by_density_feature
    out.figure_paths_by_density_feature = local_plot_density_feature_group_figures(out, params);
end

if local_should_drop_source_results(params)
    out = local_prepare_saved_output(out, params);
end
end


function out_to_save = local_prepare_saved_output(out, params)
out_to_save = out;
if local_should_drop_source_results(params)
    out_to_save.source_results = {};
    out_to_save.storage_policy = 'tables_only_no_source_results';
elseif isfield(params, 'save_full_source_results') && ~params.save_full_source_results
    out_to_save = local_strip_heavy_source_result_fields(out_to_save);
end
end


function tf = local_should_drop_source_results(params)
tf = isfield(params, 'save_source_results') && ~params.save_source_results;
end


function out = local_strip_heavy_source_result_fields(out)
if ~isfield(out, 'source_results') || isempty(out.source_results)
    return;
end

for i_src = 1:numel(out.source_results)
    if iscell(out.source_results)
        sr = out.source_results{i_src};
    else
        sr = out.source_results(i_src);
    end

    if isfield(sr, 'best_candidate') && isstruct(sr.best_candidate)
        sr.best_candidate = local_rmfield_if_present(sr.best_candidate, {'X'});
    end
    sr = local_rmfield_if_present(sr, {'best_bold_matrix'});
    if isfield(sr, 'best_maps') && isstruct(sr.best_maps)
        sr.best_maps = local_rmfield_if_present(sr.best_maps, {'corr_cube', 'valid_count'});
    end

    if iscell(out.source_results)
        out.source_results{i_src} = sr;
    else
        out.source_results(i_src) = sr;
    end
end
out.storage_policy = 'slim_source_results';
end


function S = local_rmfield_if_present(S, names)
for i = 1:numel(names)
    name = names{i};
    if isfield(S, name)
        S = rmfield(S, name);
    end
end
end


function group_files = local_export_density_feature_tables(out, params)
template = struct( ...
    'group_name', '', ...
    'field_name', '', ...
    'density_name', '', ...
    'feature_name', '', ...
    'peak_csv', '', ...
    'top_csv', '');
group_files = repmat(template, 0, 1);
if ~isfield(out, 'density_feature_group_fields') || ...
        isempty(out.density_feature_group_fields)
    return;
end

by_group_dir = fullfile(params.save_dir, 'feature');
if exist(by_group_dir, 'dir') ~= 7
    mkdir(by_group_dir);
end

n_groups = numel(out.density_feature_group_fields);
group_files = repmat(template, n_groups, 1);
for i_group = 1:n_groups
    field_name = out.density_feature_group_fields{i_group};
    feature_name = out.density_feature_feature_names{i_group};
    group_dir = local_density_feature_group_dir(by_group_dir, feature_name);
    if exist(group_dir, 'dir') ~= 7
        mkdir(group_dir);
    end
    peak_csv_i = fullfile(group_dir, sprintf('%s_peaks__%s.csv', ...
        params.save_tag, field_name));
    top_csv_i = fullfile(group_dir, sprintf('%s_top__%s.csv', ...
        params.save_tag, field_name));
    peak_table_i = out.peak_table_by_density_feature.(field_name);
    top_table_i = out.top_table_by_density_feature.(field_name);
    if ~isempty(peak_table_i)
        writetable(peak_table_i, peak_csv_i);
    end
    if ~isempty(top_table_i)
        writetable(top_table_i, top_csv_i);
    end
    group_files(i_group).group_name = out.density_feature_group_names{i_group};
    group_files(i_group).field_name = field_name;
    group_files(i_group).density_name = out.density_feature_density_names{i_group};
    group_files(i_group).feature_name = feature_name;
    group_files(i_group).peak_csv = peak_csv_i;
    group_files(i_group).top_csv = top_csv_i;
end
end


function figure_paths_by_density = local_plot_density_group_figures( ...
        out, density_group_names, density_group_fields, params)
figure_paths_by_density = repmat(struct( ...
    'density_name', '', ...
    'field_name', '', ...
    'summary_png', '', ...
    'top_curves_png', '', ...
    'top_overlay_png', ''), numel(density_group_names), 1);
if isempty(density_group_names)
    return;
end

plot_dir = fullfile(params.save_dir, 'density');
if exist(plot_dir, 'dir') ~= 7
    mkdir(plot_dir);
end

for i_group = 1:numel(density_group_names)
    density_name = density_group_names{i_group};
    field_name = density_group_fields{i_group};
    out_i = local_make_density_group_out(out, density_name, field_name);
    plot_params = params.plot;
    plot_params.save_dir = plot_dir;
    plot_params.save_tag = sprintf('%s__%s', params.save_tag, field_name);
    paths_i = plot_bold_efun_density_cross_correlation_summary(out_i, plot_params);
    figure_paths_by_density(i_group).density_name = density_name;
    figure_paths_by_density(i_group).field_name = field_name;
    figure_paths_by_density(i_group).summary_png = paths_i.summary_png;
    figure_paths_by_density(i_group).top_curves_png = paths_i.top_curves_png;
    figure_paths_by_density(i_group).top_overlay_png = paths_i.top_overlay_png;
end
end


function figure_paths_by_group = local_plot_density_feature_group_figures(out, params)
template = struct( ...
    'group_name', '', ...
    'field_name', '', ...
    'density_name', '', ...
    'feature_name', '', ...
    'summary_png', '', ...
    'top_curves_png', '', ...
    'top_overlay_png', '');
figure_paths_by_group = repmat(template, 0, 1);
if ~isfield(out, 'density_feature_group_fields') || ...
        isempty(out.density_feature_group_fields)
    return;
end

plot_root = fullfile(params.save_dir, 'feature');
if exist(plot_root, 'dir') ~= 7
    mkdir(plot_root);
end

n_groups = numel(out.density_feature_group_fields);
figure_paths_by_group = repmat(template, n_groups, 1);
for i_group = 1:n_groups
    field_name = out.density_feature_group_fields{i_group};
    density_name = out.density_feature_density_names{i_group};
    feature_name = out.density_feature_feature_names{i_group};
    out_i = local_make_density_feature_group_out( ...
        out, density_name, feature_name, field_name);
    plot_dir = local_density_feature_group_dir(plot_root, feature_name);
    if exist(plot_dir, 'dir') ~= 7
        mkdir(plot_dir);
    end
    plot_params = params.plot;
    plot_params.save_dir = plot_dir;
    plot_params.save_tag = sprintf('%s__%s', params.save_tag, field_name);
    paths_i = plot_bold_efun_density_cross_correlation_summary(out_i, plot_params);
    figure_paths_by_group(i_group).group_name = out.density_feature_group_names{i_group};
    figure_paths_by_group(i_group).field_name = field_name;
    figure_paths_by_group(i_group).density_name = density_name;
    figure_paths_by_group(i_group).feature_name = feature_name;
    figure_paths_by_group(i_group).summary_png = paths_i.summary_png;
    figure_paths_by_group(i_group).top_curves_png = paths_i.top_curves_png;
    figure_paths_by_group(i_group).top_overlay_png = paths_i.top_overlay_png;
end
end


function group_dir = local_density_feature_group_dir(root_dir, feature_name)
group_dir = fullfile(root_dir, local_feature_family(feature_name));
end


function family = local_feature_family(feature_name)
name = lower(char(string(feature_name)));
if contains(name, 'deconv')
    family = 'deconv_efun';
else
    family = 'efun';
end
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end


function out_i = local_make_density_group_out(out, density_name, field_name)
out_i = out;
if isfield(out, 'peak_table_by_density') && isstruct(out.peak_table_by_density) && ...
        isfield(out.peak_table_by_density, field_name)
    out_i.peak_table = out.peak_table_by_density.(field_name);
else
    mask = strcmp(cellstr(string(out.peak_table.density_name)), density_name);
    out_i.peak_table = out.peak_table(mask, :);
end
if isfield(out, 'top_table_by_density') && isstruct(out.top_table_by_density) && ...
        isfield(out.top_table_by_density, field_name)
    out_i.top_table = out.top_table_by_density.(field_name);
else
    out_i.top_table = out_i.peak_table(1:min(height(out_i.peak_table), out.params.top_n), :);
end
out_i.source_results = local_filter_source_results_by_density(out.source_results, density_name);
end


function out_i = local_make_density_feature_group_out( ...
        out, density_name, feature_name, field_name)
out_i = out;
if isfield(out, 'peak_table_by_density_feature') && ...
        isstruct(out.peak_table_by_density_feature) && ...
        isfield(out.peak_table_by_density_feature, field_name)
    out_i.peak_table = out.peak_table_by_density_feature.(field_name);
else
    mask = strcmp(cellstr(string(out.peak_table.density_name)), density_name) & ...
        strcmp(cellstr(string(out.peak_table.bold_feature)), feature_name);
    out_i.peak_table = out.peak_table(mask, :);
end

if isfield(out, 'top_table_by_density_feature') && ...
        isstruct(out.top_table_by_density_feature) && ...
        isfield(out.top_table_by_density_feature, field_name)
    out_i.top_table = out.top_table_by_density_feature.(field_name);
else
    out_i.top_table = out_i.peak_table( ...
        1:min(height(out_i.peak_table), out.params.top_n), :);
end
out_i.source_results = local_filter_source_results_by_density_feature( ...
    out.source_results, density_name, feature_name);
end


function source_results = local_filter_source_results_by_density(source_results_in, density_name)
if isempty(source_results_in)
    source_results = source_results_in;
    return;
end

if iscell(source_results_in)
    keep = false(size(source_results_in));
    for i_item = 1:numel(source_results_in)
        keep(i_item) = strcmp( ...
            string(source_results_in{i_item}.best_candidate.name), ...
            string(density_name));
    end
    source_results = source_results_in(keep);
else
    keep = false(size(source_results_in));
    for i_item = 1:numel(source_results_in)
        keep(i_item) = strcmp( ...
            string(source_results_in(i_item).best_candidate.name), ...
            string(density_name));
    end
    source_results = source_results_in(keep);
end
end


function source_results = local_filter_source_results_by_density_feature( ...
        source_results_in, density_name, feature_name)
if isempty(source_results_in)
    source_results = source_results_in;
    return;
end

if iscell(source_results_in)
    keep = false(size(source_results_in));
    for i_item = 1:numel(source_results_in)
        keep(i_item) = local_source_result_matches( ...
            source_results_in{i_item}, density_name, feature_name);
    end
    source_results = source_results_in(keep);
else
    keep = false(size(source_results_in));
    for i_item = 1:numel(source_results_in)
        keep(i_item) = local_source_result_matches( ...
            source_results_in(i_item), density_name, feature_name);
    end
    source_results = source_results_in(keep);
end
end


function tf = local_source_result_matches(item, density_name, feature_name)
tf = false;
if ~isstruct(item) || ~isfield(item, 'best_candidate') || ...
        ~isfield(item, 'best_feature')
    return;
end
tf = strcmp(string(item.best_candidate.name), string(density_name)) && ...
    strcmp(string(item.best_feature.name), string(feature_name));
end
