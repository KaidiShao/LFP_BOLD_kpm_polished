function group_specs = build_bold_top_xcorr_group_specs(xcorr_out, params)
%BUILD_BOLD_TOP_XCORR_GROUP_SPECS Build pipeline 8 plotting groups from xcorr tables.

if nargin < 1 || isempty(xcorr_out)
    error('xcorr_out is required.');
end
if nargin < 2 || isempty(params)
    params = struct();
end

export_combined = local_get_field(params, 'export_combined', true);
export_by_density = local_get_field(params, 'export_by_density', true);
export_by_density_feature = local_get_field(params, 'export_by_density_feature', false);
top_n = local_get_field(params, 'top_n', 5);

group_specs = repmat(local_empty_group_spec(), 0, 1);
if export_combined && isfield(xcorr_out, 'top_table') && ~isempty(xcorr_out.top_table)
    group_spec = local_empty_group_spec();
    group_spec.scope = 'combined';
    group_spec.display_name = 'combined';
    group_spec.slug = 'combined';
    group_spec.density_name = '';
    group_spec.top_table = xcorr_out.top_table;
    group_specs(end + 1, 1) = group_spec; %#ok<AGROW>
end

if export_by_density
    [density_names, field_names, top_tables] = local_density_group_top_tables(xcorr_out, top_n);
    for i_group = 1:numel(density_names)
        top_table_i = top_tables{i_group};
        if isempty(top_table_i)
            continue;
        end
        group_spec = local_empty_group_spec();
        group_spec.scope = 'density';
        group_spec.display_name = density_names{i_group};
        group_spec.slug = field_names{i_group};
        group_spec.density_name = density_names{i_group};
        group_spec.top_table = top_table_i;
        group_specs(end + 1, 1) = group_spec; %#ok<AGROW>
    end
end

if ~export_by_density_feature
    return;
end

[group_names, field_names, density_names, feature_names, top_tables] = ...
    local_density_feature_group_top_tables(xcorr_out, top_n);
for i_group = 1:numel(group_names)
    top_table_i = top_tables{i_group};
    if isempty(top_table_i)
        continue;
    end
    group_spec = local_empty_group_spec();
    group_spec.scope = 'density_feature';
    group_spec.display_name = group_names{i_group};
    group_spec.slug = field_names{i_group};
    group_spec.density_name = density_names{i_group};
    group_spec.feature_name = feature_names{i_group};
    group_spec.top_table = top_table_i;
    group_specs(end + 1, 1) = group_spec; %#ok<AGROW>
end
end


function [density_names, field_names, top_tables] = local_density_group_top_tables(xcorr_out, top_n)
density_names = {};
field_names = {};
top_tables = {};
if isfield(xcorr_out, 'top_table_by_density') && isstruct(xcorr_out.top_table_by_density) && ...
        ~isempty(fieldnames(xcorr_out.top_table_by_density))
    top_struct = xcorr_out.top_table_by_density;
    if isfield(xcorr_out, 'density_group_names') && isfield(xcorr_out, 'density_group_fields') && ...
            numel(xcorr_out.density_group_names) == numel(xcorr_out.density_group_fields)
        density_names = cellstr(string(xcorr_out.density_group_names(:)).');
        field_names = cellstr(string(xcorr_out.density_group_fields(:)).');
    else
        field_names = fieldnames(top_struct).';
        density_names = field_names;
    end
    top_tables = cell(size(field_names));
    for i_group = 1:numel(field_names)
        top_tables{i_group} = top_struct.(field_names{i_group});
    end
    return;
end

if ~isfield(xcorr_out, 'peak_table') || isempty(xcorr_out.peak_table)
    return;
end

[top_struct, density_names, field_names] = local_build_density_group_top_struct( ...
    xcorr_out.peak_table, top_n);
top_tables = cell(size(field_names));
for i_group = 1:numel(field_names)
    top_tables{i_group} = top_struct.(field_names{i_group});
end
end


function [top_struct, density_names, field_names] = ...
        local_build_density_group_top_struct(peak_table, top_n)
top_struct = struct();
density_names = {};
field_names = {};
if isempty(peak_table) || ~ismember('density_name', peak_table.Properties.VariableNames)
    return;
end
all_names = cellstr(string(peak_table.density_name));
[density_names, ~] = unique(all_names, 'stable');
density_names = density_names(:).';
field_names = cell(size(density_names));
for i_group = 1:numel(density_names)
    density_name = density_names{i_group};
    field_name = local_density_group_field_name(density_name);
    mask = strcmp(all_names, density_name);
    group_table = local_sort_peak_table(peak_table(mask, :));
    if isempty(group_table)
        top_struct.(field_name) = group_table;
    else
        top_struct.(field_name) = group_table(1:min(height(group_table), top_n), :);
    end
    field_names{i_group} = field_name;
end
end


function [group_names, field_names, density_names, feature_names, top_tables] = ...
        local_density_feature_group_top_tables(xcorr_out, top_n)
group_names = {};
field_names = {};
density_names = {};
feature_names = {};
top_tables = {};
if isfield(xcorr_out, 'top_table_by_density_feature') && ...
        isstruct(xcorr_out.top_table_by_density_feature) && ...
        ~isempty(fieldnames(xcorr_out.top_table_by_density_feature))
    top_struct = xcorr_out.top_table_by_density_feature;
    if isfield(xcorr_out, 'density_feature_group_names') && ...
            isfield(xcorr_out, 'density_feature_group_fields') && ...
            isfield(xcorr_out, 'density_feature_density_names') && ...
            isfield(xcorr_out, 'density_feature_feature_names') && ...
            numel(xcorr_out.density_feature_group_fields) == ...
            numel(xcorr_out.density_feature_group_names)
        group_names = cellstr(string(xcorr_out.density_feature_group_names(:)).');
        field_names = cellstr(string(xcorr_out.density_feature_group_fields(:)).');
        density_names = cellstr(string(xcorr_out.density_feature_density_names(:)).');
        feature_names = cellstr(string(xcorr_out.density_feature_feature_names(:)).');
    else
        field_names = fieldnames(top_struct).';
        group_names = field_names;
        density_names = field_names;
        feature_names = repmat({''}, size(field_names));
    end
    top_tables = cell(size(field_names));
    for i_group = 1:numel(field_names)
        top_tables{i_group} = top_struct.(field_names{i_group});
    end
    return;
end

if ~isfield(xcorr_out, 'peak_table') || isempty(xcorr_out.peak_table)
    return;
end

[top_struct, group_names, field_names, density_names, feature_names] = ...
    local_build_density_feature_group_top_struct(xcorr_out.peak_table, top_n);
top_tables = cell(size(field_names));
for i_group = 1:numel(field_names)
    top_tables{i_group} = top_struct.(field_names{i_group});
end
end


function [top_struct, group_names, field_names, density_names, feature_names] = ...
        local_build_density_feature_group_top_struct(peak_table, top_n)
top_struct = struct();
group_names = {};
field_names = {};
density_names = {};
feature_names = {};
if isempty(peak_table) || ~all(ismember({'density_name', 'bold_feature'}, ...
        peak_table.Properties.VariableNames))
    return;
end
all_density = cellstr(string(peak_table.density_name));
all_feature = cellstr(string(peak_table.bold_feature));
all_keys = strcat(all_density, {' | '}, all_feature);
[keys, ia] = unique(all_keys, 'stable');
group_names = keys(:).';
density_names = all_density(ia).';
feature_names = all_feature(ia).';
field_names = cell(size(group_names));
for i_group = 1:numel(group_names)
    density_name = density_names{i_group};
    feature_name = feature_names{i_group};
    field_name = local_density_feature_group_field_name(density_name, feature_name);
    mask = strcmp(all_density, density_name) & strcmp(all_feature, feature_name);
    group_table = local_sort_peak_table(peak_table(mask, :));
    if isempty(group_table)
        top_struct.(field_name) = group_table;
    else
        top_struct.(field_name) = group_table(1:min(height(group_table), top_n), :);
    end
    field_names{i_group} = field_name;
end
end


function T = local_sort_peak_table(T)
if isempty(T)
    return;
end
finite_peak = isfinite(T.peak_abs_corr);
T = [ ...
    sortrows(T(finite_peak, :), 'peak_abs_corr', 'descend'); ...
    T(~finite_peak, :)];
end


function field_name = local_density_group_field_name(name)
field_name = char(matlab.lang.makeValidName(char(string(name))));
if isempty(field_name)
    field_name = 'density_group';
end
end


function field_name = local_density_feature_group_field_name(density_name, feature_name)
field_name = char(matlab.lang.makeValidName(sprintf('%s__%s', ...
    char(string(density_name)), char(string(feature_name)))));
if isempty(field_name)
    field_name = 'density_feature_group';
end
end


function group = local_empty_group_spec()
group = struct( ...
    'scope', '', ...
    'display_name', '', ...
    'slug', '', ...
    'density_name', '', ...
    'feature_name', '', ...
    'top_table', table());
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
