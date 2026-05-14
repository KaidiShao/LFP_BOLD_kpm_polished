function info = build_bold_base_observable_info(D, n_var, source_tag, label_suffix)
%BUILD_BOLD_BASE_OBSERVABLE_INFO Build the shared base observable-info table.

if nargin < 2 || isempty(n_var)
    n_var = size(D.data, 2);
end
if nargin < 3 || isempty(source_tag)
    source_tag = 'bold';
end
if nargin < 4 || isempty(label_suffix)
    label_suffix = '';
end

observable_idx = (1:n_var).';
source = repmat({char(string(source_tag))}, n_var, 1);

if isfield(D, 'variable_labels') && numel(D.variable_labels) == n_var
    observable_label = cellstr(string(D.variable_labels(:)));
else
    observable_label = cellstr(compose("bold_var%04d", observable_idx));
end

if ~isempty(label_suffix)
    observable_label = cellstr(strcat(string(observable_label), string(label_suffix)));
end

info = table(observable_idx, source, observable_label, ...
    'VariableNames', {'observable_idx', 'source', 'observable_label'});

if isfield(D, 'variable_info') && height(D.variable_info) == n_var
    extra = D.variable_info;
    duplicate_names = intersect(info.Properties.VariableNames, extra.Properties.VariableNames);
    if ~isempty(duplicate_names)
        extra(:, duplicate_names) = [];
    end
    info = [info extra];
end
end
