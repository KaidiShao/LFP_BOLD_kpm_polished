function [component_indices, top_row_indices] = resolve_bold_dimred_top_xcorr_component_indices(top_table, n_components, top_n)
%RESOLVE_BOLD_DIMRED_TOP_XCORR_COMPONENT_INDICES Resolve unique component indices.

if nargin < 3 || isempty(top_n)
    top_n = inf;
end
component_indices = nan(1, 0);
top_row_indices = nan(1, 0);
if isempty(top_table) || ~istable(top_table)
    return;
end
if ismember('bold_component_index', top_table.Properties.VariableNames)
    values = top_table.bold_component_index;
elseif ismember('bold_mode_index', top_table.Properties.VariableNames)
    values = top_table.bold_mode_index;
else
    return;
end

for i = 1:height(top_table)
    if numel(component_indices) >= top_n
        break;
    end
    idx = round(double(values(i)));
    if ~isfinite(idx) || idx < 1 || idx > n_components
        continue;
    end
    if any(component_indices == idx)
        continue;
    end
    component_indices(end + 1) = idx; %#ok<AGROW>
    top_row_indices(end + 1) = i; %#ok<AGROW>
end
component_indices = component_indices(:).';
top_row_indices = top_row_indices(:).';
end
