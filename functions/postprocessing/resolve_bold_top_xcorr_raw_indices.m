function [raw_indices, top_row_indices] = resolve_bold_top_xcorr_raw_indices(top_table, EDMD_outputs, top_n)
%RESOLVE_BOLD_TOP_XCORR_RAW_INDICES Map ranked top-table rows to unique raw Koopman indices.

if nargin < 3 || isempty(top_n)
    top_n = inf;
end

raw_indices = nan(1, 0);
top_row_indices = nan(1, 0);
if isempty(top_table)
    raw_indices = raw_indices(:).';
    top_row_indices = top_row_indices(:).';
    return;
end

n = height(top_table);
for ii = 1:n
    if numel(raw_indices) >= top_n
        break;
    end
    mode_idx = top_table.bold_mode_index(ii);
    raw_idx = local_mode_to_raw_index(mode_idx, EDMD_outputs);
    if ~isfinite(raw_idx) || raw_idx < 1
        continue;
    end
    if raw_idx > size(EDMD_outputs.kpm_modes, 1)
        continue;
    end
    if any(raw_indices == raw_idx)
        continue;
    end
    raw_indices(end + 1) = round(raw_idx); %#ok<AGROW>
    top_row_indices(end + 1) = ii; %#ok<AGROW>
end

raw_indices = raw_indices(:).';
top_row_indices = top_row_indices(:).';
end


function raw_idx = local_mode_to_raw_index(mode_idx, EDMD_outputs)
mode_idx = round(double(mode_idx));
n_plot_modes = size(EDMD_outputs.kpm_modes, 1);
n_current_modes = NaN;
if isfield(EDMD_outputs, 'efuns') && ~isempty(EDMD_outputs.efuns)
    n_current_modes = size(EDMD_outputs.efuns, 2);
elseif isfield(EDMD_outputs, 'evalues') && ~isempty(EDMD_outputs.evalues)
    n_current_modes = numel(EDMD_outputs.evalues);
end

idx_original = NaN;
if isfield(EDMD_outputs, 'idx_final_in_original') && ...
        numel(EDMD_outputs.idx_final_in_original) >= mode_idx && ...
        ~isempty(EDMD_outputs.idx_final_in_original(mode_idx))
    idx_original = double(EDMD_outputs.idx_final_in_original(mode_idx));
end

if isfinite(n_current_modes) && n_plot_modes == n_current_modes && ...
        mode_idx >= 1 && mode_idx <= n_plot_modes
    raw_idx = mode_idx;
elseif isfinite(idx_original) && idx_original >= 1 && idx_original <= n_plot_modes
    raw_idx = idx_original;
elseif mode_idx >= 1 && mode_idx <= n_plot_modes
    raw_idx = mode_idx;
else
    raw_idx = NaN;
end
end
