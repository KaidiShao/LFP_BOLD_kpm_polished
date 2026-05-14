function [summary_result, attempted] = backfill_bold_intrinsic_roi_bar_summaries(post_file, params)
%BACKFILL_BOLD_INTRINSIC_ROI_BAR_SUMMARIES Export missing intrinsic ROI summaries from an existing BOLD_POST.

attempted = false;
summary_result = struct();

if nargin < 1 || isempty(post_file)
    return;
end
if nargin < 2 || ~isstruct(params) || ...
        ~isfield(params, 'intrinsic_roi_summary') || ~isstruct(params.intrinsic_roi_summary) || ...
        ~isfield(params.intrinsic_roi_summary, 'enabled') || ~params.intrinsic_roi_summary.enabled
    return;
end
if exist(post_file, 'file') ~= 2
    return;
end

summary_params = params.intrinsic_roi_summary;
summary_params.processed_root = params.processed_root;
summary_params.datapons_root = params.datapons_root;
summary_params.output_root = local_resolve_output_root(post_file);
attempted = true;
[summary_result, ~] = export_bold_intrinsic_roi_bar_summaries(post_file, summary_params);
end


function out_root = local_resolve_output_root(post_file)
mat_dir = fileparts(char(string(post_file)));
out_root = fileparts(mat_dir);
end
