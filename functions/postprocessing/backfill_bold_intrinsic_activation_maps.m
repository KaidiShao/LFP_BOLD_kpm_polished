function [act_result, attempted] = backfill_bold_intrinsic_activation_maps(post_file, params)
%BACKFILL_BOLD_INTRINSIC_ACTIVATION_MAPS Export missing intrinsic maps from an existing BOLD_POST artifact.

attempted = false;
act_result = struct();

if nargin < 1 || isempty(post_file)
    return;
end
if nargin < 2 || ~isstruct(params) || ...
        ~isfield(params, 'intrinsic_activation') || ~isstruct(params.intrinsic_activation) || ...
        ~isfield(params.intrinsic_activation, 'enabled') || ~params.intrinsic_activation.enabled
    return;
end
if exist(post_file, 'file') ~= 2
    return;
end

act_params = params.intrinsic_activation;
act_params.processed_root = params.processed_root;
act_params.datapons_root = params.datapons_root;
act_params.output_root = local_resolve_output_root(post_file);
attempted = true;
[act_result, ~] = export_bold_intrinsic_activation_maps(post_file, act_params);
end


function out_root = local_resolve_output_root(post_file)
mat_dir = fileparts(char(string(post_file)));
out_root = fileparts(mat_dir);
end
