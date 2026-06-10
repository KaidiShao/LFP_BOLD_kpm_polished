function params = apply_bold_cross_modal_output_mode(params, output_mode)
%APPLY_BOLD_CROSS_MODAL_OUTPUT_MODE Configure pipeline 8 output grouping.
%
% output_mode:
%   'combined'  -> only global top tables / activation maps
%   'separate'  -> only per-density top tables / activation maps
%   'both'      -> export both combined and per-density outputs

if nargin < 1 || isempty(params)
    params = build_bold_cross_modal_coupling_params();
end
if nargin < 2 || isempty(output_mode)
    output_mode = local_get_field(params, 'output_mode', 'both');
end

mode = lower(strtrim(char(string(output_mode))));
switch mode
    case 'combined'
        xcorr_combined = true;
        xcorr_by_density = false;
        act_combined = true;
        act_by_density = false;
    case 'separate'
        xcorr_combined = false;
        xcorr_by_density = true;
        act_combined = false;
        act_by_density = true;
    case 'both'
        xcorr_combined = true;
        xcorr_by_density = true;
        act_combined = true;
        act_by_density = true;
    otherwise
        error('Unsupported pipeline 8 output_mode: %s', output_mode);
end

params.output_mode = mode;
params.xcorr.export_combined = xcorr_combined;
params.xcorr.export_by_density = xcorr_by_density;
params.activation.export_combined = act_combined;
params.activation.export_by_density = act_by_density;
end


function value = local_get_field(S, name, default_value)
if isstruct(S) && isfield(S, name) && ~isempty(S.(name))
    value = S.(name);
else
    value = default_value;
end
end
