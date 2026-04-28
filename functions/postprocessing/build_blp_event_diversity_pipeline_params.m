function params = build_blp_event_diversity_pipeline_params(params)
%BUILD_BLP_EVENT_DIVERSITY_PIPELINE_PARAMS Shared defaults for the
% event-diversity quick-check branch.

if nargin < 1 || isempty(params)
    params = struct();
end

shared = build_blp_consensus_state_pipeline_params(params);

if isfield(params, 'window_params') && isstruct(params.window_params)
    window_params = params.window_params;
else
    window_params = struct();
end

window_params = local_set_window_defaults(window_params);

params = struct();
params.event_params = shared.event_params;
params.density_params = shared.density_params;
params.consensus_params = shared.consensus_params;
params.summary_params = shared.summary_params;
params.window_params = window_params;
end


function params = local_set_window_defaults(params)
if ~isfield(params, 'window_length_samples')
    params.window_length_samples = 6000;
end

if ~isfield(params, 'window_mode') || isempty(params.window_mode)
    params.window_mode = 'global';
end

if ~isfield(params, 'keep_partial_window')
    params.keep_partial_window = false;
end

if ~isfield(params, 'top_k')
    params.top_k = 10;
end

if ~isfield(params, 'save_csv')
    params.save_csv = true;
end

if ~isfield(params, 'force_recompute')
    params.force_recompute = false;
end
end
