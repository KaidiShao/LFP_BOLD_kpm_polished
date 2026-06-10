function result = compute_spkt_residual_cross_correlation(cfg, params)
%COMPUTE_SPKT_RESIDUAL_CROSS_CORRELATION Compute SPKT residual cross-correlation.

if nargin < 1 || isempty(cfg)
    error('cfg must be provided.');
end

if nargin < 2
    params = struct();
end

[ctx, params] = load_spkt_residual_cross_correlation_inputs(cfg, params);
[session_state, total_streamed_length, k_used] = align_spkt_residual_to_spike_resolution(ctx);

expected_required_length = max([session_state.lfp_end_idx]);
if total_streamed_length < expected_required_length
    error(['Stopped after %d samples, but the selected sessions require at least %d. ', ...
        'The EDMD source appears incomplete for the chosen sessions.'], ...
        total_streamed_length, expected_required_length);
end

result = build_spkt_residual_cross_correlation_result( ...
    cfg, params, ctx, session_state, total_streamed_length, k_used);
result.save_paths = empty_residual_cross_correlation_save_paths();
end
