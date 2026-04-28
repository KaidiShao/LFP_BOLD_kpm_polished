function params = build_blp_consensus_state_pipeline_params(params)
%BUILD_BLP_CONSENSUS_STATE_PIPELINE_PARAMS Apply canonical defaults for the
% mainline BLP consensus-state pipeline.

if nargin < 1 || isempty(params)
    params = struct();
end

if ~isfield(params, 'event_params') || isempty(params.event_params)
    params.event_params = struct();
end
params.event_params = local_set_event_defaults(params.event_params);

if ~isfield(params, 'density_params') || isempty(params.density_params)
    params.density_params = struct();
end
params.density_params = local_set_density_defaults(params.density_params);

if ~isfield(params, 'consensus_params') || isempty(params.consensus_params)
    params.consensus_params = struct();
end
params.consensus_params = local_set_consensus_defaults(params.consensus_params);

if ~isfield(params, 'summary_params') || isempty(params.summary_params)
    params.summary_params = struct();
end
params.summary_params = local_set_summary_defaults(params.summary_params);

if ~isfield(params, 'window_params') || isempty(params.window_params)
    params.window_params = struct();
end
params.window_params = local_set_window_defaults(params.window_params);

if ~isfield(params, 'plot_params') || isempty(params.plot_params)
    params.plot_params = struct();
end
params.plot_params = local_set_plot_defaults(params.plot_params);
end


function p = local_set_event_defaults(p)
if ~isfield(p, 'passband'), p.passband = [2, 15; 30, 90; 90, 190]; end
if ~isfield(p, 'band_labels'), p.band_labels = {'theta', 'gamma', 'ripple'}; end
if ~isfield(p, 'L_start_range'), p.L_start_range = [151, 101, 51]; end
if ~isfield(p, 'L_extract_range'), p.L_extract_range = [301, 201, 101]; end
if ~isfield(p, 'ThresRatio_range'), p.ThresRatio_range = [3.5, 4, 4]; end
if ~isfield(p, 'input_normalization'), p.input_normalization = 'zscore_per_channel'; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = local_set_density_defaults(p)
if ~isfield(p, 'bin_sec'), p.bin_sec = 2; end
if ~isfield(p, 'smooth_sigma_sec'), p.smooth_sigma_sec = 2; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = local_set_consensus_defaults(p)
if ~isfield(p, 'min_channel_count'), p.min_channel_count = []; end
if ~isfield(p, 'require_region_presence'), p.require_region_presence = false; end
if ~isfield(p, 'required_regions') || isempty(p.required_regions), p.required_regions = {'hp', 'pl'}; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = local_set_summary_defaults(p)
if ~isfield(p, 'save_csv'), p.save_csv = true; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = local_set_window_defaults(p)
if ~isfield(p, 'window_length_samples'), p.window_length_samples = 6000; end
if ~isfield(p, 'window_mode') || isempty(p.window_mode), p.window_mode = 'global'; end
if ~isfield(p, 'keep_partial_window'), p.keep_partial_window = false; end
if ~isfield(p, 'top_k'), p.top_k = 30; end
if ~isfield(p, 'save_csv'), p.save_csv = true; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = local_set_plot_defaults(p)
if ~isfield(p, 'enable'), p.enable = true; end
if ~isfield(p, 'save_png'), p.save_png = true; end
if ~isfield(p, 'close_after_save'), p.close_after_save = true; end
if ~isfield(p, 'skip_existing'), p.skip_existing = false; end
if ~isfield(p, 'fallback_size_px') || isempty(p.fallback_size_px), p.fallback_size_px = [4979, 2888]; end
if ~isfield(p, 'reference_dir'), p.reference_dir = ''; end
if ~isfield(p, 'resolution') || isempty(p.resolution), p.resolution = 220; end
if ~isfield(p, 'event_colors') || isempty(p.event_colors)
    p.event_colors = [ ...
        0.0000, 0.4470, 0.7410; ...
        0.8500, 0.3250, 0.0980; ...
        0.4660, 0.6740, 0.1880];
end
if ~isfield(p, 'freq_range_to_plot') || isempty(p.freq_range_to_plot), p.freq_range_to_plot = [0, 250]; end
if ~isfield(p, 'color_limits'), p.color_limits = []; end
if ~isfield(p, 'spec_colormap') || isempty(p.spec_colormap)
    if exist('othercolor', 'file') == 2
        p.spec_colormap = flipud(othercolor('Spectral10'));
    else
        p.spec_colormap = flipud(turbo(256));
    end
end
end
