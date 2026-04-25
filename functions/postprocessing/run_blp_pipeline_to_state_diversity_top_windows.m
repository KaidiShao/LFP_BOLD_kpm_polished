function out = run_blp_pipeline_to_state_diversity_top_windows(cfg, output_root, params)
%RUN_BLP_PIPELINE_TO_STATE_DIVERSITY_TOP_WINDOWS Run the event line to plots.

if nargin < 2 || isempty(output_root)
    output_root = get_project_processed_root();
end

if nargin < 3
    params = struct();
end

params = apply_default_params(params);

fprintf('Running state-diversity pipeline for %s\n', cfg.dataset_id);
fprintf('Raw data root:\n  %s\n', cfg.raw_data_root);
fprintf('Output root:\n  %s\n', output_root);

D = load_blp_dataset(cfg);
fprintf('Loaded concatenated raw data: [%d, %d]\n', size(D.data, 1), size(D.data, 2));

R = compute_blp_bandpass_events(D, cfg, output_root, params.event_params);
fprintf('Saved event detection result to:\n  %s\n', R.save_file);

E = compute_blp_event_density(cfg, output_root, R, params.density_params, R.save_file);
fprintf('Saved event-density result to:\n  %s\n', E.save_file);

C = compute_blp_consensus_states(cfg, output_root, R, params.consensus_params, R.save_file);
fprintf('Saved consensus-state result to:\n  %s\n', C.save_file);

S = summarize_blp_consensus_state_types(cfg, output_root, C, params.summary_params, C.save_file);
fprintf('Saved consensus-state summary to:\n  %s\n', S.save_file);

Spec = compute_blp_region_spectrograms_streamed(D, cfg, output_root, params.spectrogram_opts);
fprintf('Prepared spectrogram result in:\n  %s\n', Spec.abs_file);

W = analyze_blp_consensus_state_diversity_windows(cfg, output_root, C, params.window_params, C.save_file);
fprintf('Saved state-diversity result to:\n  %s\n', W.save_file);

P = export_top_consensus_state_diversity_window_plots(cfg, output_root, W, params.plot_params);
fprintf('Saved top-window plots to:\n  %s\n', P.save_dir);

out = struct();
out.cfg = cfg;
out.output_root = output_root;
out.D = D;
out.R = R;
out.E = E;
out.C = C;
out.S = S;
out.Spec = Spec;
out.W = W;
out.P = P;
out.params = params;
end


function params = apply_default_params(params)
if ~isfield(params, 'event_params') || isempty(params.event_params)
    params.event_params = struct();
end
params.event_params = set_event_defaults(params.event_params);

if ~isfield(params, 'density_params') || isempty(params.density_params)
    params.density_params = struct();
end
params.density_params = set_density_defaults(params.density_params);

if ~isfield(params, 'consensus_params') || isempty(params.consensus_params)
    params.consensus_params = struct();
end
params.consensus_params = set_consensus_defaults(params.consensus_params);

if ~isfield(params, 'summary_params') || isempty(params.summary_params)
    params.summary_params = struct();
end
params.summary_params = set_summary_defaults(params.summary_params);

if ~isfield(params, 'spectrogram_opts') || isempty(params.spectrogram_opts)
    params.spectrogram_opts = struct();
end
params.spectrogram_opts = set_spectrogram_defaults(params.spectrogram_opts);

if ~isfield(params, 'window_params') || isempty(params.window_params)
    params.window_params = struct();
end
params.window_params = set_window_defaults(params.window_params);

if ~isfield(params, 'plot_params') || isempty(params.plot_params)
    params.plot_params = struct();
end
end


function p = set_event_defaults(p)
if ~isfield(p, 'passband'), p.passband = [2, 15; 30, 90; 90, 190]; end
if ~isfield(p, 'band_labels'), p.band_labels = {'theta', 'gamma', 'ripple'}; end
if ~isfield(p, 'L_start_range'), p.L_start_range = [151, 101, 51]; end
if ~isfield(p, 'L_extract_range'), p.L_extract_range = [301, 201, 101]; end
if ~isfield(p, 'ThresRatio_range'), p.ThresRatio_range = [3.5, 4, 4]; end
if ~isfield(p, 'input_normalization'), p.input_normalization = 'zscore_per_channel'; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = set_density_defaults(p)
if ~isfield(p, 'bin_sec'), p.bin_sec = 2; end
if ~isfield(p, 'smooth_sigma_sec'), p.smooth_sigma_sec = 2; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = set_consensus_defaults(p)
if ~isfield(p, 'min_channel_count'), p.min_channel_count = []; end
if ~isfield(p, 'require_region_presence'), p.require_region_presence = false; end
if ~isfield(p, 'required_regions') || isempty(p.required_regions), p.required_regions = {'hp', 'pl'}; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = set_summary_defaults(p)
if ~isfield(p, 'save_csv'), p.save_csv = true; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end


function p = set_spectrogram_defaults(p)
if ~isfield(p, 'save_precision') || isempty(p.save_precision), p.save_precision = 'single'; end
if ~isfield(p, 'return_data') || isempty(p.return_data), p.return_data = false; end
if ~isfield(p, 'force_recompute') || isempty(p.force_recompute), p.force_recompute = false; end
end


function p = set_window_defaults(p)
if ~isfield(p, 'window_length_samples'), p.window_length_samples = 6000; end
if ~isfield(p, 'window_mode') || isempty(p.window_mode), p.window_mode = 'global'; end
if ~isfield(p, 'keep_partial_window'), p.keep_partial_window = false; end
if ~isfield(p, 'top_k'), p.top_k = 30; end
if ~isfield(p, 'save_csv'), p.save_csv = true; end
if ~isfield(p, 'force_recompute'), p.force_recompute = false; end
end
