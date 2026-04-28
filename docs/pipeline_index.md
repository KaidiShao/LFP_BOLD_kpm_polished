# Pipeline Index

This note defines the project using a pipeline view rather than only the
technical folder view (`io`, `preprocessing`, `postprocessing`, `plottings`).

The goal is to answer four practical questions:

1. How many real pipelines does this repo contain?
2. Which files belong to each pipeline?
3. Which files are canonical entry points?
4. Which files are shared infrastructure rather than pipeline-owned logic?

## Pipeline Count

For maintenance purposes, treat this repo as:

- 1 overall analysis framework
- 7 operational pipelines
- 1 optional comparison branch

The seven operational pipelines are:

1. BLP dictionary pipeline
2. BLP event-state pipeline
3. BOLD observable pipeline
4. ResKoopNet training/export pipeline
5. BLP eigenfunction interpretation pipeline
6. BOLD postprocessing pipeline
7. BLP-BOLD cross-modal coupling pipeline

Optional branch:

- BLP CNN comparison branch

## Pipeline Summary

| Pipeline | Purpose | Canonical entry points | Main outputs |
|---|---|---|---|
| BLP dictionary pipeline | Turn raw BLP into ResKoopNet-ready observables. | `scripts/script_preprocess_e10gb1_to_observables_streamed.m`, `docs/raw_dataset_to_reskoopnet_pipeline.md` | spectrogram MAT files, ResKoopNet dictionary MAT files, observable CSVs |
| BLP event-state pipeline | Turn raw BLP into events, densities, consensus states, and diversity windows. | `scripts/script_run_one_cfg_to_consensus_state_top_windows.m`, `scripts/script_run_cfgs_to_consensus_state_top_windows.m` | event MAT files, density MAT files, consensus-state MAT files, summary CSVs, top-window tables and plots |
| BOLD observable pipeline | Turn raw BOLD ROI time series into standardized observable packages. | `scripts/script_build_one_cfg_bold_observables.m`, `scripts/script_build_cfgs_bold_observables.m` | per-mode BOLD observable MAT files, session-aware snapshots |
| ResKoopNet training/export pipeline | Train MLP ResKoopNet models and export EDMD-style outputs. | `python_scripts/autodl/run_autodl_reskoopnet_mlp.py`, `docs/autodl_automatic_pipeline_setup_2026-04-18.md` | chunked EDMD output MAT files, summary MAT files, checkpoints, manifests |
| BLP eigenfunction interpretation pipeline | Interpret trained BLP MLP outputs using DR, top windows, and spike alignment. | `scripts/script_run_all_complete_mlp_efun_dr_pipeline.m`, `functions/postprocessing/run_eigenfunction_reduction_pipeline.m`, `functions/postprocessing/run_mlp_top_state_diversity_postprocessing_pipeline.m` | reduced eigenfunction results, state-space plots, top-window figures, spike correlation outputs |
| BOLD postprocessing pipeline | Interpret trained BOLD MLP outputs using postprocessed eigenfunctions, deconvolution, and timescales. | `scripts/script_postprocess_bold_reskoopnet_results.m`, `functions/postprocessing/postprocess_bold_reskoopnet_results.m` | BOLD post MAT files, deconvolution outputs, timescale figures, activation maps |
| BLP-BOLD cross-modal coupling pipeline | Quantify lagged coupling between BLP densities/states and BOLD eigenfunction dynamics. | `scripts/script_correlate_bold_efuns_with_densities.m`, `functions/postprocessing/compute_bold_efun_density_cross_correlation.m` | cross-correlation MAT files, peak tables, summary figures |

## Pipeline Ownership By Folder

### Shared infrastructure

These files are shared and should not be treated as belonging to only one
pipeline:

- `functions/io/+io_raw/*`
- `functions/io/+io_edmd/*`
- `functions/io/+io_results/*`
- `functions/io/+io_project/*`
- `functions/io/+io_utils/*`

The top-level `functions/io/*.m` files are migration stubs that now error on
old unqualified names.

### 1. BLP dictionary pipeline

Primary scripts:

- `scripts/script_preprocess_e10gb1_to_observables_streamed.m`
- `scripts/script_build_reskoopnet_dicts.m`

Primary preprocessing functions:

- `functions/preprocessing/compute_blp_region_spectrograms_streamed.m`
- `functions/preprocessing/compute_blp_region_spectrograms.m`
- `functions/preprocessing/build_reskoopnet_dicts.m`

Supporting plotting:

- `functions/plottings/plot_blp_segment_with_spectrogram.m`
- `functions/plottings/prepare_blp_plot_data.m`

### 2. BLP event-state pipeline

Primary scripts:

- `scripts/script_run_one_cfg_to_consensus_state_top_windows.m`
- `scripts/script_run_cfgs_to_consensus_state_top_windows.m`
- `scripts/script_run_e10gb1_to_consensus_state_top_windows.m`
- `scripts/script_run_one_cfg_to_event_diversity_windows.m`
- `scripts/script_run_cfgs_to_event_diversity_windows.m`
- `scripts/script_run_e10gb1_to_event_diversity_windows.m`
- `scripts/script_compute_one_cfg_blp_bandpass_events.m`
- `scripts/script_compute_one_cfg_blp_consensus_states.m`
- `scripts/script_compute_one_cfg_blp_event_density.m`
- `scripts/script_analyze_one_cfg_blp_consensus_state_diversity_windows.m`
- `scripts/script_plot_top_consensus_event_diversity_windows.m`
- `scripts/script_plot_top_consensus_event_diversity_windows_e10gb1.m`
- `scripts/script_plot_top_consensus_state_diversity_windows.m`
- `scripts/script_plot_top_consensus_state_diversity_windows_e10gb1.m`

Primary functions:

- `functions/preprocessing/compute_blp_bandpass_events.m`
- `functions/postprocessing/compute_blp_event_density.m`
- `functions/postprocessing/compute_blp_consensus_states.m`
- `functions/postprocessing/summarize_blp_consensus_state_types.m`
- `functions/postprocessing/analyze_blp_consensus_event_diversity_windows.m`
- `functions/postprocessing/analyze_blp_consensus_state_diversity_windows.m`

Supporting plotting:

- `functions/plottings/export_top_consensus_state_diversity_window_plots.m`
- `functions/plottings/build_blp_plot_window_cache.m`
- `functions/plottings/plot_blp_segment_with_events.m`
- `functions/plottings/plot_blp_segment_with_spectrogram_and_koopman.m`

### 3. BOLD observable pipeline

Primary scripts:

- `scripts/script_build_one_cfg_bold_observables.m`
- `scripts/script_build_cfgs_bold_observables.m`
- `scripts/script_build_e10gb1_bold_observables.m`
- `scripts/script_plot_bold_pre_reskoopnet_qc.m`
- `scripts/script_plot_bold_svd100_qc.m`

Primary preprocessing functions:

- `functions/preprocessing/preprocess_bold_sessions.m`
- `functions/preprocessing/build_bold_observables.m`
- `functions/preprocessing/build_bold_slow_band_power_observables.m`
- `functions/preprocessing/make_session_aware_snapshot_pairs.m`

Supporting plotting:

- `functions/plottings/plot_bold_pre_reskoopnet_qc.m`
- `functions/plottings/plot_bold_segment.m`

### 4. ResKoopNet training/export pipeline

Primary docs and launchers:

- `docs/autodl_automatic_pipeline_setup_2026-04-18.md`
- `python_scripts/autodl/run_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/multi_dataset_batch_controller_autodl_reskoopnet_mlp.py`

This pipeline consumes outputs from the BLP dictionary pipeline and the BOLD
observable pipeline, then produces chunked EDMD-like outputs used by the
downstream postprocessing pipelines.

### 5. BLP eigenfunction interpretation pipeline

Primary scripts:

- `scripts/script_run_all_complete_mlp_efun_dr_pipeline.m`
- `scripts/script_run_all_complete_mlp_efun_dr_umap_only.m`
- `scripts/script_run_e10gb1_full_mlp_efun_dr_pipeline.m`
- `scripts/script_run_e10gb1_eigenfunction_reduction_method_sweep.m`
- `scripts/script_run_eigenfunction_reduction_minimal.m`
- `scripts/script_run_completed_mlp_top_state_diversity_postprocessing.m`
- `scripts/script_run_completed_mlp_added_postprocessing_stages.m`
- `scripts/script_run_missing_completed_model_postprocessing.m`
- `scripts/script_compare_e10gb1_spike_residual.m`
- `scripts/script_plot_e10gb1_spike_residual_correlation.m`
- `scripts/script_plot_e10gb1_spike_residual_lagged_correlation.m`
- `scripts/script_analyze_e10gb1_efun_peaks_by_state.m`

Primary functions:

- `functions/postprocessing/run_eigenfunction_reduction_pipeline.m`
- `functions/postprocessing/reduce_eigenfunction_time_path.m`
- `functions/postprocessing/reduce_eigenfunction_spectrum_path.m`
- `functions/postprocessing/get_thresholded_density.m`
- `functions/postprocessing/get_thresholded_events.m`
- `functions/postprocessing/get_dimred_thresholded_density.m`
- `functions/postprocessing/get_dimred_thresholded_events.m`
- `functions/postprocessing/run_mlp_top_state_diversity_postprocessing_pipeline.m`
- `functions/postprocessing/run_reskoopnet_eigenfunction_density_batch.m`
- `functions/postprocessing/compute_spike_residual_comparison.m`
- `functions/postprocessing/compute_spike_residual_lagged_correlation.m`
- `functions/postprocessing/analyze_eigenfunction_component_peaks_by_consensus_state.m`
- `functions/postprocessing/fit_empirical_eigenvalues_from_efuns.m`

Supporting plotting:

- `functions/plottings/plot_eigenfunction_component_overview.m`
- `functions/plottings/plot_eigenfunction_spectrum_diagnostics.m`
- `functions/plottings/plot_eigenfunction_state_space_trajectory.m`
- `functions/plottings/plot_eigenfunction_state_space_consensus_trajectory.m`
- `functions/plottings/plot_acf_and_theoretical.m`

### 6. BOLD postprocessing pipeline

Primary scripts:

- `scripts/script_postprocess_bold_reskoopnet_results.m`
- `scripts/script_generate_missing_bold_activation_maps.m`
- `scripts/script_run_one_bold_reskoopnet_post_xcorr_activation_maps.m`
- `scripts/plot_bold_reskoopnet_mode_activation_map.m`
- `scripts/plot_bold_reskoopnet_mode_activation_map_reference_style.m`

Primary functions:

- `functions/postprocessing/postprocess_bold_reskoopnet_results.m`
- `functions/postprocessing/postprocess_EDMD_outputs.m`
- `functions/postprocessing/postprocess_EDMD_outputs_deconv_efuns.m`
- `functions/postprocessing/postprocess_EDMD_outputs_timescale.m`

### 7. BLP-BOLD cross-modal coupling pipeline

Primary scripts:

- `scripts/script_correlate_bold_efuns_with_densities.m`
- `scripts/script_run_one_bold_reskoopnet_post_xcorr_activation_maps.m`

Primary functions:

- `functions/postprocessing/compute_bold_efun_density_cross_correlation.m`
- `functions/plottings/plot_bold_efun_density_cross_correlation_summary.m`

Inputs usually come from:

- BLP event-state pipeline outputs such as event density or state diversity
- BOLD postprocessing pipeline outputs such as `BOLD_POST` MAT files

The single-run script
`script_run_one_bold_reskoopnet_post_xcorr_activation_maps.m` spans both the
BOLD postprocessing pipeline and the cross-modal coupling pipeline.

### Optional BLP CNN comparison branch

This remains a comparison branch rather than a primary pipeline. The framework
diagram includes it, but the current MATLAB organization is centered on the MLP
line.

## Practical Use

When adding or refactoring code, classify each new file by two labels:

1. pipeline owner
2. role inside that pipeline

Use these role labels:

- entry script
- stage function
- shared IO/helper
- plotting/readout
- training/export

This should be done before deciding which folder the file belongs to.

## Immediate Implication For Postprocessing

`functions/postprocessing` is not one coherent pipeline. It currently mixes at
least four different downstream lines:

- BLP event-state downstream
- BLP eigenfunction interpretation
- BOLD postprocessing
- BLP-BOLD cross-modal coupling

That is why the folder feels structurally noisy even when many individual files
are internally reasonable.
