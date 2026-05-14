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
- 8 operational pipelines
- 1 optional comparison branch

The eight operational pipelines are:

1. BLP dictionary pipeline
2. BLP event-state pipeline
3. BOLD observable pipeline
4. ResKoopNet training/export pipeline
5. BLP eigenfunction dimension-reduction pipeline
6. BLP eigenfunction postprocessing pipeline
7. BOLD postprocessing pipeline
8. BLP-BOLD cross-modal coupling pipeline

Optional branch:

- BLP CNN comparison branch

## Pipeline Summary

| Pipeline | Purpose | Canonical entry points | Main outputs |
|---|---|---|---|
| BLP dictionary pipeline | Turn raw BLP into ResKoopNet-ready observables. | `scripts/script_preprocess_e10gb1_to_observables_streamed.m`, `docs/raw_dataset_to_reskoopnet_pipeline.md` | spectrogram MAT files, ResKoopNet dictionary MAT files, observable CSVs |
| BLP event-state pipeline | Turn raw BLP into events, densities, consensus states, and diversity windows. | `scripts/script_run_one_cfg_to_consensus_state_top_windows.m`, `scripts/script_run_cfgs_to_consensus_state_top_windows.m` | event MAT files, density MAT files, consensus-state MAT files, summary CSVs, top-window tables and plots |
| BOLD observable pipeline | Turn raw BOLD ROI time series into standardized observable packages. | `scripts/script_build_one_cfg_bold_observables.m`, `scripts/script_build_cfgs_bold_observables.m` | per-mode BOLD observable MAT files, session-aware snapshots |
| ResKoopNet training/export pipeline | Train MLP ResKoopNet models and export EDMD-style outputs. Treat this as `pipeline 4.1` for BLP and `pipeline 4.2` for BOLD. | `python_scripts/autodl/run_blp_observables_mlp_reskoopnet.ps1`, `python_scripts/local/run_blp_observables_mlp_reskoopnet_local.ps1`, `python_scripts/autodl/run_bold_observables_mlp_reskoopnet.ps1`, `python_scripts/local/run_bold_observables_mlp_reskoopnet_local.ps1` | chunked EDMD output MAT files, summary MAT files, checkpoints, manifests |
| BLP eigenfunction dimension-reduction pipeline | Reduce trained BLP MLP eigenfunction outputs and derive thresholded density/event summaries from full or reduced coordinates. | `scripts/script_run_one_cfg_blp_eigenfunction_reduction.m`, `scripts/script_run_cfgs_blp_eigenfunction_reduction.m`, `functions/postprocessing/run_eigenfunction_reduction_pipeline.m` | reduced eigenfunction MAT files, thresholded density/event MAT files, state-space plots, DR summaries, peak-by-state tables/figures |
| BLP eigenfunction postprocessing pipeline | Consume reduced BLP eigenfunction outputs for top-window interpretation, residual-vs-MUA/SPKT cross-correlation, and state-linked readouts. | `scripts/script_run_one_cfg_blp_eigenfunction_postprocessing.m`, `scripts/script_run_cfgs_blp_eigenfunction_postprocessing.m` | top-window figures, residual-vs-MUA cross-correlation tables/figures, residual-vs-SPKT cross-correlation tables/figures, lagged cross-correlation tables/figures |
| BOLD postprocessing pipeline | Interpret trained BOLD MLP outputs using postprocessed eigenfunctions, deconvolution, timescales, and intrinsic mode maps. | `scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m`, `scripts/script_postprocess_bold_reskoopnet_results.m` | BOLD post MAT files, deconvolution outputs, timescale figures, intrinsic BOLD mode activation maps |
| BLP-BOLD cross-modal coupling pipeline | Quantify lagged coupling between BLP densities/states and BOLD eigenfunction dynamics. | `scripts/script_run_one_cfg_bold_cross_modal_coupling.m`, `scripts/script_correlate_bold_efuns_with_densities.m` | cross-correlation MAT files, peak tables, summary figures, xcorr-ranked BOLD activation maps |

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
- `scripts/script_plot_one_cfg_bold_pre_reskoopnet_qc.m`
- `scripts/script_plot_cfgs_bold_pre_reskoopnet_qc.m`
- `scripts/script_plot_e10gb1_bold_pre_reskoopnet_qc.m`
- `scripts/script_report_bold_observable_shapes.m`

Primary preprocessing functions:

- `functions/preprocessing/preprocess_bold_sessions.m`
- `functions/preprocessing/build_bold_observables.m`
- `functions/preprocessing/build_bold_slow_band_power_observables.m`
- `functions/preprocessing/make_session_aware_snapshot_pairs.m`

Supporting plotting:

- `functions/plottings/plot_bold_pre_reskoopnet_qc.m`
- `functions/plottings/plot_bold_segment.m`

### 4. ResKoopNet training/export pipeline

Read this family as two current branches:

- `pipeline 4.1`: BLP training/export
- `pipeline 4.2`: BOLD training/export

Primary docs and launchers:

- `docs/autodl_automatic_pipeline_setup_2026-04-18.md`
- `python_scripts/autodl/run_blp_observables_mlp_reskoopnet.ps1`
- `python_scripts/local/run_blp_observables_mlp_reskoopnet_local.ps1`
- `python_scripts/local/run_blp_observables_mlp_reskoopnet_wsl.ps1`
- `python_scripts/local/run_blp_observables_mlp_reskoopnet_wsl.sh`
- `python_scripts/autodl/run_bold_observables_mlp_reskoopnet.ps1`
- `python_scripts/local/run_bold_observables_mlp_reskoopnet_local.ps1`
- `python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.ps1`
- `python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.sh`
- `python_scripts/autodl/run_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/multi_dataset_batch_controller_autodl_reskoopnet_mlp.py`

This family consumes outputs from the BLP dictionary pipeline (`pipeline 4.1`)
and the BOLD observable pipeline (`pipeline 4.2`), then produces chunked
EDMD-like outputs used by the downstream postprocessing pipelines.

### 5. BLP eigenfunction dimension-reduction pipeline

Primary scripts:

- `scripts/script_run_one_cfg_blp_eigenfunction_reduction.m`
- `scripts/script_run_cfgs_blp_eigenfunction_reduction.m`
- `scripts/script_run_e10gb1_blp_eigenfunction_reduction.m`
- `scripts/script_analyze_e10gb1_efun_peaks_by_state.m`
- `scripts/script_plot_e10gb1_efun_svd_top30_state_diversity_windows.m`

Primary functions:

- `functions/postprocessing/build_blp_eigenfunction_reduction_params.m`
- `functions/postprocessing/build_blp_eigenfunction_reduction_method_specs.m`
- `functions/postprocessing/load_blp_eigenfunction_reduction_source.m`
- `functions/postprocessing/build_blp_eigenfunction_reduction_cfg.m`
- `functions/postprocessing/plot_blp_eigenfunction_reduction_outputs.m`
- `functions/postprocessing/run_eigenfunction_reduction_pipeline.m`
- `functions/postprocessing/reduce_eigenfunction_time_path.m`
- `functions/postprocessing/reduce_eigenfunction_spectrum_path.m`
- `functions/postprocessing/get_thresholded_density.m`
- `functions/postprocessing/get_thresholded_events.m`
- `functions/postprocessing/get_dimred_thresholded_density.m`
- `functions/postprocessing/get_dimred_thresholded_events.m`
- `functions/postprocessing/fit_empirical_eigenvalues_from_efuns.m`
- `functions/postprocessing/analyze_eigenfunction_component_peaks_by_consensus_state.m`

Supporting plotting:

- `functions/plottings/plot_eigenfunction_component_overview.m`
- `functions/plottings/plot_eigenfunction_spectrum_diagnostics.m`
- `functions/plottings/plot_eigenfunction_state_space_trajectory.m`
- `functions/plottings/plot_eigenfunction_state_space_consensus_trajectory.m`
- `functions/plottings/plot_acf_and_theoretical.m`

### 6. BLP eigenfunction postprocessing pipeline

Primary scripts:

- `scripts/script_run_one_cfg_blp_eigenfunction_postprocessing.m`
- `scripts/script_run_cfgs_blp_eigenfunction_postprocessing.m`
- `scripts/script_run_e10gb1_blp_eigenfunction_postprocessing.m`
- `scripts/script_run_completed_mlp_top_state_diversity_postprocessing.m`
- `scripts/script_run_completed_mlp_added_postprocessing_stages.m`
- `scripts/script_run_missing_completed_model_postprocessing.m`
- `scripts/script_compare_e10gb1_mua_residual_cross_correlation.m`
- `scripts/script_plot_e10gb1_mua_residual_cross_correlation.m`
- `scripts/script_compare_e10gb1_spkt_residual_cross_correlation.m`
- `scripts/script_plot_e10gb1_spkt_residual_cross_correlation.m`
- `scripts/script_plot_e10gb1_spkt_residual_lagged_cross_correlation.m`

Primary functions:

- `functions/postprocessing/compute_mua_residual_cross_correlation.m`
- `functions/postprocessing/compute_spkt_residual_cross_correlation.m`
- `functions/postprocessing/compute_spkt_residual_lagged_cross_correlation.m`

Supporting plotting:

- `functions/plottings/plot_eigenfunction_state_space_consensus_trajectory.m`
- `functions/plottings/plot_acf_and_theoretical.m`

Main outputs in practice:

- top-window postprocessing figures (`main`, `timescale`, `deconv`, `deconv_localwin_norm`)
- residual-vs-MUA comparison MAT/CSV outputs and overview figures
- residual-vs-SPKT comparison MAT/CSV outputs and overview figures
- SPKT-residual lagged cross-correlation MAT/CSV outputs and summary figures

Combined manifests still exist in some runners, but they should be treated as
bookkeeping rather than the primary scientific outputs of this pipeline.

### 7. BOLD postprocessing pipeline

Primary scripts:

- `scripts/script_run_one_cfg_bold_reskoopnet_postprocessing.m`
- `scripts/script_postprocess_bold_reskoopnet_results.m`

Primary functions:

- `functions/postprocessing/build_bold_reskoopnet_postprocessing_params.m`
- `functions/postprocessing/discover_completed_bold_reskoopnet_runs.m`
- `functions/postprocessing/prepare_bold_reskoopnet_postprocessing_run.m`
- `functions/postprocessing/run_bold_reskoopnet_postprocessing_stages.m`
- `functions/postprocessing/postprocess_bold_reskoopnet_results.m`
- `functions/postprocessing/export_bold_intrinsic_activation_maps.m`
- `functions/postprocessing/build_bold_activation_plot_context.m`
- `functions/postprocessing/plot_bold_activation_map_reference_style.m`
- `functions/postprocessing/postprocess_EDMD_outputs.m`
- `functions/postprocessing/postprocess_EDMD_outputs_deconv_efuns.m`
- `functions/postprocessing/postprocess_EDMD_outputs_timescale.m`

Convenience wrappers:

- `scripts/script_run_one_bold_reskoopnet_postprocess.m`
- `functions/postprocessing/run_one_bold_reskoopnet_postprocessing_core.m`

### 8. BLP-BOLD cross-modal coupling pipeline

Primary scripts:

- `scripts/script_run_one_cfg_bold_cross_modal_coupling.m`
- `scripts/script_correlate_bold_efuns_with_densities.m`

Primary functions:

- `functions/postprocessing/build_bold_cross_modal_coupling_params.m`
- `functions/postprocessing/discover_completed_bold_cross_modal_coupling_runs.m`
- `functions/postprocessing/process_bold_cross_modal_coupling_runs.m`
- `functions/postprocessing/run_one_bold_cross_modal_coupling_core.m`
- `functions/postprocessing/compute_bold_efun_density_cross_correlation.m`
- `functions/postprocessing/export_bold_top_xcorr_activation_maps.m`
- `functions/plottings/plot_bold_efun_density_cross_correlation_summary.m`

Inputs usually come from:

- BLP event-state pipeline outputs such as event density or state diversity
- BOLD postprocessing pipeline outputs such as `BOLD_POST` MAT files

Compatibility wrappers still exist for the previous mixed entry points:

- `scripts/script_run_one_bold_cross_modal_coupling.m`
- `scripts/script_run_one_bold_reskoopnet_post_xcorr_activation_maps.m`
- `scripts/script_generate_missing_bold_activation_maps.m`

These wrappers now delegate to the separated pipeline 7 and pipeline 8
implementations, but they should not be treated as canonical pipeline-owned
entry points.

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
least five different downstream lines:

- BLP event-state downstream
- BLP eigenfunction dimension reduction
- BLP eigenfunction postprocessing
- BOLD postprocessing
- BLP-BOLD cross-modal coupling

That is why the folder feels structurally noisy even when many individual files
are internally reasonable.
