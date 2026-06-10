# Current Pipeline Inventory, Short Version

Last updated: 2026-05-22

This is a compact map of what the project currently contains. It is meant as
a readable orientation document, not a full mathematical specification.

The detailed reference remains:

```text
docs/project_mainline_pipeline_reference_2026-05-16.md
```

Status snapshot, 2026-05-22:

- Completed: invalid/stale BLP P4 registry is in place and used by P5/P6/P11 discovery.
- Completed: P11 resolves the current P7 mainline from the BOLD P4 best checkpoint, then applies the same provenance chain to P8/P9/P10.
- Completed: P2 event-diversity/window outputs are minor/nonblocking for default completeness.
- Completed: P8/P10 default runs are numeric-first; MAT/CSV xcorr outputs are the completeness source, while activation/ROI/browse figures are optional/backfill.
- Completed: P8 first-stage cross-session consistency outputs exist for top5 strength, method x k, preference, lag, and ranking tables.
- Updated: P11 should be run from the single wrapper `scripts/p11_current_mainline_check.py`; lower-level audit/scorecard scripts are implementation details.
- Updated: P11 wrapper now writes public P1-P10 outputs with stable names:
  `P1-P10_current_status.csv`, `P1-P10_recompute_plan.csv`,
  `P1-P10_flat_figure_manifest.csv`, `P1-P10_legacy_candidates.csv`, and
  `summary.md`.
- Updated: P4/P7 checks now treat `raw` and `standardize` as explicit parameter modes. For each dataset and current observable/BLP condition, P11 should report whether both versions exist and should select best checkpoints separately within each mode.
- Updated: P8/P10 interpretation tables must split BOLD `efun` and `deconv_efun` hits and report top raw LFP eigenfunction index/timescale plus top dimred LFP component process labels.
- Updated 2026-05-28: raw/observable/component trace browsing figures should
  no longer hard-clip display amplitude.  This is a plotting/QC provenance
  change, not a numeric pipeline change.  Old raw-trace browsing figures made
  before this policy can look saturated because of display clipping and should
  be treated as stale for visual-amplitude interpretation.  Details:
  `docs/plotting_no_trace_clipping_archive_2026-05-28.md`.
- Implemented in code: P5 thresholded raw/dimred density now defaults to
  `abs` activity magnitude and writes activity metadata, raw efun index,
  eigenvalue-derived timescale/frequency, and dimred component weighted
  timescale metadata. The preferred slow/fast-mode timescale is now the
  discrete-eigenvalue definition `-dt/log(abs(lambda_d))`, with frequency from
  `abs(angle(lambda_d))/(2*pi*dt)`; bilinear continuous-time values are kept as
  legacy/provenance columns only. P11 marks older P5 density MATs as stale
  against this new provenance rule.
- Still pending: legacy path cleanup, P6 MUA all-channel rerun, P5
  activity/envelope peak-statistics regeneration, P5 selective-component 3D
  trajectory exporter, P3 QC build-path integration, P9 figure-stage
  normalization, P8/P10 display simplification, and missing data backfill
  reported by P11 audit.

## Current Scope

Current downstream datasets:

```text
e10gb1
e10fV1
e10gh1
f12m01
```

`e10gW1` exists in the project but should be treated as not ready for the
current P7-P10 downstream set while its P4/BOLD training is still being
resolved.

Current downstream BOLD observables:

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

Other historical BOLD observable options exist, but they are not the current
minimal mainline unless explicitly re-enabled.

Current main processed-data root:

```text
E:\DataPons_processed\<dataset>\
```

Current P11 summary root:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
```

Rule of thumb: dataset-local pipeline folders are the source of truth.
`summary_figures` is a browsing/summary layer and should not be used to decide
whether a scientific result exists.

## Big Picture

The project has two input modalities:

- BLP/LFP: spectrograms, events, consensus states, learned eigenfunctions, and
  thresholded density sources.
- BOLD: ROI time series, BOLD observables, BOLD ResKoopNet outputs, BOLD
  eigenfunctions, and BOLD spatial profiles.

The later pipelines ask how learned BLP/LFP dynamics relate to learned BOLD
dynamics across sessions.

## Pipeline Summary

| Pipeline | Short name | Main question | Current source/output |
|---|---|---|---|
| P1 | BLP spectrogram/dictionary | Turn raw BLP into spectrograms and ResKoopNet dictionaries. | `pipeline1_spectrograms`, `pipeline1_reskoopnet_dictionary` |
| P2 | BLP events/states | Detect events, event density, consensus states, and diversity windows. | `pipeline2_event_detection`, `pipeline2_event_density`, `pipeline2_consensus_states`, diversity-window folders |
| P3 | BOLD observables | Build BOLD observable MATs and QC figures. | `pipeline3_bold_observables`, `pipeline3_figures_bold_pre_reskoopnet_qc` |
| P4 | ResKoopNet training/export | Train/export BLP and BOLD ResKoopNet/EDMD outputs. | external AutoDL/training output roots; downstream uses provenance |
| P5 | BLP eigenfunction reduction | Reduce BLP eigenfunctions and build density/event/component summaries. | `pipeline5_eigenfunction_reduction` plus P5 density/peak outputs |
| P6 | BLP top-window postprocess | Postprocess top state-diversity windows and residual channel coupling. | `pipeline6_top_state_diversity_postprocessing`, SPKT/MUA xcorr folders |
| P7 | BOLD postprocess | Postprocess current best BOLD checkpoints into BOLD_POST, timescales, deconv, ROI/activation summaries. | `pipeline7_bold_reskoopnet_postprocessing` |
| P8 | BOLD efun vs LFP density xcorr | Test whether BOLD efun/deconv streams couple to LFP density sources. | `pipeline8_xcorr`; optional `pipeline8_top_maps` |
| P9 | BOLD eigenfunction reduction | Reduce BOLD eigenfunction/deconv streams into BOLD components. | `pipeline9_bold_eigenfunction_reduction` |
| P10 | P9 components vs LFP density xcorr | Test whether reduced BOLD components couple to LFP density sources. | `pipeline10_dimred_xcorr`; optional `pipeline10_dimred_top_maps` |
| P11 | Current analysis organizer | Audit, flatten, summarize, and make current first-stage figures. | `summary_figures\pipeline11_current_analysis_summary`, `results\pipeline11_*` |

## P1: BLP Spectrogram And Dictionary

Purpose:

- Convert raw BLP/LFP into region spectrogram features.
- Build ResKoopNet-ready dictionary/observable packages.

Canonical folders:

```text
E:\DataPons_processed\<dataset>\pipeline1_spectrograms\
E:\DataPons_processed\<dataset>\pipeline1_reskoopnet_dictionary\
```

Typical figures:

- Raw BLP/spectrogram QC.
- Region spectrogram checks.
- These are useful for input sanity checks, but they are not currently the
  main analysis endpoint.

## P2: BLP Events, Density, Consensus States

Purpose:

- Detect band/event-like BLP events.
- Build event density in binned time.
- Build consensus states and diversity windows.

Canonical folders:

```text
E:\DataPons_processed\<dataset>\pipeline2_event_detection\
E:\DataPons_processed\<dataset>\pipeline2_event_density\
E:\DataPons_processed\<dataset>\pipeline2_consensus_states\
E:\DataPons_processed\<dataset>\pipeline2_consensus_state_diversity_windows\
E:\DataPons_processed\<dataset>\pipeline2_event_diversity_windows\
```

Current note:

- Implemented in P11: event-diversity/window products are useful but treated
  as minor for current completeness checks unless explicitly needed.
- P8/P10 can use P2 event density as one of the LFP density sources.

## P3: BOLD Observables

Purpose:

- Convert BOLD ROI time series into standardized observable packages for BOLD
  ResKoopNet.
- Current downstream observable set is four modes:

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

Canonical folders:

```text
E:\DataPons_processed\<dataset>\pipeline3_bold_observables\
E:\DataPons_processed\<dataset>\pipeline3_figures_bold_pre_reskoopnet_qc\
```

Current role in P11:

- Flatten current P3 QC figures into P11 summary.
- Do not use old copied QC folders as completeness proof.

## P4: ResKoopNet Training And Export

Purpose:

- Train BLP and BOLD ResKoopNet/MLP models.
- Export EDMD-like outputs that P5/P6/P7 consume.

Important point:

- P4 output lives outside the dataset-local `E:\DataPons_processed` pipeline
  folders.
- Downstream P7/P8/P9/P10 must resolve the current best BOLD checkpoint first,
  then treat outputs from older checkpoints as stale/legacy.

Current BOLD rule:

- Implemented in P11: for each dataset x observable, choose the best checkpoint
  by validation metric, not by newest folder timestamp.
- P11 must do this separately for `raw` and `standardize` P4 parameter modes.
  A single mixed best-checkpoint search is not acceptable because it hides
  whether a downstream P7/P8/P9/P10 result came from raw or time-standardized
  training.
- P7 is a deterministic postprocessing snapshot of that chosen checkpoint.
- If P4 best checkpoint changes later, P7/P8/P9/P10 should be re-audited or
  regenerated.

## P5: BLP Eigenfunction Reduction

Purpose:

- Reduce learned BLP eigenfunctions/components with multiple methods.
- Build density sources used later by P8/P10.
- Analyze component activity around consensus states/events.

Current main grid:

```text
condition: abs, csplit
method: svd, nmf, mds, umap
k: 03:08
```

Canonical reduction folder:

```text
E:\DataPons_processed\<dataset>\pipeline5_eigenfunction_reduction\<condition_tag>\<method_tag>\
```

This folder saves the actual dimred efun/component result, not only a figure
cache.  The canonical MAT under `mat\*_efun__*.mat` stores `result.core`
component time courses and mode weights plus `quality`, `summary`, and
provenance.  The default `compact` payload omits the largest raw source arrays,
but the saved reduction files can still be hundreds of MB.  P5 dimred density,
peak/selectivity outputs, trajectories, and P8/P10 density-source analyses
should trace back to these saved reduction MATs.

Important P5 derived folders:

```text
E:\DataPons_processed\<dataset>\pipeline5_raw_thresholded_density\
E:\DataPons_processed\<dataset>\pipeline5_dimred_thresholded_density\
E:\DataPons_processed\<dataset>\pipeline5_eigenfunction_peaks_by_state\
E:\DataPons_processed\<dataset>\pipeline5_eigenfunction_peaks_by_state_maxabs\
```

Current P5 notes:

- Scatter trajectory figures are the desired current visual style.
- Old line/surface trajectory figures should be treated as historical unless
  regenerated intentionally.
- Selective dimred-component 3D consensus-state trajectory plotting is planned
  but not urgent. The trajectory should not blindly use the first three
  dimensions; it should first choose the selective dimred components and then
  plot consensus-state trajectories in that selected 3D space.
- Planned P5 process-label table: generate `dimred_efun_process_labels.csv`
  for every dimred component by reusing the current theta/ripple/gamma
  selectivity logic on the new activity-magnitude peak statistics
  (`abs(component)` first, `rms_envelope` / eigenvalue-window-aware envelope
  later).
- Implemented density metadata: raw efun density is now an activity-magnitude
  product using `abs(phi)` with `lfp_activity_transform=abs_magnitude` and
  `lfp_activity_window_policy=samplewise_abs_no_envelope`. Optional
  `rms_envelope(phi)` remains a later refinement.
- Implemented density metadata: raw density MATs save `mode_metadata`
  including `raw_efun_index`, source eigenvalue, continuous-time real/imag
  part, preferred discrete-log timescale/frequency, bilinear legacy
  timescale/frequency, activity transform, and window policy.
- Implemented density metadata: dimred density MATs save
  `component_timescale_metadata` with component weighted mean/median
  timescale, weighted frequency, top source modes, top raw efun indices, and
  activity policy.
- Current `maxabs` selectivity summaries are prototype/reference outputs only;
  mainline P8/P10 labels should come from the activity/envelope branch.
- P5 summary/copy folders are caches, not source of truth.
- The current activity-density rerun keeps `top30_window_plots` off by
  default because the top-window branch can still reference historical
  consensus-window paths. Top-window plots remain optional/backfill browsing.

Useful current scripts:

```text
scripts/script_run_current_p5_activity_density_unblocked.m
scripts/check_current_pipeline5_completeness.py
scripts/script_run_current_pipeline5_peak_state_all_datasets.m
scripts/script_redraw_current_pipeline5_scatter_figures_from_mat.m
scripts/summarize_peak_state_method_consistency.py
scripts/summarize_peak_event_family_component_activity.py
```

## P6: BLP Top-Window Postprocessing

Purpose:

- Postprocess top state-diversity windows from BLP ResKoopNet outputs.
- Build run-level timescale diagnostics.
- Compare residual channels with SPKT/MUA.

Canonical folders:

```text
E:\DataPons_processed\<dataset>\pipeline6_top_state_diversity_postprocessing\
E:\DataPons_processed\<dataset>\pipeline6_spkt_residual_cross_correlation\
E:\DataPons_processed\<dataset>\pipeline6_mua_residual_cross_correlation\
E:\DataPons_processed\<dataset>\pipeline6_figures_timescale_diagnostics\
E:\DataPons_processed\<dataset>\pipeline6_figures_spkt_residual_cross_correlation\
E:\DataPons_processed\<dataset>\pipeline6_figures_mua_residual_cross_correlation\
```

Default P11 P6 flat views:

- Run-level timescale diagnostics.
- SPKT residual xcorr overview.
- MUA residual xcorr overview.
- Top-window EDMD postprocess paired with top-window window-normalized deconv.

Current note:

- `window_figures` should be treated as an old/copy cache, not canonical.
- Lagged SPKT branch is optional unless explicitly enabled.
- P6 residual-channel QC should compare residual modes with all available channels for both SPKT and MUA/MUA-proxy. Current code defaults SPKT to all channels but MUA to selected/observable channels; that is a planned fix, and selected-channel MUA outputs should be treated as stale for the all-channel QC.
- Known stale/corrupted figure to backfill: f12m01 `projected_vlambda_abs` run-level timescale PNG from `2026-05-14` is black at the canonical source, so P11 only copied a bad source image. Re-export this figure during the next P6 backfill, and treat it as stale until then.
- P11 should add a basic figure-quality check for required P6 run-level PNGs so blank/black renderer failures do not count as complete outputs.

## P7: BOLD ResKoopNet Postprocessing

Purpose:

- Convert current best BOLD ResKoopNet checkpoints into BOLD_POST artifacts.
- Produce BOLD eigenfunction, deconv, timescale, ROI, and activation summaries.

Canonical folder:

```text
E:\DataPons_processed\<dataset>\pipeline7_bold_reskoopnet_postprocessing\<run>\
```

Current mainline rule:

- P7 must be selected from the current P4 best checkpoint.
- P11 should not decide mainline P7 by scanning existing P7 folders alone.
- Old P7 folders from non-best checkpoints are legacy/stale for current
  completeness.

Figure/output families:

- Main BOLD EDMD/eigenfunction summary.
- Deconv / Koopman residual eigenfunction summary.
- Timescale diagnostics.
- Intrinsic ROI bar summary.
- Intrinsic activation maps.

Current first-stage preference:

- For broad method selection, P7 activation maps and ROI maps are secondary.
- They become important after narrowing candidates.

## P8: BOLD Efun Vs LFP Density Cross-Correlation

Purpose:

- Ask whether a BOLD observable's eigenfunction/deconv streams are stably
  coupled to LFP/BLP density sources across sessions.

Canonical numeric folder:

```text
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\
```

Optional spatial/browse folder:

```text
E:\DataPons_processed\<dataset>\pipeline8_top_maps\<run_tag>\
```

Current P8 density source grid:

```text
event_density
raw_abs_q070
raw_csplit_q070
dim_abs_<svd/nmf/mds/umap><k03:k08>_q070
dim_csplit_<svd/nmf/mds/umap><k03:k08>_q070
```

This is 51 density sources:

- 1 event density
- 2 raw densities
- 48 dimred densities

Current P8 numeric output levels:

- Combined xcorr top tables.
- Per-density top tables.
- Per-density-feature top tables.

P8 first-stage analysis should focus on:

1. Strength consistency:
   dataset x 51 density source, split by `all`, `efun`, `deconv_efun`.
2. Method x k consistency:
   dimred only, split by `abs/csplit` and `efun/deconv_efun`.
3. Preference consistency:
   per-session winner for all-density, dimred-only, abs-only, csplit-only.
4. Lag consistency:
   median lag and lag-sign agreement.
5. Top-hit interpretation:
   split BOLD `efun` and `deconv_efun`, then report top raw LFP efun
   index/timescale and top dimred LFP component process labels.

Current P8 first-stage outputs:

```text
results\pipeline8_cross_session_consistency_current\p8_main_consistency_ranking.csv
results\pipeline8_cross_session_consistency_current\p8_strength_cross_session_summary.csv
results\pipeline8_cross_session_consistency_current\p8_method_k_consistency_summary.csv
results\pipeline8_cross_session_consistency_current\p8_preference_winners_by_dataset.csv
results\pipeline8_cross_session_consistency_current\p8_preference_agreement_summary.csv
results\pipeline11_parameter_selection_scorecard_current\top_hit_interpretation_summary.csv
```

Current P8 first-stage figures:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\cross_session_consistency\strength\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\cross_session_consistency\method_k\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\cross_session_consistency\preference\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\cross_session_consistency\lag\
```

Current P8 scripts:

```text
scripts/summarize_pipeline8_cross_session_consistency.py
scripts/export_p8_p10_roi_profile_consistency_sources.m
scripts/plot_p8_p10_roi_profile_consistency.py
```

Current policy:

- Do not default to activation maps or ROI maps for first-stage P8 selection.
- Use activation/ROI only after a small number of stable candidates has been
  selected.
- Implemented first-stage consistency uses `topN = 5` and writes strength,
  method x k, preference, lag, and compact ranking outputs.
- Planned interpretation layer: join P5 activity-magnitude
  `dimred_efun_process_labels.csv` onto dimred-density P8 top hits so strong
  couplings can be read as theta-selective, ripple-selective, theta/ripple
  mixed, broad-active, or unlabeled LFP subprocess candidates.
- Planned raw-efun interpretation layer: for raw-density P8 top hits, preserve
  `raw_efun_index` and preferred discrete-log timescale metadata, then summarize
  whether strong xcorr is concentrated in fast, intermediate, or slow LFP
  eigenfunction modes. Default P11 figure should include top xcorr raw efun
  index distribution histograms and normalized index-quantile distributions
  using `topN = 5`, plus raw-vs-dimred mean/median timescale comparison
  heatmaps.

## P9: BOLD Eigenfunction Reduction

Purpose:

- Apply a P5-like dimension reduction to BOLD eigenfunction/deconv streams.
- Create BOLD component time series and component spatial loadings for P10.

Canonical folder:

```text
E:\DataPons_processed\<dataset>\pipeline9_bold_eigenfunction_reduction\<run_tag>\<feature>\<method_k>\
```

Current main feature:

```text
efun_real
```

Current method/k grid:

```text
method: svd, nmf, mds, umap
k: 03:08
```

Current note:

- P9 figures are co-located with MAT results under the same folder.
- A separate `pipeline9_figures_bold_eigenfunction_reduction` layout has been
  discussed, but current P11 should still read the co-located current outputs.
- P9 must trace back to current P7 BOLD_POST provenance.

## P10: P9 BOLD Components Vs LFP Density Cross-Correlation

Purpose:

- Compare reduced BOLD components from P9 with the same LFP density source grid
  used by P8.

Canonical numeric folder:

```text
E:\DataPons_processed\<dataset>\pipeline10_dimred_xcorr\<run_tag>__<p9_feature>__<method_k>\
```

Optional spatial/browse folder:

```text
E:\DataPons_processed\<dataset>\pipeline10_dimred_top_maps\<p10_tag>\
```

Current grid:

```text
run_tag: pv_gsvd100, pv_gsvd100_ds, pv_hp100, pv_roi
p9_feature: efun_real
p9_method_k: svd/nmf/mds/umap x k03:k08
lfp_density_sources: same 51-source P8/P10 grid
```

Current P10 cross-session outputs:

```text
results\pipeline10_cross_session_consistency_current\p10_top_xcorr_hits_readable.csv
results\pipeline10_cross_session_consistency_current\p10_xcorr_density_feature_scores.csv
results\pipeline10_cross_session_consistency_current\p10_cross_session_density_feature_summary.csv
```

Current P10 figures:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p10\cross_session_consistency\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p10\roi_profile_consistency\
```

Current P10 scripts:

```text
scripts/summarize_pipeline10_cross_session_consistency.py
scripts/export_p8_p10_roi_profile_consistency_sources.m
scripts/plot_p8_p10_roi_profile_consistency.py
```

Current policy:

- P10 inherits the P7 -> P9 provenance chain.
- Existing P10 figure/cache folders should not be used as completeness proof.
- Like P8, P10 should first be analyzed by numeric consistency before drawing
  activation/ROI maps.
- P10 should mirror P8's top-hit interpretation: split by BOLD component-side
  feature family, report top raw LFP efun index/timescale, and join dimred LFP
  process labels for each top dimred density hit.
- Implemented numeric-first/default provenance checks exist, but P10 is still
  less complete than P8 because source-component interpretability,
  theta/ripple-specific selectivity, and some data backfill remain pending.
- Planned interpretation layer: use the same P5 activity-magnitude dimred
  component labels as P8 for dimred LFP density sources, while keeping this
  separate from P10 source-component interpretability for the BOLD dimred side.
- Planned raw-efun interpretation layer: mirror P8 by adding raw LFP efun index
  and timescale metadata to raw-density P10 top hits, then summarize whether
  stable BOLD dimred couplings prefer fast raw LFP eigenfunction modes. Include
  the same top xcorr raw efun index distribution views as P8, plus the matched
  raw-vs-dimred component timescale comparison.

## P11: Current Analysis Organizer

Purpose:

- Audit current outputs.
- Flatten current figures into a single browsing folder.
- Regenerate derived summaries.
- Keep historical/cache/legacy outputs out of the current view unless requested.

Main P11-like scripts/results:

```text
scripts/p11_current_analysis_audit.py
scripts/p11_current_mainline_check.py
scripts/script_backfill_current5_minimal_pipeline_set.m
scripts/audit_p1_p4_readiness.py
scripts/build_pipeline5_dimred_component_process_labels.py
scripts/summarize_p11_parameter_selection_scorecard.py
results\pipeline11_current_analysis_audit\
results\pipeline11_p1_p4_readiness_current\
results\pipeline5_dimred_component_process_labels_current\
results\pipeline11_parameter_selection_scorecard_current\
results\pipeline11_current_mainline_check\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\
```

Current P11 should treat these as derived current-analysis products:

- Public P1-P10 wrapper outputs:
  - `results\pipeline11_current_mainline_check\P1-P10_current_status.csv`
  - `results\pipeline11_current_mainline_check\P1-P10_recompute_plan.csv`
  - `results\pipeline11_current_mainline_check\P1-P10_flat_figure_manifest.csv`
  - `results\pipeline11_current_mainline_check\P1-P10_legacy_candidates.csv`
  - `results\pipeline11_current_mainline_check\summary.md`
- P3 QC flat views.
- P5 peak statistics, method consistency, event-family activity/selectivity.
- P6 diagnostic flat views.
- P7 current-mainline flat views and possible intrinsic ROI consistency.
- P8 first-stage consistency views.
- P9 current BOLD-dimred summaries.
- P10 cross-session and ROI-profile consistency views.
- Parameter-selection scorecard for choosing:
  `abs/csplit`, dimred method, BOLD observable type, and dimred component count.
- Initial scorecard implementation exists repo-locally and P11 audit supports
  `--dataset-scope all9` for all nine cfg datasets.
- Initial P1-P4 readiness audit exists repo-locally and is now registered in
  P11 audit, so front-end gaps are not hidden by downstream P5-P10 checks.

Important P11 policy:

- Daily/default entry point:

```text
python scripts\p11_current_mainline_check.py --dataset-scope all9
```

- Run headless by default.
- Do not pop up MATLAB figures during backfill.
- Do not delete legacy files automatically.
- Move questionable historical outputs to legacy/quarantine first.
- Completeness should be decided from canonical pipeline folders, not summary
  figure caches.

## Current First-Stage Analysis Philosophy

For deciding what to trust first, the current project should not start from
activation maps or ROI maps.

Recommended first-stage order:

1. Confirm data completeness and current-best provenance.
2. Compare numeric cross-session strength.
3. Check method/k stability.
4. Check winner/preference consistency.
5. Check lag sign/range.
6. Check P5 theta/ripple process-label agreement and raw-vs-dimred timescale
   agreement.
7. Build the parameter-selection scorecard.
8. Only then inspect ROI maps and activation maps for a small candidate set.

This applies especially to P8 and P10.

The scorecard should combine top-N xcorr strength, dataset coverage, label
agreement, lag consistency, ROI/profile consistency when available,
raw-efun-index/timescale concentration, raw-vs-dimred timescale agreement, and a
complexity preference for simpler methods or smaller `k` when the scientific
scores are similar.

Default scorecard figures:

- Metric matrix for top parameter candidates.
- `abs` vs `csplit` condition comparison.
- Method x k heatmaps.
- BOLD observable x metric profile heatmap.
- k sensitivity curves.
- Process-label composition stacked bars.
- Raw-vs-dimred timescale agreement plot.
- Compact candidate decision panel.

## Known Confusing/Legacy Areas

Historical folders still exist. They should not be deleted until checked, but
they should not define the current result.

Examples:

- Old P5 trajectory styles: line/surface instead of current scatter.
- P5 `pipeline5_summary_figures`: useful cache, not canonical.
- P6 `window_figures`: copied browsing cache.
- P8 old folders under some datasets, such as old cross-correlation or flat
  top-map folders.
- P10 old summary/cache folders.
- Any non-P11 `summary_figures\pipeline8*` or `summary_figures\pipeline10*`
  folder may be stale.

## Minimal Information Needed For A New Analysis Design

To design a clean analysis rather than another large pile of figures, define:

```text
1. Selection target:
   P5 method? P8 density source? P10 BOLD dimred method? observable?

2. Ranking unit:
   e.g. observable x feature_family x condition x method x k

3. Inclusion/exclusion:
   raw density allowed?
   event_density allowed?
   HP_svd100 included?
   efun only or deconv_efun too?

4. Success criteria:
   n_datasets = 4/4
   high mean score
   low std/range
   winner method/k agreement
   lag sign agreement
   ROI profile consistency

5. Output size:
   one ranking table?
   three figures?
   per-dataset panels?
```

Without this contract, the project naturally produces too many plausible but
hard-to-interpret summary figures.
