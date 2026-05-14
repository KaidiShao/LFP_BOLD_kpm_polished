# Remaining Non-Pipeline Processed Directories

Date: `2026-04-30`

This note records the remaining processed-data directories under
`E:\DataPons_processed` that do **not** follow the canonical `pipeline*`
naming rule yet.

The goal is to preserve a short checklist for later cleanup.

## Root-Level Shared Directories

These are still intentionally kept at the processed-data root and are **not**
old per-pipeline stage directories:

- `E:\DataPons_processed\postprocessing_manifests`
- `E:\DataPons_processed\summary_figures`
- `E:\DataPons_processed\window_figures`

Interpretation:

- `postprocessing_manifests` is a shared manifest root used by multiple
  pipelines.
- `summary_figures` is a shared cross-dataset summary figure root.
- `window_figures` is a shared cross-pipeline window comparison root.

These should not be treated as leftover dataset-local pipeline stage folders.

## Dataset-Level Remaining Non-Pipeline Directories

### `e10fV1`

Remaining non-`pipeline*` directories:

- `E:\DataPons_processed\e10fV1\postprocessing`

Current interpretation:

- `postprocessing`
  - old `pipeline 5` residual directory
  - mainly old `eigenfunction_density_scan` history

### `e10gb1`

Remaining non-`pipeline*` directories:

- `E:\DataPons_processed\e10gb1\postprocessing`
- `E:\DataPons_processed\e10gb1\koopman_postprocessing`
- `E:\DataPons_processed\e10gb1\eigenfunction_reduction`

Current interpretation:

- `postprocessing`
  - old `pipeline 5` density-scan history
- `koopman_postprocessing`
  - old `pipeline 5` reduction/top30/peaks tree
- `eigenfunction_reduction`
  - old `pipeline 5` naming residue

### `e10gh1`

Remaining non-`pipeline*` directories:

- `E:\DataPons_processed\e10gh1\postprocessing`

Current interpretation:

- `postprocessing`
  - old `pipeline 5` density-scan history

### `f12m01`

Remaining non-`pipeline*` directories:

- `E:\DataPons_processed\f12m01\postprocessing`

Current interpretation:

- `postprocessing`
  - old `pipeline 5` density-scan history

## Recommended Follow-Up

When returning to this cleanup, review these in order:

1. `E:\DataPons_processed\e10gb1\eigenfunction_reduction`
   - verify whether it is fully redundant with:
     - `pipeline5_eigenfunction_reduction`
     - `pipeline5_efun_dimred_top30`
     - `pipeline5_eigenfunction_peaks_by_state`

2. `E:\DataPons_processed\e10gb1\koopman_postprocessing`
   - likely old `pipeline 5`
   - verify redundancy before deleting

3. `E:\DataPons_processed\<dataset>\postprocessing`
   - currently mainly old `pipeline 5` density-scan history
   - note that a copy of this history has already been preserved under:
     - `pipeline5_raw_thresholded_density_scan`

## Important Context

The old `pipeline 5` density-scan history was already copied into canonical
archive-style directories:

- `E:\DataPons_processed\e10fV1\pipeline5_raw_thresholded_density_scan`
- `E:\DataPons_processed\e10gb1\pipeline5_raw_thresholded_density_scan`
- `E:\DataPons_processed\e10gh1\pipeline5_raw_thresholded_density_scan`
- `E:\DataPons_processed\f12m01\pipeline5_raw_thresholded_density_scan`

So the later decision is mainly whether the old `postprocessing` copies can be
deleted after one more check.

## Update

The old `E:\DataPons_processed\e10fV1\bold_reskoopnet_single_run_analysis`
directory was later migrated into canonical `pipeline 8` directories:

- `E:\DataPons_processed\e10fV1\pipeline8_bold_efun_density_cross_correlation`
- `E:\DataPons_processed\e10fV1\pipeline8_figures_bold_top_xcorr_activation_maps`
