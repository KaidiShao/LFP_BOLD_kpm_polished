# Pipeline 5 Cleanup Follow-up

Date: `2026-04-29`

This note records the remaining manual cleanup decisions for `pipeline 5`.

## Current State

The canonical `pipeline 5` outputs now live under:

- `E:\DataPons_processed\<file_stem>\pipeline5_eigenfunction_reduction`
- `E:\DataPons_processed\<file_stem>\pipeline5_efun_dimred_top30`
- `E:\DataPons_processed\<file_stem>\pipeline5_raw_thresholded_density`
- `E:\DataPons_processed\<file_stem>\pipeline5_raw_thresholded_events`
- `E:\DataPons_processed\<file_stem>\pipeline5_dimred_thresholded_density`
- `E:\DataPons_processed\<file_stem>\pipeline5_dimred_thresholded_events`
- `E:\DataPons_processed\<file_stem>\pipeline5_eigenfunction_peaks_by_state`

Old threshold-scan history has also been copied into:

- `E:\DataPons_processed\e10fV1\pipeline5_raw_thresholded_density_scan`
- `E:\DataPons_processed\e10gb1\pipeline5_raw_thresholded_density_scan`
- `E:\DataPons_processed\e10gh1\pipeline5_raw_thresholded_density_scan`
- `E:\DataPons_processed\f12m01\pipeline5_raw_thresholded_density_scan`

The following non-canonical leftovers were already removed:

- `E:\DataPons_processed\e10fV1\postprocessing\logs`
- `E:\DataPons_processed\e10gh1\postprocessing\logs`
- `E:\DataPons_processed\f12m01\postprocessing\logs`
- `E:\DataPons_processed\e10gh1\postprocessing\e10gh1_raw_sampling_audit.csv`

## Pending Decisions

### 1. Old density-scan directories

These old directories were intentionally kept for now:

- `E:\DataPons_processed\e10fV1\postprocessing\eigenfunction_density_scan`
- `E:\DataPons_processed\e10gb1\postprocessing\eigenfunction_density_scan`
- `E:\DataPons_processed\e10gh1\postprocessing\eigenfunction_density_scan`
- `E:\DataPons_processed\f12m01\postprocessing\eigenfunction_density_scan`

Reason:

- They contain historical multi-threshold sweep outputs.
- A copy now exists in `pipeline5_raw_thresholded_density_scan`.

Later task:

- Decide whether the old `postprocessing\eigenfunction_density_scan` trees can be deleted now that the `pipeline5_raw_thresholded_density_scan` copies exist.

### 2. Old e10gb1 koopman postprocessing tree

This old tree was also intentionally kept for now:

- `E:\DataPons_processed\e10gb1\koopman_postprocessing\eigenfunction_reduction`

Reason:

- It appears to overlap with the current canonical `pipeline5_*` outputs, but it was not yet deleted until the overlap is reviewed one more time.

The likely modern counterparts are:

- `E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_reduction`
- `E:\DataPons_processed\e10gb1\pipeline5_efun_dimred_top30`
- `E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_peaks_by_state`

Later task:

- Recheck whether the old `koopman_postprocessing\eigenfunction_reduction` content is fully redundant with the current `pipeline5_*` directories.
- If yes, delete the old tree.

## Recommended Next Step

When returning to this cleanup:

1. Compare old `e10gb1\koopman_postprocessing\eigenfunction_reduction` files against:
   - `pipeline5_eigenfunction_reduction`
   - `pipeline5_efun_dimred_top30`
   - `pipeline5_eigenfunction_peaks_by_state`
2. If that old tree is confirmed redundant, delete it.
3. Then delete the old `postprocessing\eigenfunction_density_scan` trees, since the `pipeline5_raw_thresholded_density_scan` copies already preserve the history.
