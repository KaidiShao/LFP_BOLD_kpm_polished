# Current Storage Usage and Migration Plan

Date: 2026-06-05

This document records a read-only storage audit for the current
LFP/BOLD Koopman project and proposes a non-destructive migration plan.
Nothing has been moved or deleted during this audit.

## 1. Current Drive Status

| Drive | Total GB | Used GB | Free GB | Notes |
|---|---:|---:|---:|---|
| E:\ | 7451.91 | 6864.77 | 587.14 | Main active project/result drive. Too full for many more large reruns. |
| G:\ | 3725.99 | 3068.61 | 657.38 | Current external/raw-data drive. Also nearly full. |

Main conclusion:

- E:\ is no longer safe as both an active compute drive and a historical archive.
- G:\ can hold raw data or temporary new input, but it does not have enough free space to absorb the current E:\ legacy/archive burden.
- A third large archive drive is strongly recommended before deleting anything.

## 2. Project-Related Roots on E:

| Root | Size GB | Files | Dirs | Current interpretation |
|---|---:|---:|---:|---|
| E:\DataPons_processed | 2219.84 | 483637 | 34278 | Main processed data/results root. Contains both current and legacy outputs. |
| E:\DataPons | 1362.39 | 14742 | 495 | Raw/input data root. Can live on an external drive if paths are configured cleanly. |
| E:\autodl_results | 972.55 | 41394 | 328 | Older/default training-output root. Strong archive candidate if current code uses autodl_results_new. |
| E:\autodl_results_new | 491.40 | 12650 | 288 | Current BLP P4 training-output root. Keep active for now. |
| E:\autodl_results_local | 63.21 | 42464 | 1722 | BOLD/local output root. Review before moving. |

These five roots alone occupy about 5.11 TB.

## 3. Current Minimal Mainline Definition

The current first-stage analysis should be kept narrow:

```text
BLP condition:
    complex_split_projected_vlambda_standardize

BLP activity transform:
    adaptive_envelope primary
    abs as QC/control

P5 label:
    P2 theta/gamma/ripple band-event strict label only

BOLD side:
    first stage: roi_mean only
    second stage: global_svd100 / gsvd100_ds / HP_svd100

P8/P10:
    efun and deconv_efun analyzed separately
    first-stage question: whether roi_mean deconv_efun is more RG-like than roi_mean efun
```

This means the storage plan should preserve the current standardized complex-split
P5/P8/P10 route and avoid keeping every historical density/event/trajectory
variant on the active E:\ drive.

## 4. Main Processed Root: E:\DataPons_processed

Top-level size distribution:

| Directory | Size GB | Status | Notes |
|---|---:|---|---|
| derived_autodl_results_standardize | 412.15 | REVIEW_DEPENDENCY | Very large. Likely standardized P4/P5 export/provenance source. Check before moving. |
| e10gb1 | 268.36 | KEEP_CURRENT | Current/important processed dataset. |
| e10fV1 | 206.54 | KEEP_CURRENT | Current/important processed dataset. |
| f12m01 | 192.17 | KEEP_CURRENT | Current/important processed dataset. |
| legacy_pipeline4plus_20260511 | 184.84 | ARCHIVE_CANDIDATE | Historical. Move to archive drive first, do not delete immediately. |
| e10gh1 | 184.59 | KEEP_CURRENT | Current/important processed dataset. |
| e10gw1 | 137.17 | KEEP_CURRENT | Current/important processed dataset. |
| _legacy_quarantine | 92.95 | ARCHIVE_CANDIDATE | Historical quarantine. Move to archive drive first. |
| e10aw1 | 83.80 | REVIEW_DATASET_SCOPE | Has P1-P3 only; not in current minimal P5/P8/P10 story unless intentionally included. |
| e10bv1 | 83.33 | REVIEW_DATASET_SCOPE | Has P1-P3 only; not in current minimal P5/P8/P10 story unless intentionally included. |
| k13m17 | 79.50 | KEEP_CURRENT / QC_CAUTION | Current processed dataset, but interpretation has known caution. |
| k13m19 | 71.40 | KEEP_CURRENT | Current K13 processed dataset. |
| k13m18 | 67.26 | KEEP_CURRENT | Current K13 processed dataset. |
| k13m20 | 67.14 | KEEP_CURRENT | Current K13 processed dataset. |
| k13m23 | 50.63 | KEEP_CURRENT | Current K13 processed dataset. |
| k13m21 | 34.44 | KEEP_CURRENT | Current K13 processed dataset. |
| summary_figures | 2.80 | KEEP_QC | Small enough to keep active. |
| window_figures | 0.20 | ARCHIVE_CANDIDATE | Historical/cache-like summary. |
| derived_autodl_results_standardize_smoke | 0.20 | ARCHIVE_CANDIDATE | Smoke/test output. |
| derived_autodl_results_bold_full_export | 0.18 | REVIEW_DEPENDENCY | BOLD export/provenance. Small, keep until confirmed. |

Immediate low-risk archive candidates inside DataPons_processed:

```text
E:\DataPons_processed\legacy_pipeline4plus_20260511
E:\DataPons_processed\_legacy_quarantine
E:\DataPons_processed\window_figures
E:\DataPons_processed\derived_autodl_results_standardize_smoke
```

These would free roughly 278 GB from E:\ if moved elsewhere.

## 5. Pipeline-Level Usage in Current Dataset Folders

Across the currently relevant processed dataset folders
`e10gb1, e10fV1, e10gh1, e10gw1, f12m01, k13m17, k13m18, k13m19,
k13m20, k13m21, k13m23`:

| Pipeline directory | Total GB | Datasets | Status | Notes |
|---|---:|---:|---|---|
| pipeline5_eigenfunction_reduction | 385.28 | 11 | KEEP_CURRENT | Core P5 reduction MAT/results. Do not move unless path layer supports it. |
| pipeline1_spectrograms | 323.77 | 11 | KEEP_CURRENT / COLD_TIER_LATER | Needed for provenance/rebuilds. Can move to slower storage only after downstream is frozen. |
| pipeline1_reskoopnet_dictionary | 186.32 | 11 | KEEP_CURRENT / COLD_TIER_LATER | Needed for P4/efun provenance. |
| pipeline8_xcorr | 108.74 | 11 | KEEP_CURRENT | Main P8 coupling evidence. |
| pipeline10_dimred_xcorr | 104.03 | 11 | KEEP_CURRENT | Main P10 coupling evidence if P10 remains in the story. |
| pipeline5_dimred_thresholded_density | 84.47 | 11 | REVIEW_DEPENDENCY | May be read by P8/P10. Check before archiving. |
| pipeline10_dimred_top_maps | 66.79 | 4 | REGENERATABLE / REVIEW | Large maps/figures. Likely not first-stage minimal evidence. |
| pipeline6_spkt_residual_cross_correlation | 21.07 | 5 | SECONDARY | Keep if P6 remains active; not first-stage minimal story. |
| pipeline6_mua_residual_cross_correlation | 20.25 | 5 | SECONDARY | Keep if P6 remains active; not first-stage minimal story. |
| pipeline9_bold_eigenfunction_reduction | 15.11 | 11 | KEEP_CURRENT | Needed for P10 provenance. |
| pipeline8_top_maps | 13.43 | 4 | REGENERATABLE / REVIEW | Large maps/figures. Likely second-stage. |
| pipeline7_bold_reskoopnet_postprocessing | 6.52 | 11 | KEEP_CURRENT | Needed for BOLD efun/deconv and ROI provenance. |
| pipeline3_bold_observables | 6.36 | 11 | KEEP_CURRENT | BOLD observable input. |
| pipeline5_raw_thresholded_density | 5.96 | 11 | REVIEW_DEPENDENCY | Needed if raw density xcorr is active. |
| pipeline5_efun_dimred_top30 | 3.08 | 3 | REGENERATABLE / QC | Useful QC, but not the core numeric source. |
| pipeline2_event_detection | 1.55 | 11 | KEEP_CURRENT | Small and essential for P5 band-event labels. |
| pipeline5_eigenfunction_peaks_by_state_rmsenv_adaptive | 0.79 | 11 | KEEP_QC | Current adaptive-envelope P5 peak results. |
| pipeline5_event_diagnostic_plates | 0.24 | 11 | KEEP_QC | Current P5 explanation figures. |
| pipeline5_paired_subprocess_traces | 0.24 | 11 | KEEP_QC | Current P5 explanation figures. |
| pipeline5_component_activity_top30_windows | 0.11 | 11 | KEEP_QC | Current P5 top-window QC. |

Important interpretation:

- The current mainline is large even without legacy: P1/P5/P8/P10 are the real heavy parts.
- However, the most removable active-drive burden is still historical output, especially old `autodl_results` and legacy/quarantine folders.
- `pipeline10_dimred_top_maps` and `pipeline8_top_maps` are big enough to matter, but they should be treated as regeneratable/review rather than immediate deletion targets.

## 6. Raw/Input Data: E:\DataPons

Top raw/input folders:

| Directory | Size GB | Status | Notes |
|---|---:|---|---|
| K13.m23 | 189.39 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Very large raw/input dataset. |
| E10.bv1 | 110.04 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Not first-stage current unless intentionally included. |
| E10.aW1 | 107.64 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Not first-stage current unless intentionally included. |
| other | 93.86 | REVIEW_RAW | Mixed source; inspect before moving/deleting. |
| E10.gW1 | 85.14 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current/important input. |
| E10.gH1 | 81.51 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current/important input. |
| E10.fV1 | 75.37 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current/important input. |
| E10.gb1 | 69.84 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current/important input. |
| F12.m04 | 67.02 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Future/optional. |
| F12.m05 | 64.83 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Future/optional. |
| K13.m17 | 61.98 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current but QC caution. |
| K13.m19 | 59.15 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current K13 input. |
| K13.m20 | 58.51 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current K13 input. |
| K13.m18 | 49.57 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current K13 input. |
| F12.m01 | 48.44 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current/important input. |
| K13.m21 | 25.54 | RAW_INPUT_CAN_LIVE_ON_G_OR_H | Current K13 input. |

Raw data can be stored on an external drive, because it is mostly read-only.
The risk is speed and path stability, not scientific correctness.

Recommended rule:

```text
Raw input:
    can live on G:\ or new H:\

Active processed output:
    keep on E:\ unless the pipeline explicitly supports a different processed root

Archive / legacy:
    move to new H:\ archive drive
```

## 7. Training Output Roots

### E:\autodl_results

| Directory | Size GB | Status |
|---|---:|---|
| e10gb1 | 466.97 | ARCHIVE_CANDIDATE if superseded by autodl_results_new |
| e10fV1 | 266.75 | ARCHIVE_CANDIDATE if superseded by autodl_results_new |
| f12m01 | 128.08 | ARCHIVE_CANDIDATE if superseded by autodl_results_new |
| e10gh1 | 110.61 | ARCHIVE_CANDIDATE if superseded by autodl_results_new |

This root is the single largest obvious candidate for migration.
If all current BLP P4 scripts use `E:\autodl_results_new`, then moving
`E:\autodl_results` to an archive drive would free about 973 GB.

### E:\autodl_results_new

| Directory | Size GB | Status |
|---|---:|---|
| e10fV1 | 127.53 | KEEP_CURRENT |
| e10gw1 | 97.16 | KEEP_CURRENT |
| e10gb1 | 92.31 | KEEP_CURRENT |
| f12m01 | 88.46 | KEEP_CURRENT |
| e10gh1 | 82.30 | KEEP_CURRENT |
| K13/current small entries | < 1 each | KEEP_CURRENT / VERIFY_COMPLETENESS |

This is the current intended BLP training-output root and should remain active
until all downstream exports are confirmed.

### E:\autodl_results_local

Status: `REVIEW_DEPENDENCY`.

This may contain BOLD/local outputs. Do not move until P7/P8/P10 path discovery
is checked.

## 8. Storage Classes

Use these labels for future cleanup:

```text
KEEP_CURRENT
    Required by the current first-stage mainline.

KEEP_QC
    Small or high-value QC/summary figures that make interpretation easier.

REVIEW_DEPENDENCY
    Potentially important dependency/provenance source. Inspect scripts before moving.

REGENERATABLE
    Derived figures/maps/cache that can be recreated from current numeric outputs.

ARCHIVE_CANDIDATE
    Not current mainline. Move to archive drive first; do not delete immediately.

RAW_INPUT_CAN_LIVE_ON_G_OR_H
    Raw/source data. Can live on external storage if paths are configured.

SECONDARY
    Scientifically useful but outside the current minimal P5/P8/P10 story.
```

## 9. Recommended Migration Plan

### Phase 0: Do not delete anything yet

Before moving large folders:

1. Record source path, destination path, size, file count, and move date.
2. Copy or move to an archive drive.
3. Verify size and file count after transfer.
4. Keep the archive mounted until the next audit confirms no active path depends on it.

### Phase 1: Add a new archive drive

Recommended:

```text
12 TB or 16 TB external HDD
```

Use it as:

```text
H:\LFP_BOLD_archive\
```

Suggested structure:

```text
H:\LFP_BOLD_archive\2026-06-05_legacy_from_E\
    autodl_results\
    DataPons_processed_legacy_pipeline4plus_20260511\
    DataPons_processed_legacy_quarantine\
    DataPons_processed_window_figures\
```

Do not use the archive HDD as the main compute output root. It is for cold
storage, not active analysis.

### Phase 2: Move low-risk archive candidates

First migration batch:

```text
E:\autodl_results
E:\DataPons_processed\legacy_pipeline4plus_20260511
E:\DataPons_processed\_legacy_quarantine
E:\DataPons_processed\window_figures
E:\DataPons_processed\derived_autodl_results_standardize_smoke
```

Estimated E:\ space recovered:

```text
about 1.25 TB
```

This would raise E:\ free space from about 0.59 TB to about 1.84 TB.

### Phase 3: Review large dependency candidates

Inspect before moving:

```text
E:\DataPons_processed\derived_autodl_results_standardize
E:\DataPons_processed\<dataset>\pipeline5_dimred_thresholded_density
E:\DataPons_processed\<dataset>\pipeline5_raw_thresholded_density
E:\DataPons_processed\<dataset>\pipeline8_top_maps
E:\DataPons_processed\<dataset>\pipeline10_dimred_top_maps
E:\autodl_results_local
```

Questions to answer:

1. Does any current P5/P8/P10 script read this folder directly?
2. Is it a numeric source, or only a figure/cache output?
3. Can it be regenerated from kept MAT/CSV files?
4. Does moving it require path config changes or only archive bookkeeping?

### Phase 4: Raw data relocation

Raw data can live on G:\ or H:\, but G:\ is already nearly full.

If a new archive drive is available:

```text
H:\DataPons_raw\
```

can hold older or less active raw inputs, while G:\ can hold incoming/new data.

Recommended split:

```text
G:\DataPons_new_or_active_raw\
H:\DataPons_raw_archive\
E:\DataPons_processed_active\
```

### Phase 5: Keep E:\ as the active working drive

After migration:

```text
E:\DataPons_processed
E:\autodl_results_new
E:\autodl_results_local, if still active
```

should be the only large active roots on E:\.

The goal is to keep at least 1.5-2.0 TB free on E:\ before starting large
multi-dataset reruns.

## 10. What Not To Move Yet

Do not move these until the current P5/P8/P10 analysis is stable:

```text
E:\DataPons_processed\<current_dataset>\pipeline5_eigenfunction_reduction
E:\DataPons_processed\<current_dataset>\pipeline8_xcorr
E:\DataPons_processed\<current_dataset>\pipeline10_dimred_xcorr
E:\DataPons_processed\<current_dataset>\pipeline7_bold_reskoopnet_postprocessing
E:\DataPons_processed\<current_dataset>\pipeline9_bold_eigenfunction_reduction
E:\autodl_results_new
```

Also keep:

```text
E:\DataPons_processed\<dataset>\pipeline2_event_detection
E:\DataPons_processed\<dataset>\pipeline5_event_diagnostic_plates
E:\DataPons_processed\<dataset>\pipeline5_paired_subprocess_traces
E:\DataPons_processed\<dataset>\pipeline5_component_activity_top30_windows
E:\DataPons_processed\summary_figures
```

These are small relative to the whole drive and directly support interpretation.

## 11. Practical Recommendation

Buy one more large drive.

Preferred role:

```text
H:\ = archive / legacy / cold raw backup
```

Keep:

```text
E:\ = active processed outputs and current training outputs
G:\ = active/new raw input, if it stays connected
H:\ = archive and raw cold storage
```

This gives a safer storage hierarchy:

```text
fast/active:
    E:\

connected raw input:
    G:\

cold archive:
    H:\
```

The key principle is: move historical mass off E:\, but keep the current numeric
mainline on E:\ until the pipeline has explicit path support for external roots.

