# E10gb1 Zscore Discussion, Analysis, and Current Decision

Date: 2026-04-17

## Purpose

This note records the current discussion and decisions around whether
`zscore` should be applied in the `E10gb1` pipeline, especially for:

- event detection
- observable generation

The goal is to keep one stable written reference instead of relying on
 chat history.

## Short Conclusion

- For `event detection`, keep **global zscore per channel on the
  concatenated trace**.
- Do **not** switch `event detection` to `session-wise zscore` at this
  point.
- For `observable generation`, make **no change for now**.
- Current reasoning:
  - the sessions use the same electrodes over only a few hours
  - per-session mean drift is very small
  - per-session std drift exists, but is not a clean uniform long-term
    drift across all channels
  - session-wise zscore changes detected events materially and may remove
    real across-session activity differences

## What Was Investigated

### 1. Whether the old pipeline used extra normalization

The old prepared-MAT source for `test23` was traced back to the script:

- [sr_pl_e10gb1_koopman_dataprep_new.m](/D:/Onedrive/ICPBR/Alberta/koopman_events/code/data_preparation/sr_pl_e10gb1_koopman_dataprep_new.m:1)

Important finding:

- that preparation script does **not** apply zscore
- it simply concatenates `blp.dat(:, [sr_ch, pl_ch], 1)` across sessions

So the old prepared MAT was not normalized in that script.

### 2. Why the old MAT looked normalized while the current root did not

Different raw-data roots were found for the same dataset:

- old normalized-scale copy:
  - `D:\Onedrive\experimental_event_data\ripple_monkeys\E10.gb1\blp\`
  - `E:\DataPons\E10.gb1.old\blp\`
- current repo default root:
  - `E:\DataPons\E10.gb1\blp\`

Key finding:

- the old root and old prepared MAT match exactly
- the current repo default root is a different file version

This explains why the old prepared MAT had values around `std ~ 1`, while
the current root had much larger raw-channel scales.

### 3. Whether global zscore changes event detection on the current root

This was tested directly on the current `cfg_E10gb1()` raw-data root by
comparing:

- `input_normalization = 'none'`
- `input_normalization = 'zscore_per_channel'`

Result:

- they produced **exactly the same detected peaks**
- total peaks:
  - `none = 26742`
  - `global zscore = 26742`
- per band / per channel overlap:
  - `Jaccard = 1.0` everywhere

Interpretation:

- for the current detector
  - `bandpass -> mean + k*std threshold`
- applying one fixed linear scaling per channel does not change the
  event-detection result on the current root

This is consistent with the detector form and was also confirmed
empirically.

### 4. Whether session-wise zscore changes event detection

This was tested using:

- [script_compare_e10gb1_zscore_scopes.m](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/tests/test24_detection_alignment_checks/script_compare_e10gb1_zscore_scopes.m:1)

Comparison:

- `global_zscore`
- `sessionwise_zscore`

Aggregate result:

- `global total peaks = 26742`
- `sessionwise total peaks = 26465`
- `shared peaks = 24776`
- `recall_vs_global = 0.92648`
- `precision_vs_global = 0.93618`
- `mean per-channel/band Jaccard = 0.86456`

Largest mismatches included:

- `gamma / ch4`
- `ripple / ch4`

Interpretation:

- `session-wise zscore` is not a tiny cosmetic change
- it materially changes the detector output
- switching to it means redefining the detector, not just cleaning up scale

### 5. Whether per-session mean/std drift is large on the current root

This was diagnosed with:

- [script_diagnose_e10gb1_session_scale_drift.m](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/tests/test24_detection_alignment_checks/script_diagnose_e10gb1_session_scale_drift.m:1)

Outputs:

- [e10gb1_session_std_by_channel.png](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/tests/test24_detection_alignment_checks/outputs_e10gb1_session_scale_drift/e10gb1_session_std_by_channel.png)
- [e10gb1_session_std_summary.png](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/tests/test24_detection_alignment_checks/outputs_e10gb1_session_scale_drift/e10gb1_session_std_summary.png)
- [trend_summary_table.csv](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/tests/test24_detection_alignment_checks/outputs_e10gb1_session_scale_drift/trend_summary_table.csv)

Main numeric summary on the current root:

- `mean_absmax_overall = 0.047405`
- `mean_range_overall = 0.064349`
- `std_ratio_max_overall = 1.260740`

Interpretation:

- session means are already close to zero
- the more noticeable drift is in session std, not session mean

## Important Caveat About Session Std

Session std is **not** a pure recording-gain metric.

It can reflect:

- baseline/background variation
- slow drift
- artifacts
- real state differences
- event abundance and event amplitude

This matters because a session with more large events can naturally have a
higher overall raw-signal std even if the recording scale did not drift.

Therefore:

- observed session std differences do **not** automatically imply that
  session-wise zscore is justified

## Drift Pattern Interpretation

The current per-session std pattern does **not** look like one clean
global recording drift shared by every channel.

Observed pattern:

- some `hp` channels show fairly strong positive trends with session order
- some `pl` channels are mixed or even trend downward within the `spont`
  block

This suggests:

- there may be some drift-like behavior in part of the recording
- but the pattern is not consistent enough to justify forcing every session
  into its own z-scored coordinate system

## Discussion for Observable Generation

The same zscore question also affects observable generation.

Current reasoning:

- the same electrodes are used
- the full recording spans only a few hours
- for observables, changing each session into a separate normalized space
  risks changing the shared state-space geometry

Therefore:

- **do not** adopt `session-wise zscore` for observable generation
- leave the observable pipeline unchanged for now

In other words:

- `event detection` and `observable generation` should not be forced to use
  the same normalization decision
- at the moment, observable generation should stay as-is

## Current Decision

### Event detection

Use:

- `input_normalization = 'zscore_per_channel'`

Meaning:

- concatenate included sessions first
- then zscore each selected channel once on the full concatenated trace

Do not use:

- `session-wise zscore`

Reason:

- same electrodes, short total recording window
- small session mean drift
- session std differences may partly reflect real activity differences
- session-wise zscore materially changes event results

### Observable generation

Use:

- no change for now

Do not adopt:

- `session-wise zscore`

Reason:

- it would make each session live in a different transformed coordinate
  system and may distort the shared state space

## Current Code State

The event detector was updated so that the default event-detection input
normalization is global per-channel zscore:

- [compute_blp_bandpass_events.m](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/functions/preprocessing/compute_blp_bandpass_events.m:1)

Related scripts explicitly set:

- [script_compute_one_cfg_blp_bandpass_events.m](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/scripts/script_compute_one_cfg_blp_bandpass_events.m:1)
- [script_run_e10gb1_to_event_diversity_windows.m](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/scripts/script_run_e10gb1_to_event_diversity_windows.m:1)

The official current processed event file for `E10gb1` was refreshed with
this metadata:

- `E:\DataPons_processed\e10gb1\event_detection\e10gb1_bandpass_events_3bands.mat`

## Recommended Next Step

Do not spend more time on zscore unless a later modeling result clearly
shows a failure that points back to normalization.

The zscore question is currently settled well enough to proceed with:

- consensus-state definition
- manual ground-truth comparison
- downstream model/observable evaluation
