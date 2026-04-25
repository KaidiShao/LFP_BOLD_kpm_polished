# Test24 Detection Alignment Checks

This test branch isolates why `event_detection` results may differ between:

- the older `test23_*bandpass_threshold_pipeline*.m` scripts
- the current `compute_blp_bandpass_events.m` pipeline

Current scope:

- `e10gb1` only
- compare legacy prepared MAT input vs current `load_blp_dataset(cfg_E10gb1())`
- compare fixed `Fs = 667` vs `Fs = 1 / median(session_dx)`
- compare one-shot concatenated detection vs session-split detection

Canonical entry point:

- `script_compare_e10gb1_detection_variants.m`
- `script_compare_e10gb1_zscore_scopes.m`

Outputs are written into:

- `tests/test24_detection_alignment_checks/outputs_e10gb1`

Current finding for `e10gb1` on 2026-04-16:

- The old `test23`-style prepared MAT does have an exact upstream source:
  loading from
  `D:\Onedrive\experimental_event_data\ripple_monkeys\E10.gb1\blp\`
  with `sr+pl`, session ids
  `[1:5, 6:7, 9:13, 20:25]`, and `blp.dat(:,:,1)` reproduces
  `Nikos_ripples_e10gb1_sr_pl.mat` exactly (`maxabs = 0`).
- Therefore the earlier mismatch is not because the comparison code was
  broken, and not because `test23` applied an extra normalization step in
  the detector script. The mismatch comes from using a different raw-data
  root.
- The current repo default root
  `E:\DataPons\E10.gb1\blp\`
  is a different file version from the older copies in
  `D:\Onedrive\experimental_event_data\...`
  and `E:\DataPons\E10.gb1.old\...`.
  Even the same session file has a different size on disk.
- `loader_data_segmented_dxfs_current` reproduces the saved current result
  exactly, so this test folder now provides a faithful readout of the
  current detector behavior.
- The largest mismatch versus the legacy `test23` reference comes from the
  input source / channel layout difference:
  legacy prepared MAT uses the older normalized-scale `sr + pl` blp copy,
  while the current cfg/loader uses the newer `E:\DataPons` `hp + pl`
  copy.
- Session splitting by itself is a small effect:
  `legacy_mat_segmented_fixed667` still matches `legacy_mat_flat_fixed667`
  with recall `0.9893`.
- Changing from fixed `Fs = 667` to `Fs = 1 / median(session_dx)` for the
  current loader data adds a noticeable but smaller shift on top of the
  data-source mismatch.

Useful outputs:

- `aggregate_match_table.csv`
- `data_alignment_table.csv`
- `saved_result_alignment_table.csv`
- `variant_summary_table.csv`
- `outputs_e10gb1_zscore_scope/aggregate_table.csv`
- `outputs_e10gb1_zscore_scope/threshold_table.csv`

Additional finding on 2026-04-17:

- On the current `cfg_E10gb1()` raw-data root, session means are already
  close to zero, but per-session standard deviations drift by up to about
  `1.26x` across sessions.
- Comparing detector input normalization scopes:
  `global_zscore` versus `sessionwise_zscore`
  changes event detection materially even when all other detector settings
  are held fixed.
- Aggregate peak overlap on the current root:
  global total = `26742`, sessionwise total = `26465`,
  shared = `24776`, mean per-channel/band Jaccard = `0.86456`.
- The largest differences concentrate in a few band/channel pairs, such as
  `gamma/ch4` and `ripple/ch4`.
