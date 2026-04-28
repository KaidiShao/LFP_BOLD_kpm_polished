# BLP Event Pipeline: Storage And Flow

This note summarizes the current event-detection pipeline organization for processed BLP datasets and records the intended storage layout.

## Scope

This document covers the event-analysis line:

1. raw BLP loading
2. bandpass event detection
3. event density
4. consensus states
5. consensus-state summary
6. diversity / variability window analysis
7. top-window plotting

This is separate from the ResKoopNet / spectrogram / Koopman postprocessing line, except when later visualization compares them.

Current interpretation of the two downstream window branches:

- `event_diversity_windows` is an early diagnostic / quick-check branch
- `consensus_state_diversity_windows` is the primary long-term downstream branch

For the primary consensus-state branch, spectrogram generation should be treated
as plotting support for top-window figures, not as part of the core
event/state logic. In practice, the canonical plotting path should reuse saved
spectrogram files when they already exist and should fail immediately if the
required plotting-support spectrogram is missing.

## Related Notes

- E10gb1 zscore discussion and current decision:
  [e10gb1_zscore_discussion_and_decision_2026-04-17.md](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/docs/e10gb1_zscore_discussion_and_decision_2026-04-17.md:1)

## Canonical Storage Root

Processed event-line outputs should live under:

- `E:\DataPons_processed\<dataset>\...`

Current datasets:

- `E:\DataPons_processed\f12m01`
- `E:\DataPons_processed\e10gb1`

## Canonical Directory Layout

For a dataset such as `f12m01` or `e10gb1`, the event-analysis outputs should use:

- `event_detection`
- `event_density`
- `consensus_states`
- `consensus_state_summary`
- `event_diversity_windows`
- `consensus_state_diversity_windows` (optional; currently used for `e10gb1`)

Examples:

- `E:\DataPons_processed\f12m01\event_detection`
- `E:\DataPons_processed\f12m01\event_diversity_windows`
- `E:\DataPons_processed\e10gb1\consensus_state_diversity_windows`

## Current Dataset State

### `f12m01`

Currently present:

- `event_detection`
- `event_density`
- `consensus_states`
- `consensus_state_summary`
- `event_diversity_windows`
- `koopman_postprocessing`
- `reskoopnet_dictionary`
- `spectrograms`

Event-line files currently present:

- `event_detection\f12m01_bandpass_events_3bands.mat`
- `event_density\f12m01_event_density_2s.mat`
- `consensus_states\f12m01_consensus_states_min4ch.mat`
- `consensus_state_summary\f12m01_consensus_state_type_summary.mat`
- `consensus_state_summary\f12m01_consensus_state_type_summary.csv`
- `event_diversity_windows\f12m01_event_diversity_windows_5000samp.mat`
- `event_diversity_windows\f12m01_event_diversity_windows_5000samp.csv`
- `event_diversity_windows\f12m01_event_diversity_windows_5000samp_top.csv`
- `event_diversity_windows\top_window_plots\...`

### `e10gb1`

Currently present:

- `event_detection`
- `event_density`
- `consensus_states`
- `consensus_state_summary`
- `event_diversity_windows`
- `consensus_state_diversity_windows`
- `reskoopnet_dictionary`
- `spectrograms`

Event-line files currently present:

- `event_detection\e10gb1_bandpass_events_3bands.mat`
- `event_density\e10gb1_event_density_2s.mat`
- `consensus_states\e10gb1_consensus_states_min5ch.mat`
- `consensus_state_summary\e10gb1_consensus_state_type_summary.mat`
- `consensus_state_summary\e10gb1_consensus_state_type_summary.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_5000samp.mat`
- `event_diversity_windows\e10gb1_event_diversity_windows_5000samp.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_5000samp_top.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp.mat`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_top.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_globalwin.mat`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_globalwin.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_globalwin_top.csv`
- `event_diversity_windows\top_window_plots\...`
- `consensus_state_diversity_windows\e10gb1_consensus_state_diversity_windows_6000samp_globalwin.mat`
- `consensus_state_diversity_windows\e10gb1_consensus_state_diversity_windows_6000samp_globalwin.csv`
- `consensus_state_diversity_windows\e10gb1_consensus_state_diversity_windows_6000samp_globalwin_top.csv`
- `consensus_state_diversity_windows\top_window_plots\...`

## Pipeline Flow

### Step 1. Raw BLP

Load raw concatenated BLP:

- `io_raw.load_blp_dataset(cfg)`

This is the only required raw-data input stage for the event line.

### Step 2. Event Detection

Run:

- `compute_blp_bandpass_events(D, cfg, output_root, params)`

Output:

- `event_detection\<dataset>_bandpass_events_3bands.mat`

This file is the main saved event-detection result and is later reloaded as `R`.

For the current E10gb1 normalization discussion and the decision to keep
global per-channel zscore for event detection, see:

- [e10gb1_zscore_discussion_and_decision_2026-04-17.md](/D:/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/docs/e10gb1_zscore_discussion_and_decision_2026-04-17.md:1)

### Step 3. Event Density

Load saved event result:

- `io_results.load_event_results(cfg, output_root, event_input)`

Then run:

- `compute_blp_event_density(cfg, output_root, R, params, source_event_file)`

Output:

- `event_density\<dataset>_event_density_2s.mat`

### Step 4. Consensus States

Load saved event result:

- `io_results.load_event_results(cfg, output_root, event_input)`

Then run:

- `compute_blp_consensus_states(cfg, output_root, R, params, source_event_file)`

Output:

- `consensus_states\<dataset>_consensus_states_<threshold>.mat`

Examples:

- `f12m01_consensus_states_min4ch.mat`
- `e10gb1_consensus_states_min5ch.mat`

### Step 5. Consensus-State Summary

Load saved consensus-state result:

- `io_results.load_consensus_state_results(cfg, output_root, consensus_input)`

Then run:

- `summarize_blp_consensus_state_types(cfg, output_root, C, params, source_consensus_file)`

Outputs:

- `consensus_state_summary\<dataset>_consensus_state_type_summary.mat`
- `consensus_state_summary\<dataset>_consensus_state_type_summary.csv`

### Step 6. Diversity Windows

Load saved consensus-state result:

- `io_results.load_consensus_state_results(cfg, output_root, consensus_input)`

Then run:

- `analyze_blp_consensus_event_diversity_windows(cfg, output_root, C, params, source_consensus_file)`

Outputs:

- `event_diversity_windows\<dataset>_event_diversity_windows_5000samp.mat`
- `event_diversity_windows\<dataset>_event_diversity_windows_5000samp.csv`
- `event_diversity_windows\<dataset>_event_diversity_windows_5000samp_top.csv`

For `e10gb1`, the canonical diversity output is now:

- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_globalwin.mat`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_globalwin.csv`
- `event_diversity_windows\e10gb1_event_diversity_windows_6000samp_globalwin_top.csv`

Top-window plots should live under:

- `event_diversity_windows\top_window_plots`

For `e10gb1`, these top-window plots and matching reference figures should be interpreted relative to the canonical `6000samp_globalwin` ranking, not the older session-local rankings.

### Step 7. Consensus-State Diversity Windows

This is currently implemented for `e10gb1`:

- `analyze_blp_consensus_state_diversity_windows(...)`

Outputs:

- `consensus_state_diversity_windows\<dataset>_consensus_state_diversity_windows_6000samp_globalwin.mat`
- `consensus_state_diversity_windows\<dataset>_consensus_state_diversity_windows_6000samp_globalwin.csv`
- `consensus_state_diversity_windows\<dataset>_consensus_state_diversity_windows_6000samp_globalwin_top.csv`
- `consensus_state_diversity_windows\top_window_plots\...`

This should be treated as the preferred consensus-state window-selection line when the goal is to find windows containing the largest variety of derived state labels.

## Main Script Entry Points

Single-stage scripts:

- `scripts/script_compute_one_cfg_blp_bandpass_events.m`
- `scripts/script_compute_one_cfg_blp_event_density.m`
- `scripts/script_compute_one_cfg_blp_consensus_states.m`
- `scripts/script_summarize_one_cfg_blp_consensus_state_types.m`

Canonical quick-check branch entries:

- `scripts/script_run_one_cfg_to_event_diversity_windows.m`
- `scripts/script_run_cfgs_to_event_diversity_windows.m`
- `scripts/script_run_e10gb1_to_event_diversity_windows.m`

Fixed-dataset examples:

- `scripts/script_run_e10gb1_to_consensus_state_top_windows.m`
- `scripts/script_run_e10gb1_to_event_diversity_windows.m`

## Current Consistency Check

The on-disk dataset layouts are already mostly consistent.

Consistent items:

- both datasets use the same core event-line directory names
- both datasets keep diversity outputs inside `event_diversity_windows`
- both datasets keep top-window diversity plots inside `event_diversity_windows\top_window_plots`

Known differences that are currently expected:

- `f12m01` uses `min4ch` in the saved consensus filename
- `e10gb1` uses `min5ch` in the saved consensus filename
- `consensus_state_diversity_windows` currently exists only for `e10gb1`

Known legacy / non-canonical diversity outputs that may still remain on disk:

- older `5000samp` diversity files
- `6000samp` files generated in session-local mode

For `e10gb1`, prefer the `6000samp_globalwin` diversity outputs when comparing top windows or matching reference figures.

## Current Naming / Path Convention

The codebase now resolves the processed-data root through:

- `io_project.get_project_processed_root()`

Current default behavior:

- Windows default: `E:\DataPons_processed\`
- environment override: `LFP_BOLD_KPM_PROCESSED_ROOT`

This means single-stage scripts and helper functions should now follow the same processed-data root by default.

## Recommended Next Step

For future event-line work, prefer this interpretation order:

1. event detection
2. consensus states
3. diversity windows
4. variability windows
5. top-window plotting

The likely next improvement is to promote variability-based window selection into the main event-analysis flow, rather than treating it as a side branch.
