# Plotting Trace Clipping Removal Archive

Date: 2026-05-28

## Reason

K13m17 top-window inspection suggested that many theta-event raw traces looked
saturated.  The first QC pass did not show large-scale hard raw-data clipping,
but the plotting code was still visually clipping traces.  That made it hard to
tell whether a signal was truly saturated or only flattened by the display
layer.

Therefore the project plotting policy is now:

```text
Raw / observable / component trace plots should not hard-clip amplitude by
default.  Large excursions should remain visible, and the y-axis should expand
to show them.
```

## Scope

This is a visualization/provenance change.  It does not change P1-P10 numeric
outputs, event detection, eigenfunction reduction, density construction, peak
statistics, or cross-correlation values.

It changes how browsing/QC figures display time-series traces.

## What Changed

### BLP top-window and segment plots

Changed files:

```text
functions/plottings/prepare_blp_plot_data.m
functions/plottings/build_blp_plot_window_cache.m
functions/plottings/plot_blp_segment_with_spectrogram.m
functions/plottings/plot_blp_segment_with_spectrogram_and_koopman.m
```

Behavior:

- `trace_clip` now defaults to `Inf`.
- Old `cfg.plot.trace_clip` values are kept for backward-compatible structs but
  ignored by the shared BLP plot-data path.
- The old display operation
  `max(min(trace, trace_clip), -trace_clip)` was removed.
- Raw trace-panel y-limits now expand from the actual un-clipped displayed
  trace values.

### BOLD pre-ResKoopNet QC traces

Changed file:

```text
functions/plottings/plot_bold_pre_reskoopnet_qc.m
```

Behavior:

- Representative BOLD observable traces are still robust-zscored for layout.
- They are no longer clipped to `[-4, 4]`.

### BOLD observable/eigenfunction dimred summary

Changed file:

```text
functions/postprocessing/plot_bold_observable_efun_dimred_summary.m
```

Behavior:

- The helper no longer clips z-scored observable traces to `[-5, 5]`.

### P2/P5 top30 diagnostic script

Changed file:

```text
scripts/plot_p2_p5_top30_windows.py
```

Behavior:

- P5 component-activity traces are still robust-scaled for comparability within
  a window.
- They are no longer clipped to `[0, 1]`.
- Per-panel y-limits now expand to show the un-clipped trace.

### Script entry defaults

The following plotting entry scripts no longer set `cfg.plot.trace_clip = 4`;
they set `Inf` instead:

```text
scripts/script_plot_blp_segment_with_events.m
scripts/script_plot_blp_segment_with_spectrogram.m
scripts/script_plot_blp_segment_with_spectrogram_and_koopman.m
scripts/script_plot_top_consensus_state_diversity_windows.m
scripts/script_plot_top_consensus_event_diversity_windows.m
scripts/script_plot_top_consensus_event_diversity_windows_e10gb1.m
scripts/script_rerun_e10gh1_cfg_to_observables_and_events.m
```

## What Did Not Change

The following are not display trace clipping and were intentionally left alone:

- Heatmap color limits (`clim`, `caxis`) for correlations, spectrograms,
  activation maps, state-space labels, or density maps.
- Adaptive envelope window bounds such as
  `np.clip(win_sec, ENVELOPE_MIN_SEC, ENVELOPE_MAX_SEC)`.
- Boundary clipping of event windows to a valid time range.
- Colormap/RGB bounds.
- Algorithmic options whose names contain `clip`, such as `clip_zero`.

## Regenerated K13m17 Figures

After removing trace clipping, K13m17 top-window figures were regenerated.

P2 official top-window figures:

```text
E:\DataPons_processed\k13m17\pipeline2_figures_consensus_state_top_window_plots\
```

P2 lightweight event-raster top30 diagnostics:

```text
E:\DataPons_processed\k13m17\pipeline2_figures_consensus_state_top30_event_raster\
```

P5 component-activity top30 diagnostics:

```text
E:\DataPons_processed\k13m17\pipeline5_component_activity_top30_windows\complex_split_projected_vlambda_standardize\
```

## Interpretation Rule

Figures generated before this change may visually flatten large events even
when the raw data are not hard-saturated.  For raw trace shape, top-window
inspection, and K13m17 event-quality review, prefer figures generated on or
after 2026-05-28 with the no-trace-clipping policy.

P11 or any future figure audit should treat pre-2026-05-28 raw/trace browsing
figures as potentially stale for visual-amplitude interpretation.  They can
still be useful for event timing, labels, and provenance if the underlying
numeric source is current.

## Verification

Checks performed:

```text
python -m py_compile scripts/plot_p2_p5_top30_windows.py scripts/qc_blp_raw_saturation_by_p2_events.py
MATLAB checkcode on changed plotting functions
K13m17 P2/P5 top-window PNG nonblank checks
```

Remaining code references to `clip` after this change are non-display uses as
listed above.
