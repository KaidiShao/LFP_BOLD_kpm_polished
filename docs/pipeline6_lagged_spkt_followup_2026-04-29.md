# Pipeline 6 Lagged SPKT Follow-Up

Date: 2026-04-29

## Decision

For the delivery target of 2026-04-30 16:00, do not treat lagged
SPKT-residual cross-correlation as required for Pipeline 6 completeness.

## Why It Was Deferred

- The canonical SPKT stage computes same-bin correlation after alignment; it
  does not sweep over positive/negative lags.
- The lagged analysis is a separate optional branch implemented by
  `functions/postprocessing/compute_spkt_residual_lagged_cross_correlation.m`.
- The current lagged entry is still dataset-specific:
  `scripts/script_plot_e10gb1_spkt_residual_lagged_cross_correlation.m`.
- Shipping two datasets with stable top-window, timescale, SPKT zero-lag, and
  MUA outputs is higher priority than folding lagged SPKT into the canonical
  runner right now.

## Current Pipeline 6 Completion Standard

For the near-term run target, treat a dataset as complete when it has:

1. top-window postprocessing outputs
2. timescale diagnostics
3. SPKT zero-lag residual cross-correlation outputs
4. MUA residual cross-correlation outputs

Lagged SPKT remains optional until the canonical Pipeline 6 runner is updated.

## Follow-Up Task

When Pipeline 6 cleanup resumes, add lagged SPKT as an optional canonical
stage:

1. add a `do_spkt_lagged_cross_correlation` flag to the canonical runner
2. run lagged analysis from the current run's SPKT result, not from a
   "latest MAT" directory guess
3. write lagged MAT/CSV/overview paths into the Pipeline 6 manifest
4. decide whether lagged summary figures should also sync into
   `summary_figures`

## Revisit Trigger

Revisit this after at least two datasets have been regenerated successfully
through the current canonical Pipeline 6 path.
