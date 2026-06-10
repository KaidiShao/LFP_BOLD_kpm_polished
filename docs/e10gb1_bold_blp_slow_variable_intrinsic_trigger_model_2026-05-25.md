# E10gb1 BOLD-BLP Slow Variable vs Intrinsic Trigger Model

Date: 2026-05-25

Status: working interpretation from the E10gb1 condition-comparison probe.  This
is not yet a cross-session conclusion.

## Core Question

The current scientific question is whether BLP Koopman eigenfunction-derived
densities can be interpreted as subprocesses that couple to different aspects
of BOLD dynamics.

The working distinction is:

- `BOLD efun x BLP density`: coupling to slow BOLD latent-state variables.
- `BOLD deconv efun x BLP density`: coupling to residual/intrinsic BOLD
  perturbations.

The important new interpretation is that these two BOLD-side targets may map
onto two different BLP subprocesses:

- slow theta/state-like variable
- fast ripple-gamma/intrinsic-trigger-like activity

This model became much clearer after plotting actual top-xcorr time series,
not just rank tables.

## Four BLP Conditions Compared

The four E10gb1 BLP conditions are:

```text
abs_projected_vlambda
complex_split_projected_vlambda
abs_projected_vlambda_standardize
complex_split_projected_vlambda_standardize
```

Conceptually this is:

```text
BLP observable type: abs vs complex split
BLP normalization: raw vs standardized
```

The standardized branches were exported from standardized BLP EDMD outputs and
then run through P5.  The currently most important condition is:

```text
complex_split_projected_vlambda_standardize
```

For the activity-density and selectivity analyses, the current mainline
activity definition is RMS-envelope with adaptive eigenvalue-derived window:

```text
*_rmsenv_adaptive
*_standardize_rmsenv_adaptive
```

## Data Scope

Dataset:

```text
e10gb1
```

What is relatively complete:

- P5 selectivity comparison across the four BLP conditions.
- P8/P10 standardized csplit probe.
- csplit nonstandard vs csplit standardized raw-density xcorr comparison.
- raw BLP efun density top-xcorr time-series comparison.
- dimred BLP efun density top-xcorr time-series comparison.

What is not complete yet:

- Equal-depth P8/P10 topN and time-series comparison for all four BLP
  conditions.
- Cross-session validation across all datasets.
- Final P11 one-click report integration for this exact model.

## Experiment 1: P5 Selectivity Across Four Conditions

Purpose:

Determine whether P5 dimred BLP eigenfunction components can separate a
theta-like subprocess from a ripple-gamma-like subprocess.

Main outputs:

```text
results\e10gb1_p5_two_subprocess_candidates\summary.md
results\e10gb1_p5_two_subprocess_candidates\method_k_two_subprocess_candidates.csv
results\e10gb1_p5_two_subprocess_candidates\component_subprocess_flags.csv
```

Main scripts:

```text
scripts\export_e10gb1_standardized_edmd_outputs.py
scripts\script_run_e10gb1_standardized_p5_ripple_probe.m
scripts\compare_e10gb1_p5_standardize_ripple_selectivity.py
scripts\compare_e10gb1_p5_relaxed_ripple_gamma.py
```

Key design change:

Strict ripple-only selectivity was too restrictive.  The current biological
target is two subprocesses:

- theta subprocess
- ripple-gamma subprocess

Gamma co-activity is allowed for the ripple-gamma subprocess.  Sharp-wave-ripple
is not treated as pure ripple-only because it can contain theta/gamma structure.
Pure theta activity is the main exclusion criterion for a ripple-gamma component.

Main result:

`complex_split_projected_vlambda_standardize` is the strongest E10gb1 condition
so far for separating theta and ripple-gamma subprocesses.

Examples:

- standardized csplit, SVD k03:
  - theta component = 1
  - ripple-gamma no-pure-theta component = 2
- standardized csplit, NMF k03:
  - theta component = 2
  - ripple-gamma no-pure-theta component = 1

Broader E10gb1 standardized csplit candidates:

- SVD: k03-k08
- NMF: k03-k08
- MDS: k04-k08
- UMAP: k05-k08

Interpretation:

Standardized complex split appears to make the fast/ripple-gamma component more
visible while preserving a theta-like component in the same low-k reduction.

## Experiment 2: Raw BLP Efun Density Xcorr Distributions

Purpose:

Compare how raw BLP efun density couples to BOLD efun versus BOLD deconv efun,
especially under nonstandard vs standardized csplit.

Main output tables:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_xcorr_distribution_std_vs_nonstandard_long.csv
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_xcorr_distribution_std_vs_nonstandard_summary.csv
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_xcorr_topN_compact_interpretation.csv
```

Main figures:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_ecdf_std_vs_nonstandard_top50.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_topN_median_std_vs_nonstandard.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_topN_delta_std_minus_nonstandard.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_threshold_fraction_by_topN.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_xcorr_exact_rank_profile_std_vs_nonstandard.png
```

Important correction:

Raw efun timescale interpretation now uses the preferred discrete eigenvalue
decay timescale:

```text
tau = -dt / log(abs(lambda_discrete))
```

Earlier bilinear continuous estimates are legacy/fallback provenance only.

Main result:

- `P8 efun`: standardized csplit is weaker than nonstandard at top ranks.
- `P8 deconv_efun`: standardized and nonstandard are similar in strength.
- `P10 efun`: standardized csplit becomes stronger.
- `P10 deconv_efun`: standardized csplit is weaker in raw-density xcorr.

Interpretation:

The xcorr magnitude alone is not enough.  The target class matters:

- efun-side hits tend to be slow/state-like.
- deconv-side hits tend to be fast/activity-like.

## Experiment 3: Raw Top5 Xcorr Time-Series Inspection

Purpose:

Check what the strongest raw BLP efun density xcorr hits actually look like in
time.  This was the key step that made the interpretation clear.

Main script:

```text
scripts\plot_e10gb1_raw_csplit_top5_xcorr_time_series.m
```

Output table:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\raw_csplit_top5_xcorr_timeseries_rows.csv
```

Output directory:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\raw_csplit_top5_xcorr_timeseries\
```

Figures:

```text
raw_csplit_top5_xcorr_timeseries__nonstandard__P8__efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P8__efun.png
raw_csplit_top5_xcorr_timeseries__nonstandard__P8__deconv_efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P8__deconv_efun.png
raw_csplit_top5_xcorr_timeseries__nonstandard__P10__efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P10__efun.png
raw_csplit_top5_xcorr_timeseries__nonstandard__P10__deconv_efun.png
raw_csplit_top5_xcorr_timeseries__standardized__P10__deconv_efun.png
```

Plot convention:

- blue = raw BLP efun density
- red = BOLD-side target shifted to the peak lag
- both are z-scored
- each figure shows top5 xcorr rows for one condition x pipeline x feature
  family

Key visual result:

`BOLD efun x raw BLP density` often captures slow, session-scale structure.
This looks like a slow state variable or session-wise drift.

`BOLD deconv efun x raw BLP density` captures faster activity and burst-like
structure.  This looks more like residual/intrinsic perturbation.

Interpretation:

This supports separating the two BOLD-side targets:

- BOLD efun: slow state variable
- BOLD deconv efun: intrinsic trigger / residual perturbation

## Experiment 4: Dimred BLP Efun Density Labels

Purpose:

Check whether the raw efun slow/fast distinction is preserved after BLP
dimension reduction.

Main output tables:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\top_dimred_rows.csv
results\e10gb1_standardized_csplit_p8_p10_probe\best_dimred_method_k_rows.csv
results\e10gb1_standardized_csplit_p8_p10_probe\all_top_rows.csv
```

Important existing figure:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\p8_p10_dimred_two_subprocess_label_counts_by_family_top50.png
```

Key label result:

Top dimred hits show a strong feature-family separation:

- `efun` top hits are mixed but include strong `theta_strict` components.
- `deconv_efun` top hits are dominated by
  `ripple_gamma_no_pure_theta`.

This is important because it means the two-subprocess structure is not just a
raw-mode artifact.  It is visible in dimred BLP density as well.

## Experiment 5: Dimred Top5 Xcorr Time-Series Inspection

Purpose:

Directly inspect whether dimred BLP efun density top hits visually match the
slow-variable vs intrinsic-trigger model.

Main script:

```text
scripts\plot_e10gb1_dimred_csplit_top5_xcorr_time_series.m
```

Output table:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\dimred_csplit_top5_xcorr_timeseries_rows.csv
```

Output directory:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\dimred_csplit_top5_xcorr_timeseries\
```

Figures:

```text
dimred_csplit_top5_xcorr_timeseries__P8__efun.png
dimred_csplit_top5_xcorr_timeseries__P8__deconv_efun.png
dimred_csplit_top5_xcorr_timeseries__P10__efun.png
dimred_csplit_top5_xcorr_timeseries__P10__deconv_efun.png
```

Top5 observations:

- P8 efun top5:
  - all are `mds k08 component 3`
  - label = `theta_strict`
  - visually slow/state-like
- P10 efun top5:
  - first three are `mds k08 component 3`
  - label = `theta_strict`
  - remaining two are `nmf k03 component 1`
  - label = `ripple_gamma_no_pure_theta`
  - interpretation: efun side is theta/state-dominant but can mix in
    ripple-gamma components
- P8 deconv_efun top5:
  - dominated by `umap k03/k04 component 2` and `mds k03 component 3`
  - label = `ripple_gamma_no_pure_theta`
- P10 deconv_efun top5:
  - dominated by `umap k03/k04 component 2`
  - label = `ripple_gamma_no_pure_theta`

Interpretation:

Dimred BLP efun density supports the same distinction as raw BLP efun density:

```text
BOLD efun        -> slow/theta/state-like subprocess
BOLD deconv efun -> fast/ripple-gamma/intrinsic-trigger subprocess
```

The deconv side is cleaner than the efun side.

## Working Scientific Model

The current model is:

```text
theta-like slow BLP subprocess
    -> tracks or modulates BOLD latent state variable
    -> visible through BOLD efun x BLP density

ripple-gamma-like fast BLP subprocess
    -> aligns with residual BOLD perturbation / intrinsic trigger
    -> visible through BOLD deconv efun x BLP density
```

In words:

The BOLD efun target appears to capture a slow state variable.  The BLP
subprocess aligned with it is slow and often theta/state-like.

The BOLD deconv efun target is closer to a residual after accounting for slow
BOLD dynamics.  The BLP subprocess aligned with it is fast and
ripple-gamma-like, which makes it a plausible intrinsic trigger or perturbation
signal.

This is the cleanest interpretation so far.

## Why This Is Better Than A Single Xcorr Ranking

Combining efun and deconv efun into one ranking mixes two different mechanisms:

- slow latent-state coupling
- fast residual-trigger coupling

Therefore P8/P10 should not be summarized only as "which density has the
largest xcorr".  The BOLD-side target must be part of the interpretation.

Recommended first-level split:

```text
efun        -> slow variable / theta-state analysis
deconv_efun -> intrinsic trigger / ripple-gamma analysis
```

Only after this split should one compare:

- BLP observable type: abs vs complex split
- BLP normalization: raw vs standardized
- dimred method: SVD/NMF/MDS/UMAP
- component count: k03-k08
- BOLD observable type
- raw vs dimred BLP density

## Current Caveats

1. This is still mainly an E10gb1 result.

2. The deepest P8/P10 inspection currently focuses on csplit nonstandard versus
   csplit standardized.  The P5 selectivity comparison covers all four BLP
   conditions more completely than the P8/P10 visual inspection.

3. Some top5 rows are near-duplicates across BOLD modes or method/k settings.
   This is useful for seeing stability, but future summaries should also include
   de-duplicated views by density component.

4. High xcorr does not imply causality.  The current interpretation is
   mechanistic/phenomenological: slow variable versus intrinsic trigger.

5. The standardized csplit branch looks promising, but it still needs
   cross-session validation before it becomes the main analysis choice.

## Recommended Next Checks

1. Repeat the same dimred top5 time-series and label analysis across all
   available datasets.

2. For each dataset, explicitly test whether:

```text
efun top hits are enriched for theta/state labels
deconv_efun top hits are enriched for ripple-gamma labels
```

3. Add a P11 section that reports:

- label composition by BOLD-side target
- topN xcorr strength by BOLD-side target
- topN raw efun timescale distribution by BOLD-side target
- top5 or de-duplicated top component time-series snapshots

4. Extend the same analysis from standardized csplit to all four BLP
   conditions, but keep the output grouped by BOLD-side target.

5. Use this model to guide parameter choice:

- choose BLP observable/normalization that separates theta-state and
  ripple-gamma-trigger components;
- choose dimred method/k that gives stable theta and ripple-gamma labels;
- choose BOLD observable type separately for slow-variable and intrinsic-trigger
  interpretations.

## Experiment 6: BOLD Observable-Separated Check

### Motivation

The previous slow-variable versus intrinsic-trigger model is useful, but it
cannot be interpreted without separating BOLD observable type.  The `pv` prefix
in P8/P10 run tags is shorthand for `projected_vlambda`; it is not itself a new
biological observable.  The current observable modes are:

```text
global_svd100
gsvd100_ds
HP_svd100
roi_mean
```

`HP_svd100` should be treated separately from the more global/state-like
observables because it is a local/high-pass BOLD representation.

### Script And Outputs

Main script:

```text
scripts\plot_e10gb1_observable_separated_subprocess_probe.m
```

Input table:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\all_top_rows.csv
```

Output tables:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\observable_separated_dimred_label_summary_top50.csv
results\e10gb1_standardized_csplit_p8_p10_probe\observable_separated_raw_timescale_summary_top50.csv
results\e10gb1_standardized_csplit_p8_p10_probe\observable_separated_dimred_top5_rows.csv
```

Output figures:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\observable_separated_subprocess_probe\observable_separated_dimred_label_fraction_top50.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\observable_separated_subprocess_probe\observable_separated_dimred_xcorr_strength_top50.png
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\observable_separated_subprocess_probe\observable_separated_raw_timescale_top50.png
```

### Main Findings

`BOLD deconv_efun -> ripple-gamma / intrinsic trigger` mostly still holds,
especially for:

```text
P8 deconv_efun:
    gsvd100_ds
    HP_svd100
    roi_mean

P10 deconv_efun:
    HP_svd100
```

These conditions are dominated by `ripple_gamma_no_pure_theta` dimred BLP
density labels among the high-xcorr rows.  This supports the interpretation
that deconv BOLD eigenfunctions are closer to residual/intrinsic perturbation
targets.

However, `P8 global_svd100 deconv_efun` is an exception in the current E10gb1
probe; it has more theta/other structure.  Therefore deconv-side results still
need observable-specific checks.

The simpler statement `BOLD efun -> slow theta/state variable` only holds for
some observable modes.  It is clearest for:

```text
global_svd100 efun
```

`HP_svd100 efun` does not look like a slow state variable in the same sense.
It is more local/fast and ripple-gamma-like.  `roi_mean` and `gsvd100_ds` are
mixed: they can show theta/state-like top hits, but top50 and P10 summaries can
shift toward ripple-gamma.

### Revised Model

The more accurate model is now:

```text
global/state-like BOLD observable efun
    -> slow theta/state variable

local/HP BOLD observable efun
    -> not necessarily slow-state;
       can reflect local fast/ripple-gamma-like activity

BOLD deconv_efun
    -> often residual / intrinsic perturbation / trigger,
       but still requires observable-specific validation
```

Therefore P8/P10 summaries should be split at least by:

```text
pipeline: P8 vs P10
BOLD feature family: efun vs deconv_efun
BOLD observable: global_svd100 / gsvd100_ds / HP_svd100 / roi_mean
```

Only after this split should BLP observable type, normalization, dimred method,
and component count be compared.

### Current Interpretation

The best current E10gb1 interpretation is not simply "efun is slow and deconv is
fast."  It is:

```text
Different BOLD observables expose different levels of BOLD dynamics.

Global/state-like BOLD efun is the best target for slow theta/state subprocesses.

HP/local BOLD efun and many deconv_efun targets are better matched to
ripple-gamma / intrinsic-trigger subprocesses.

complex_split_projected_vlambda_standardize remains the most promising E10gb1
BLP condition because it can expose both theta and ripple-gamma components in
P5 and then separate the corresponding BOLD-side targets in P8/P10.
```
