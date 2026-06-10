# E10gb1 P5 Standardize Condition Comparison Archive

Date: 2026-05-24

Purpose: archive the e10gb1 raw-vs-standardized BLP P4/P5 probe before using
the standardized complex-split branch in P8/P10.

## Scope

Dataset:

- `e10gb1`

Compared BLP conditions:

- raw abs: `abs_projected_vlambda`
- raw complex split: `complex_split_projected_vlambda`
- standardized abs: `abs_projected_vlambda_standardize`
- standardized complex split: `complex_split_projected_vlambda_standardize`

The standardized EDMD outputs were exported from the current standardized BLP
best checkpoints into:

```text
E:\DataPons_processed\derived_autodl_results_standardize\e10gb1\mlp\outputs\
```

The full standardized P5 reduction outputs are under:

```text
E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_reduction\abs_projected_vlambda_standardize\
E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_reduction\complex_split_projected_vlambda_standardize\
```

The standardized P5 peak statistics use the RMS-envelope adaptive activity
policy and are under:

```text
E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_peaks_by_state_rmsenv_adaptive\abs_projected_vlambda_standardize_rmsenv_adaptive\
E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_peaks_by_state_rmsenv_adaptive\complex_split_projected_vlambda_standardize_rmsenv_adaptive\
```

## Scripts And Result Tables

Main scripts added or used:

- `scripts/export_e10gb1_standardized_edmd_outputs.py`
- `scripts/script_run_e10gb1_standardized_p5_ripple_probe.m`
- `scripts/compare_e10gb1_p5_standardize_ripple_selectivity.py`
- `scripts/compare_e10gb1_p5_relaxed_ripple_gamma.py`

Result directories:

```text
results\pipeline5_dimred_component_process_labels_e10gb1_standardize_rmsenv_adaptive_all_components\
results\e10gb1_p5_standardize_ripple_probe_v2\
results\e10gb1_p5_standardize_ripple_probe_v3_allow_gamma\
results\e10gb1_p5_relaxed_ripple_gamma_probe\
results\e10gb1_p5_two_subprocess_candidates\
```

Useful summary files:

```text
results\e10gb1_p5_two_subprocess_candidates\summary.md
results\e10gb1_p5_two_subprocess_candidates\method_k_two_subprocess_candidates.csv
results\e10gb1_p5_two_subprocess_candidates\component_subprocess_flags.csv
```

## Main Findings

Strict ripple-only selectivity became worse after standardization. The earlier
raw abs branch had four strict ripple-selective components, while the
standardized branches did not improve that strict definition.

This strict definition is probably not the right biological target. The current
working target is two subprocesses:

- theta subprocess
- ripple-gamma subprocess

Gamma co-activity with ripple is acceptable. Pure theta activity is the main
thing that should disqualify a ripple-gamma component. Sharp-wave-ripple is not
treated as pure ripple-only, because it can include theta/gamma structure.

Using this relaxed two-subprocess view, standardized complex split is the most
useful e10gb1 branch so far.

## Strong E10gb1 Two-Subprocess Candidates

The criterion here is: within the same P5 method/k setting, there is at least
one theta component and at least one ripple-gamma component with pure theta
inactive.

Standardized complex split satisfies this for many method/k settings:

- SVD: k03-k08
- NMF: k03-k08
- MDS: k04-k08
- UMAP: k05-k08

Compact candidates:

- `complex_split_projected_vlambda_standardize`, SVD k03:
  theta component = 1, ripple-gamma no-pure-theta component = 2
- `complex_split_projected_vlambda_standardize`, NMF k03:
  theta component = 2, ripple-gamma no-pure-theta component = 1

Interpretation: for e10gb1 alone, standardized complex split looks promising
because it can separate a theta-like component from a ripple-gamma-like
component inside the same low-k reduction.

## Current Caveat

The standardized P5 probe initially generated full reduction and peak/selectivity
outputs, but did not generate the P5 raw/dimred thresholded density MAT files.
Those density files are required before P8/P10 can test whether the standardized
complex-split branch also gives stronger or more interpretable BLP-BOLD
cross-correlation.

The next P8/P10 check should therefore first generate:

```text
E:\DataPons_processed\e10gb1\pipeline5_raw_thresholded_density\complex_split_projected_vlambda_standardize_rmsenv_adaptive\
E:\DataPons_processed\e10gb1\pipeline5_dimred_thresholded_density\complex_split_projected_vlambda_standardize_rmsenv_adaptive\
```

Then P8/P10 should use:

```text
blp_density_condition_suffix = 'standardize_rmsenv_adaptive'
density_source_kinds = {'raw_complex_split_density', 'dimred_complex_split_density'}
```

This keeps the standardized probe isolated from the current raw mainline while
making it directly comparable in P8/P10.

## 2026-05-24 P8/P10 Probe Update

The standardized complex-split density branch was generated for P8/P10:

```text
E:\DataPons_processed\e10gb1\pipeline5_raw_thresholded_density\complex_split_projected_vlambda_standardize_rmsenv_adaptive\
E:\DataPons_processed\e10gb1\pipeline5_dimred_thresholded_density\complex_split_projected_vlambda_standardize_rmsenv_adaptive\
```

UMAP k03-k08 initially failed while saving the dimred density MAT files because
the full Windows path length reached the MATLAB/HDF5 `-v7.3` boundary.  The
failed 0-byte files were moved to:

```text
E:\DataPons_processed\_legacy_quarantine\20260524_umap_dimred_density_zero_byte\e10gb1\
```

The density-only script now uses a shorter dimred density filename stem, and
the P8/P10 density resolver ignores 0-byte MAT artifacts.  After rerunning UMAP,
all 24 standardized complex-split dimred density MAT files are present and
non-empty.

P8 was run for current-best e10gb1 BOLD P7 observables using:

```text
xcorr_save_tag = 'xcorr_csplit_standardize_rmsenv_adaptive'
density_source_kinds = {'event_density','raw_complex_split_density','dimred_complex_split_density'}
blp_density_condition_suffix = 'standardize_rmsenv_adaptive'
top_n = 50
```

P10 was run with the same BLP density sources using:

```text
xcorr_save_tag = 'dimred_xcorr_csplit_standardize_rmsenv_adaptive'
feature_names = {'efun_real','deconv_real'}
P9 BOLD method tags = SVD/NMF/MDS/UMAP k05-k08 where available
top_n = 50
```

Probe summary outputs:

```text
results\e10gb1_standardized_csplit_p8_p10_probe\
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8_p10_e10gb1_standardized_csplit_probe\
```

First-pass interpretation:

- P8: dimred standardized complex-split density wins all e10gb1 contexts.
- P10: dimred wins the majority of contexts; raw standardized complex-split
  density also wins a meaningful subset; event density wins only some efun
  contexts.
- The strongest dimred efun hits include both theta-strict and
  ripple-gamma-no-pure-theta components.
- The top deconv-efun dimred hits are dominated by
  ripple-gamma-no-pure-theta components in both P8 and P10.

Raw BLP efun density distribution check:

Important correction on 2026-05-25:

- The first raw-timescale plot used `mode_metadata.timescale_sec`, which older
  P5 density MAT files defined from the bilinear continuous eigenvalue
  transform. That is not the preferred slow/fast-mode definition for the
  P8/P10 raw efun comparison.
- The corrected v3 plot uses `timescale_sec_preferred`, computed from the
  discrete Koopman eigenvalue as `-dt/log(abs(lambda))`; frequency uses
  `abs(angle(lambda))/(2*pi*dt)`.
- Ignore the old figure
  `raw_blp_efun_timescale_ecdf_by_topN.png` for interpretation. Use
  `raw_blp_efun_timescale_ecdf_by_topN_v3_discrete_log.png` and
  `raw_blp_efun_distribution_summary_v3_discrete_log.csv`.

- The raw BLP efun density hits split very strongly by BOLD feature family.
- For top10/top20 raw hits, BOLD `efun` coupling is concentrated in low-index,
  slow raw BLP efuns, mainly raw efun indices 13/14 with median timescale about
  70 s under the corrected discrete-log decay definition.
- For top10/top20 raw hits, BOLD `deconv_efun` coupling is concentrated in
  high-index, very fast raw BLP efuns, around raw efun indices 255-273 with
  sub-millisecond-to-millisecond-scale discrete-log decay estimates.
- At top50, BOLD `efun` becomes more mixed, especially in P10, but
  `deconv_efun` remains strongly biased toward high-index fast raw BLP efuns.

This supports using standardized complex split as a serious candidate branch,
but it is still e10gb1-only and therefore not a cross-session conclusion.
