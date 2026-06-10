# RMS Envelope Adaptive P5/P8/P10 Run Log - 2026-05-23

This document records what was implemented and run in the recent
`rmsenv_adaptive` pass.  It is a provenance log, not a scientific conclusion.
The strongest top-xcorr rows are only one small diagnostic view; they are not
the final analysis design.

## Motivation

The immediate problem was eigenfunction polarity/sign ambiguity.  If a raw or
dimred eigenfunction is thresholded directly, a sign flip can reverse what is
called active.  The temporary working rule for this pass was:

- Use magnitude-based activity for all LFP/BLP eigenfunction density and peak
  statistics.
- Use an RMS envelope rather than samplewise absolute value.
- Choose envelope window length adaptively from eigenvalue-derived timescale.
- Preserve all existing P5 outputs; write a separate version with
  `_rmsenv_adaptive` suffix.

## Activity Transform Implemented

New helper:

```text
functions/postprocessing/compute_eigenfunction_activity.m
```

Activity transform used in this pass:

```text
activity(t, mode) = sqrt(movmean(abs(phi(t, mode)).^2, window_samples))
```

Window policy:

```text
activity_transform      = rms_envelope
activity_window_policy  = eigenvalue_adaptive_rms_envelope
envelope_policy         = eigenvalue_adaptive_rms
envelope_alpha          = 0.35
envelope_min_window_sec = 0.03
envelope_max_window_sec = 1.0
envelope_fallback_sec   = 0.10
```

For raw eigenfunctions, the window is derived from each mode eigenvalue.  In
continuous time, when the real part is negative:

```text
tau_sec = -dt / log(abs(lambda_discrete))
window_sec = clamp(0.35 * tau_sec, 0.03, 1.0)
```

This was corrected on 2026-05-25. Earlier notes and first-pass raw-timescale
plots used the bilinear continuous value `-1/real(lambda_continuous)` in some
places. The current preferred timescale for slow/fast-mode interpretation is
the discrete eigenvalue decay timescale above; bilinear values are legacy
provenance/fallback only. If a valid discrete eigenvalue/timescale cannot be
found, the fallback window is `0.10 sec`.

For dimred components, the component does not have its own direct Koopman
eigenvalue.  The current implementation estimates component timescale from the
source raw eigenfunction weights:

- source eigenvalue timescales are computed for raw modes;
- component weights are transformed with `abs`;
- weighted component timescale metadata are saved;
- the adaptive envelope window uses the component weighted median timescale.

The dimred metadata includes top contributing raw eigenfunction indices and
weighted timescale summaries where available.

## Code Changes

Core activity/density changes:

```text
functions/postprocessing/compute_eigenfunction_activity.m
functions/postprocessing/get_thresholded_density.m
functions/postprocessing/get_dimred_thresholded_density.m
```

P5 parameter propagation:

```text
functions/postprocessing/build_blp_eigenfunction_reduction_params.m
functions/postprocessing/build_blp_raw_eigenfunction_threshold_cfg.m
functions/postprocessing/build_blp_eigenfunction_reduction_cfg.m
functions/postprocessing/build_blp_eigenfunction_reduction_defaults.m
configs/cfg_eigenfunction_reduction_minimal.m
scripts/script_run_one_cfg_blp_eigenfunction_reduction.m
```

P5 peak-statistics activity support:

```text
functions/postprocessing/analyze_eigenfunction_component_peaks_by_consensus_state.m
scripts/script_run_current_pipeline5_peak_state_all_datasets.m
```

P8/P10 density source suffix support:

```text
functions/postprocessing/resolve_bold_cross_modal_density_sources.m
functions/postprocessing/build_bold_cross_modal_coupling_params.m
functions/postprocessing/build_bold_dimred_cross_modal_coupling_params.m
scripts/script_run_one_cfg_bold_cross_modal_coupling.m
scripts/script_run_one_cfg_bold_dimred_cross_modal_coupling.m
```

New runners/reports:

```text
scripts/script_run_current_p5_adaptive_envelope_density_unblocked.m
scripts/script_run_current_p5_adaptive_envelope_density_from_existing_reductions_unblocked.m
scripts/script_run_current_p5_adaptive_envelope_density_full_unblocked.m
scripts/script_run_current_p5_adaptive_envelope_raw_density_unblocked.m
scripts/script_run_current_p5_adaptive_envelope_dimred_density_unblocked.m
scripts/script_run_current_p5_adaptive_envelope_dimred_density_force_unblocked.m
scripts/script_run_current_p5_adaptive_envelope_peak_state_unblocked.m
scripts/script_run_current_p8_rmsenv_adaptive_xcorr_unblocked.m
scripts/script_run_current_p10_rmsenv_adaptive_xcorr_unblocked.m
scripts/report_current_p8_p10_top_xcorr_density.py
```

Summary-script updates:

```text
scripts/summarize_pipeline8_cross_session_consistency.py
scripts/summarize_pipeline10_cross_session_consistency.py
```

## Dataset/Condition Scope

Datasets used in this run:

```text
e10gb1
e10fV1
e10gh1
e10gw1
f12m01
```

New P5 condition tags:

```text
abs_projected_vlambda_rmsenv_adaptive
complex_split_projected_vlambda_rmsenv_adaptive
```

Important exclusion:

```text
e10fV1 | complex_split_projected_vlambda
```

This condition was already downgraded as invalid/stale because previous P6
top-window QC found selected top-30 eigenfunctions were time-constant in many
top windows.  It was skipped in P5 density, peak statistics, P8, and P10.

## P5 Density Execution

Initial mistake:

- I first started a full P5 adaptive rerun that began recomputing the
  dimension-reduction artifact itself under an `_rmsenv_adaptive` condition
  folder.
- That was not the intended final design for this pass.  In this pass,
  `rmsenv_adaptive` means: keep the existing SVD/NMF/MDS/UMAP component
  definitions, then recompute activity/density/peak statistics from the
  magnitude RMS envelope.  Under that design, the component definitions do not
  need to be refit.
- The mistaken route did write one orphan reduction MAT before I stopped it:

```text
E:\DataPons_processed\e10gb1\pipeline5_eigenfunction_reduction\abs_projected_vlambda_rmsenv_adaptive\svd_k03\mat\e10gb1_abs_projected_vlambda_rmsenv_adaptive_efun__time__abs__svd__svd_k03__20260523_102045.mat
```

- Size: about `317 MB`.
- This file is not part of the intended final density-only reuse strategy and
  should be treated as orphan/stale provenance unless we explicitly decide to
  create a new branch where the dimension reductions themselves are fit on RMS
  envelope activity.
- I did not delete it.
- After stopping that route, I switched to density-only reuse of the existing
  unsuffixed P5 reduction MAT files.

Important design distinction:

- If envelope is treated only as an activity transform for density, peak
  statistics, selectivity, P8, and P10, then existing dimred components can be
  reused.
- If envelope is treated as a new input feature space for fitting SVD/NMF/MDS/
  UMAP, then the dimension reductions do change and the full P5 reduction step
  must be rerun, producing another large set of MAT files.

Final P5 density strategy:

- Raw density: recomputed from current BLP EDMD outputs with
  `rms_envelope`.
- Dimred density: reused existing unsuffixed P5 reduction results and recomputed
  only thresholded density with `rms_envelope`.
- Existing P5 results were preserved.
- New outputs use `_rmsenv_adaptive`.

Final raw density runner:

```text
scripts/script_run_current_p5_adaptive_envelope_raw_density_unblocked.m
```

Final dimred density runner:

```text
scripts/script_run_current_p5_adaptive_envelope_dimred_density_unblocked.m
```

Then a forced short-file-name dimred refresh was run:

```text
scripts/script_run_current_p5_adaptive_envelope_dimred_density_force_unblocked.m
```

Reason for forced refresh:

- One generated dimred MAT path was exactly 260 characters:

```text
E:\DataPons_processed\e10gb1\pipeline5_dimred_thresholded_density\complex_split_projected_vlambda_rmsenv_adaptive\umap_k03\mat\e10gb1_dimred_thresholded_density__quantile__ratio_070__complex_split_projected_vlambda_rmsenv_adaptive_umap_k03__20260523_110054.mat
```

- MATLAB could not reliably `load`/`save` this long path on Windows.
- I changed adaptive dimred density save stems from long
  `*_dimred_thresholded_density*` to short `*_dtd*`.
- Old long-name files were not deleted.
- New short-name files are newer, so the resolver picks them.
- Example short file was verified with MATLAB `load(..., 'D')`.

P5 density output roots:

```text
E:\DataPons_processed\<dataset>\pipeline5_raw_thresholded_density\<condition>_rmsenv_adaptive\mat\...
E:\DataPons_processed\<dataset>\pipeline5_dimred_thresholded_density\<condition>_rmsenv_adaptive\<method_k>\mat\...
```

P5 density/peak completion counts:

| dataset | condition | raw MAT | short dimred MAT | peak stats |
|---|---:|---:|---:|---:|
| e10gb1 | abs_projected_vlambda_rmsenv_adaptive | 3 | 24 | 24 |
| e10gb1 | complex_split_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| e10fV1 | abs_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| e10fV1 | complex_split_projected_vlambda_rmsenv_adaptive | 0 | 0 | 0 |
| e10gh1 | abs_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| e10gh1 | complex_split_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| e10gw1 | abs_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| e10gw1 | complex_split_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| f12m01 | abs_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |
| f12m01 | complex_split_projected_vlambda_rmsenv_adaptive | 1 | 24 | 24 |

Notes:

- `e10gb1 abs` has 3 raw MAT files because of early smoke/partial attempts plus
  the final output.  Downstream resolver uses the latest matching file.
- The short dimred MAT count is the important one after the path-length fix.

## P5 Peak Statistics and Selectivity

Peak statistics runner:

```text
scripts/script_run_current_p5_adaptive_envelope_peak_state_unblocked.m
```

Important implementation detail:

- Source reduction variant remains unsuffixed:

```text
abs_projected_vlambda
complex_split_projected_vlambda
```

- Output peak-stat variant is suffixed:

```text
abs_projected_vlambda_rmsenv_adaptive
complex_split_projected_vlambda_rmsenv_adaptive
```

- This avoids recomputing reduction while still producing envelope-based peak
  statistics.

Peak output root:

```text
E:\DataPons_processed\<dataset>\pipeline5_eigenfunction_peaks_by_state_rmsenv_adaptive\<condition>_rmsenv_adaptive\<method_k>\...
```

Peak run summary:

```text
Discovered tasks after invalid exclusion: 216
OK: 178
Skipped existing: 38
Failed: 0
```

Selectivity summary command output:

```text
Runs found          : 216
Missing stats files : 24
Component rows      : 3564
```

The 24 missing files are the intentionally excluded `e10fV1 complex_split`
condition.

Selectivity outputs:

```text
results\peak_event_family_component_activity_current_rmsenv_adaptive\
  best_component_by_run_family.csv
  component_family_activity.csv
  method_family_summary.csv
  missing_stats_files.csv
  summary.md
```

Dimred process-label outputs:

```text
results\pipeline5_dimred_component_process_labels_current_rmsenv_adaptive\
```

Label command output:

```text
Input rows : 648
Label rows : 556
Missing    : 240
```

The label table is still based on the existing selectivity heuristic, now fed
with envelope-based peak statistics.  It is not yet the final subprocess-label
model.

## P8 Adaptive Xcorr

P8 runner:

```text
scripts/script_run_current_p8_rmsenv_adaptive_xcorr_unblocked.m
```

P8 settings:

```text
blp_density_condition_suffix = rmsenv_adaptive
xcorr_save_tag               = xcorr_rmsenv_adaptive
current_best_p7_only          = true
output_mode                  = separate
top_n                        = 5
make_xcorr_figures           = false
make_activation_maps         = false
make_roi_summaries           = false
require_all_density_sources  = false
```

P8 was limited to the current 4 main BOLD observable modes:

```text
HP_svd100
global_svd100
gsvd100_ds
roi_mean
```

P8 density sources:

- event density: `blp_evt`
- raw abs density: `raw_abs_q070_rmsenv_adaptive`
- dimred abs density:
  `dim_abs_<method><k>_q070_rmsenv_adaptive`, method = SVD/NMF/MDS/UMAP,
  k = 3:8
- raw complex-split density: `raw_csplit_q070_rmsenv_adaptive`
- dimred complex-split density:
  `dim_csplit_<method><k>_q070_rmsenv_adaptive`, method = SVD/NMF/MDS/UMAP,
  k = 3:8

For normal datasets this gives 51 sources:

```text
1 event + 1 raw_abs + 24 dim_abs + 1 raw_csplit + 24 dim_csplit = 51
```

For `e10fV1`, complex-split sources were skipped, giving 26 sources:

```text
1 event + 1 raw_abs + 24 dim_abs = 26
```

P8 output root:

```text
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\xcorr_rmsenv_adaptive.mat
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\density\xcorr_rmsenv_adaptive_top__*.csv
E:\DataPons_processed\<dataset>\pipeline8_xcorr\<run_tag>\feature\...\xcorr_rmsenv_adaptive_top__*.csv
```

P8 completion:

```text
p8_top_csv = 4600
```

Expected count:

```text
e10gb1  : 4 runs * 255 CSV = 1020
e10fV1  : 4 runs * 130 CSV = 520
e10gh1  : 4 runs * 255 CSV = 1020
e10gw1  : 4 runs * 255 CSV = 1020
f12m01  : 4 runs * 255 CSV = 1020
total   : 4600
```

P8 cross-session summary:

```text
results\pipeline8_cross_session_consistency_current_rmsenv_adaptive\
```

P8 summary command output:

```text
Score rows       : 5327
Readable hit rows: 23000
P8 topN          : 5
Figures          : 80
```

P8 figures were also written to:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p8\cross_session_consistency_rmsenv_adaptive\
```

## P10 Adaptive Xcorr

P10 runner:

```text
scripts/script_run_current_p10_rmsenv_adaptive_xcorr_unblocked.m
```

P10 settings:

```text
blp_density_condition_suffix = rmsenv_adaptive
xcorr_save_tag               = dimred_xcorr_rmsenv_adaptive
current_best_p7_only          = true
output_mode                  = separate
top_n                        = 5
make_xcorr_figures           = false
make_activation_maps         = false
make_roi_summaries           = false
require_all_density_sources  = false
```

P10 was also limited to the 4 main BOLD observable modes:

```text
HP_svd100
global_svd100
gsvd100_ds
roi_mean
```

After realizing the default P10 scope was too wide, I changed P10 to a
first-pass subset:

```text
P9 feature_names = efun_real, deconv_real
P9 method_tags   = svd_k05, svd_k08,
                   nmf_k05, nmf_k08,
                   mds_k05, mds_k08,
                   umap_k05, umap_k08
```

Important caveat:

- Before this narrowing, a broader partial P10 run had already written some
  `dimred_xcorr_rmsenv_adaptive` outputs.
- Because the broad partial run and first-pass run used the same save tag, the
  raw P10 output tree contains a mixture of first-pass and partial broad files.
- The clean first-pass top-density report below filters P10 to
  `efun_real/deconv_real` and the k05/k08 method tags.
- The P10 cross-session summary script currently does not have a method-tag
  filter, so its summary directory may include extra partial P10 artifacts if
  matching files exist.  Treat it as provisional until we either clean the tag
  or add method-tag filtering there too.

P10 output root:

```text
E:\DataPons_processed\<dataset>\pipeline10_dimred_xcorr\<p10_tag>\dimred_xcorr_rmsenv_adaptive.mat
E:\DataPons_processed\<dataset>\pipeline10_dimred_xcorr\<p10_tag>\density\dimred_xcorr_rmsenv_adaptive_top__*.csv
E:\DataPons_processed\<dataset>\pipeline10_dimred_xcorr\<p10_tag>\feature\...\dimred_xcorr_rmsenv_adaptive_top__*.csv
```

P10 top CSV count in output tree:

```text
p10_top_csv = 25948
```

This count includes the same-tag partial broad files mentioned above.

P10 cross-session summary:

```text
results\pipeline10_cross_session_consistency_current_rmsenv_adaptive\
```

P10 summary command output:

```text
Score rows       : 24112
Readable hit rows: 120560
P10 topN         : 5
Figures          : 26
```

Because of the same-tag partial-run caveat, treat this summary as provisional.

P10 figures were also written to:

```text
E:\DataPons_processed\summary_figures\pipeline11_current_analysis_summary\p10\cross_session_consistency_rmsenv_adaptive\
```

## Top-Xcorr Density Reports

Initial broad report:

```text
results\p8_p10_top_xcorr_density_current_rmsenv_adaptive\
```

Clean first-pass report:

```text
results\p8_p10_top_xcorr_density_current_rmsenv_adaptive_firstpass_clean\
```

The clean report applies P10 filters:

```text
p10_features    = efun_real, deconv_real
p10_method_tags = svd_k05, svd_k08,
                  nmf_k05, nmf_k08,
                  mds_k05, mds_k08,
                  umap_k05, umap_k08
```

Main file to inspect:

```text
results\p8_p10_top_xcorr_density_current_rmsenv_adaptive_firstpass_clean\p8_p10_top_xcorr_density_hits.csv
```

Top rows from the clean report:

```text
P8  e10gw1 pv_roi     raw_csplit_q070_rmsenv_adaptive mode150 r=0.8686 lag=2
P10 e10gb1 pv_gsvd100 raw_csplit_q070_rmsenv_adaptive mode147 r=0.8659 lag=8
P8  e10gw1 pv_roi     raw_csplit_q070_rmsenv_adaptive mode150 r=0.8620 lag=10
P10 e10gb1 pv_gsvd100 raw_csplit_q070_rmsenv_adaptive mode147 r=0.8611 lag=8
P10 e10gw1 pv_hp100   raw_csplit_q070_rmsenv_adaptive mode150 r=0.8567 lag=10
P8  e10gb1 pv_gsvd100 raw_csplit_q070_rmsenv_adaptive mode147 r=0.8456 lag=8
P10 e10gb1 pv_gsvd100 raw_csplit_q070_rmsenv_adaptive mode147 r=0.8402 lag=10
P8  e10gb1 pv_roi     raw_abs_q070_rmsenv_adaptive    mode162 r=0.8374 lag=6
```

This says only that the strongest individual top rows are raw-density hits,
especially raw complex-split modes.  It does not answer selectivity,
cross-session stability, method choice, or whether dimred density gives a more
interpretable subprocess label.

## Validation Performed

Static/smoke checks:

- MATLAB `checkcode` was run on new/modified runners.
- Python AST checks were run on modified summarizers/report scripts.
- A MATLAB toy smoke test passed for raw and dimred adaptive envelope density.
- A real e10gb1/svd_k03 smoke run passed before the broader run.
- The short-name `e10gb1 complex_split umap_k03` MAT was explicitly loaded in
  MATLAB after the path-length fix.

Run logs:

```text
results\logs\p5_rmsenv_adaptive_raw_only_20260523_105412.out.log
results\logs\p5_rmsenv_adaptive_dimred_force_short_20260523_121957.out.log
results\logs\p5_rmsenv_adaptive_peak_state_20260523_114125.out.log
results\logs\p8_rmsenv_adaptive_xcorr_main4_resume_20260523_131617.out.log
results\logs\p10_rmsenv_adaptive_xcorr_firstpass_20260523_131617.out.log
```

Aborted/partial logs also exist and should not be treated as final provenance:

```text
results\logs\p5_rmsenv_adaptive_unblocked_20260523_102214.*
results\logs\p5_rmsenv_adaptive_unblocked_20260523_102316.*
results\logs\p5_rmsenv_adaptive_unblocked_20260523_102851.*
results\logs\p5_rmsenv_adaptive_density_only_20260523_103630.*
results\logs\p5_rmsenv_adaptive_dimred_only_20260523_105412.*
results\logs\p5_rmsenv_adaptive_peak_state_20260523_113338.*
results\logs\p8_rmsenv_adaptive_xcorr_20260523_120108.*
results\logs\p10_rmsenv_adaptive_xcorr_20260523_120108.*
results\logs\p8_rmsenv_adaptive_xcorr_main4_20260523_124438.*
results\logs\p10_rmsenv_adaptive_xcorr_main4_20260523_124438.*
```

## Known Problems / Things To Fix Before Treating This As Mainline

1. The top-hit table is not the analysis you actually want.

   It is dominated by individual raw modes.  It does not answer cross-session
   interpretability or subprocess selectivity.

2. P10 output tag was reused during partial broad and first-pass runs.

   Clean first-pass top report is filtered, but the raw P10 output tree and
   provisional P10 summary may contain same-tag partial broad artifacts.  A
   cleaner rerun should use a distinct tag, for example:

   ```text
   dimred_xcorr_rmsenv_adaptive_firstpass_k05k08
   ```

3. P8/P10 still report strongest raw eigenfunction density, but we have not yet
   connected those raw indices to timescale/selectivity labels.

   This is required before deciding whether raw hits are meaningful or merely
   fast/noisy mode dominance.

4. The dimred component process labels exist, but they are still based on the
   older selectivity heuristic.

   They now use envelope peak statistics as input, but the label logic itself
   has not been redesigned.

5. No file deletion was performed.

   Old long-name adaptive density files and aborted/partial xcorr files remain
   in place.  This was deliberate, but it means audit scripts must prefer
   latest valid/short files or use explicit tags.

6. P8/P10 summary scripts need tighter filtering/provenance.

   P8 has `--xcorr-save-tag`; P10 now has `--xcorr-save-tag`, but P10 still
   needs method-tag filtering for clean first-pass summary generation.

7. This pass used five datasets, not the broader nine-dataset target.

   The current run matched the unblocked/current P5/P7/P8/P10 set:

   ```text
   e10gb1, e10fV1, e10gh1, e10gw1, f12m01
   ```

8. P8/P10 were run with figures/activation/ROI summaries disabled.

   Numeric xcorr tables and summary heatmaps were generated, but no activation
   maps or ROI maps were regenerated in this pass.

## What This Run Is Actually Good For

This run is useful for checking:

- whether adaptive envelope density can be generated end to end;
- whether P8/P10 can consume all `abs/csplit x method x k03:k08` P5 density
  sources;
- whether path-length and invalid-condition handling are now visible;
- which density sources dominate naive top-xcorr rows;
- whether P5 envelope peak/selectivity tables can be generated.

It is not sufficient for deciding:

- best `abs` vs `complex_split` observable;
- best dimred method;
- best component number;
- best BOLD observable type;
- whether a dimred component has stable theta/ripple selectivity;
- whether the strongest xcorr source is biologically meaningful.

Those require a redesigned analysis view, probably centered on:

- raw efun index/timescale distribution among top xcorr hits;
- dimred component process labels among top xcorr hits;
- cross-session agreement of density family/method/k;
- separation of P8 BOLD efun vs deconv_efun;
- separate raw-density and dimred-density conclusions.
