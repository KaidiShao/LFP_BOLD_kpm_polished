# P6 BLP Standardized SPKT/MUA Comparison

Date: 2026-05-26

This note records the Pipeline 6 probe comparing non-standardized and
standardized BLP ResKoopNet outputs on residual coupling to SPKT and MUA.

## Purpose

After standardized P4 BLP training appeared smoother and more stable, we tested
whether the downstream P6 residual-channel coupling changed in a meaningful
way. The immediate question was:

- Does standardized `abs` or `complex_split` change SPKT/MUA residual xcorr?
- Is any change a single pooled peak, or a broader session-level effect?

## Scope

Datasets:

- `e10gb1`
- `e10gh1`

Observables:

- `abs`
- `complex_split`

Variants:

- raw/non-standardized current P4 outputs
- standardized P4 outputs:
  - `abs`: `stdT`
  - `complex_split`: `std_complex_pair`

P6 stages run:

- `pipeline6_spkt_residual_cross_correlation`
- `pipeline6_mua_residual_cross_correlation`

P6 stages intentionally skipped for this probe:

- top-window postprocessing figures
- timescale diagnostics
- deconvolution window plots

This made the run focused on residual-channel coupling rather than visual
eigenfunction inspection.

## Inputs

Standardized full EDMD chunks were read from:

```text
E:\DataPons_processed\derived_autodl_results_standardize\<dataset>\mlp\outputs\
```

At the start of this probe, three standardized full-chunk exports already
existed:

- `e10gb1 abs stdT`: 600 chunks
- `e10gb1 complex_split std_complex_pair`: 600 chunks
- `e10gh1 complex_split std_complex_pair`: 882 chunks

The missing input was exported before P6:

- `e10gh1 abs stdT`: 882 chunks

Export helper updated:

```text
scripts/export_e10gh1_standardized_edmd_outputs.py
```

The update allowed `--run abs_std` in addition to `complex_std`.

## P6 Runner

New focused P6 runner:

```text
scripts/script_run_standardized_blp_p6_spkt_mua_probe_20260525.m
```

It sets:

```matlab
autodl_root = 'E:\DataPons_processed\derived_autodl_results_standardize';
make_main_plot = false;
make_timescale_plot = false;
make_deconv_plot = false;
make_deconv_window_norm_plot = false;
make_spkt_residual_cross_correlation = true;
make_mua_residual_cross_correlation = true;
spkt_cross_skip_existing = false;
mua_cross_skip_existing = false;
```

Temporary PowerShell driver used for the background MATLAB run:

```text
tmp/run_p6_std_spkt_mua_20260525.ps1
```

Run log:

```text
tmp/run_logs/p6_std_spkt_mua_20260525.out.log
tmp/run_logs/p6_std_spkt_mua_20260525.err.log
tmp/run_logs/p6_std_spkt_mua_20260525.status.log
```

The MATLAB P6 run finished successfully:

```text
END 2026-05-26T00:14:26.1240479+08:00 EXIT=0
```

## Output Locations

Canonical P6 outputs were written under:

```text
E:\DataPons_processed\e10gb1\pipeline6_spkt_residual_cross_correlation\
E:\DataPons_processed\e10gb1\pipeline6_mua_residual_cross_correlation\
E:\DataPons_processed\e10gh1\pipeline6_spkt_residual_cross_correlation\
E:\DataPons_processed\e10gh1\pipeline6_mua_residual_cross_correlation\
```

Summary and figures are here:

```text
results\pipeline6_std_vs_raw_spkt_mua_20260525\
```

Important files:

```text
results\pipeline6_std_vs_raw_spkt_mua_20260525\README.md
results\pipeline6_std_vs_raw_spkt_mua_20260525\std_vs_raw_delta.csv
results\pipeline6_std_vs_raw_spkt_mua_20260525\summary_by_variant.csv
results\pipeline6_std_vs_raw_spkt_mua_20260525\session_feature_top_channels.csv
results\pipeline6_std_vs_raw_spkt_mua_20260525\pooled_top_channels.csv
results\pipeline6_std_vs_raw_spkt_mua_20260525\figures\mean_session_top_abs_corr.png
results\pipeline6_std_vs_raw_spkt_mua_20260525\figures\metrics_raw_vs_standardized.png
```

Summary script:

```text
scripts/summarize_pipeline6_std_vs_raw_spkt_mua.py
```

Note: standardized `complex_split` P6 save tags exceeded Windows path length in
the MATLAB saver, so it fell back to compact names such as
`e10gb1_spkt_spkt_all_top20_*`. The summary script explicitly maps those compact
files back to `complex_split standardized` for this probe.

## Metrics

For session-level metrics, the script first takes the row with maximum
`abs_corr` within each:

```text
dataset x modality x observable x variant x session x residual feature
```

It then summarizes those session-feature top correlations.

Metrics reported:

- `mean session-top |corr|`
- `median session-top |corr|`
- `max session-top |corr|`
- `pooled top |corr|`, defined as the maximum `abs_corr` in
  `*_pooled_top_xcorr.csv`

## Results

Mean session-top absolute correlation:

| dataset | modality | observable | raw | standardized | std/raw |
|---|---|---:|---:|---:|---:|
| e10gb1 | MUA | abs | 0.01861 | 0.01932 | 1.04 |
| e10gb1 | MUA | complex_split | 0.02235 | 0.02388 | 1.07 |
| e10gb1 | SPKT | abs | 0.05221 | 0.05263 | 1.01 |
| e10gb1 | SPKT | complex_split | 0.07252 | 0.07016 | 0.967 |
| e10gh1 | MUA | abs | 0.03382 | 0.04159 | 1.23 |
| e10gh1 | MUA | complex_split | 0.04570 | 0.04596 | 1.01 |
| e10gh1 | SPKT | abs | 0.06901 | 0.09253 | 1.34 |
| e10gh1 | SPKT | complex_split | 0.12175 | 0.11984 | 0.984 |

Median, max session-top, and pooled top:

| dataset | modality | observable | metric | raw | standardized | std/raw |
|---|---|---|---:|---:|---:|---:|
| e10gb1 | MUA | abs | median | 0.01722 | 0.01799 | 1.04 |
| e10gb1 | MUA | abs | max session-top | 0.04370 | 0.04377 | 1.00 |
| e10gb1 | MUA | abs | pooled top | 0.01982 | 0.02155 | 1.09 |
| e10gb1 | MUA | complex_split | median | 0.02234 | 0.02414 | 1.08 |
| e10gb1 | MUA | complex_split | max session-top | 0.04538 | 0.04543 | 1.00 |
| e10gb1 | MUA | complex_split | pooled top | 0.02329 | 0.02180 | 0.936 |
| e10gb1 | SPKT | abs | median | 0.05104 | 0.05411 | 1.06 |
| e10gb1 | SPKT | abs | max session-top | 0.10850 | 0.10920 | 1.01 |
| e10gb1 | SPKT | abs | pooled top | 0.06007 | 0.06306 | 1.05 |
| e10gb1 | SPKT | complex_split | median | 0.07273 | 0.07167 | 0.985 |
| e10gb1 | SPKT | complex_split | max session-top | 0.12090 | 0.12970 | 1.07 |
| e10gb1 | SPKT | complex_split | pooled top | 0.06366 | 0.06305 | 0.990 |
| e10gh1 | MUA | abs | median | 0.02641 | 0.03474 | 1.32 |
| e10gh1 | MUA | abs | max session-top | 0.11350 | 0.12330 | 1.09 |
| e10gh1 | MUA | abs | pooled top | 0.03847 | 0.04085 | 1.06 |
| e10gh1 | MUA | complex_split | median | 0.04330 | 0.04397 | 1.02 |
| e10gh1 | MUA | complex_split | max session-top | 0.12500 | 0.12760 | 1.02 |
| e10gh1 | MUA | complex_split | pooled top | 0.04106 | 0.04292 | 1.05 |
| e10gh1 | SPKT | abs | median | 0.06462 | 0.09040 | 1.40 |
| e10gh1 | SPKT | abs | max session-top | 0.19570 | 0.20970 | 1.07 |
| e10gh1 | SPKT | abs | pooled top | 0.09935 | 0.09880 | 0.994 |
| e10gh1 | SPKT | complex_split | median | 0.11197 | 0.11124 | 0.994 |
| e10gh1 | SPKT | complex_split | max session-top | 0.21200 | 0.21210 | 1.00 |
| e10gh1 | SPKT | complex_split | pooled top | 0.10030 | 0.10030 | 1.00 |

## Interpretation

The clearest effect is in `e10gh1 abs standardized`.

- SPKT mean session-top increased from 0.06901 to 0.09253.
- SPKT median session-top increased from 0.06462 to 0.09040.
- MUA mean and median also increased.
- SPKT pooled top did not increase: 0.09935 to 0.09880.

This suggests the standardized `abs` branch is not merely creating one larger
pooled peak. Instead, the coupling appears stronger across more sessions or
session-feature slices.

For `complex_split`, the standardized branch is broadly similar to raw.

- `e10gb1 SPKT complex_split` mean decreased slightly, but max session-top
  increased.
- `e10gh1 SPKT complex_split` is nearly unchanged across median, max, and pooled
  top.
- MUA complex_split changes are small.

The recurrent pooled top channels are stable:

- `e10gb1`: mostly `ch10_hp`
- `e10gh1`: mostly `ch11_pl`

So standardization changes the distribution of coupling strength more than it
changes the leading pooled channel identity.

## Reproduction

Run P6 standardized SPKT/MUA:

```powershell
powershell -NoProfile -ExecutionPolicy Bypass -File tmp\run_p6_std_spkt_mua_20260525.ps1
```

Regenerate summary and figures from Windows Python, with SVG fallback if
matplotlib is missing:

```powershell
python scripts\summarize_pipeline6_std_vs_raw_spkt_mua.py --datasets e10gb1 e10gh1 --output-dir results\pipeline6_std_vs_raw_spkt_mua_20260525
```

Regenerate PNG figures with WSL matplotlib:

```powershell
wsl.exe -d Ubuntu-22.04 bash -lc "cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished && python3 scripts/summarize_pipeline6_std_vs_raw_spkt_mua.py --processed-root /mnt/e/DataPons_processed --datasets e10gb1 e10gh1 --output-dir results/pipeline6_std_vs_raw_spkt_mua_20260525"
```

