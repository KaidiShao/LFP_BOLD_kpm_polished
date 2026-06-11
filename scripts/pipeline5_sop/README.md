# Pipeline 5 SOP: Band-Selectivity and Theta/RG Density Checks

This folder contains the formal, minimal P5 SOP code for standardized
complex-split BLP dimred eigenfunction analysis.

The default SOP question is:

```text
Can a dataset's P5 dimred eigenfunction space separate two LFP/BLP subprocesses?

1. theta-like slow/state subprocess
2. ripple/gamma-no-theta fast subprocess
```

The default scope is intentionally narrow:

```text
condition = complex_split_projected_vlambda_standardize
transform = adaptive_envelope
k = 03:16
labels = strict P2 theta/gamma/ripple band-event selectivity
```

The legacy `abs` transform is still available for sensitivity checks, but it is
not the default SOP.

## One-Command Runner

Run from WSL/Anaconda or any Python environment with `numpy`, `pandas`, `h5py`,
and `matplotlib`:

```bash
python scripts/pipeline5_sop/run_p5_sop.py \
  --datasets e10gb1 f12m05 \
  --workspace /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished \
  --processed-root /mnt/e/DataPons_processed
```

Optional cross-dataset theta/RG density-correlation refresh:

```bash
python scripts/pipeline5_sop/run_p5_sop.py \
  --datasets e10gb1 e10fV1 e10gh1 e10gw1 f12m01 f12m02 f12m03 f12m05 k13m18 k13m23 \
  --refresh-density-correlation
```

Dry-run only:

```bash
python scripts/pipeline5_sop/run_p5_sop.py --datasets e10gb1 --dry-run
```

## Default Outputs

Per dataset:

```text
results/<dataset>_p5_p2_band_event_response_v2_selective_envelope_20260528/
    component_band_event_response_v2.csv
    README_p5_band_event_response_v2.md
    figures/
        02_adaptive_envelope_strict_label_grid_by_method_k.png
        03_adaptive_envelope_strict_label_composition_by_method_k.png
        04_adaptive_envelope_strict_theta_vs_ripple_gamma_effect_scatter.png
        05_adaptive_envelope_strict_two_subprocess_candidate_map.png
```

Diagnostic plates for method-k settings that pass the two-subprocess gate:

```text
E:/DataPons_processed/<dataset>/
    pipeline5_event_diagnostic_plates/
        complex_split_projected_vlambda_standardize/
            adaptive_envelope/
                <dataset>_adaptive_envelope_<method>_k<kk>_thetaC<cc>_rgC<cc>_sem_diagnostic_plate.png
                diagnostic_plate_manifest_sem.csv
```

The diagnostic plate is generated only for method-k settings with both:

```text
theta_selective component
RG-like component:
    ripple_gamma_no_theta, gamma_selective, or ripple_selective
```

Cross-dataset theta/RG density correlation, if requested:

```text
results/standardized_csplit_k03_k16_all_current_20260607/
    p5_theta_rg_density_correlation/
        theta_rg_density_all_pair_correlations.csv
        theta_rg_density_canonical_pair_correlations.csv
        theta_rg_density_correlation_dataset_summary.csv
        sop_figures/
            01_main_dataset_methodk_canonical_theta_rg_density_corr_heatmap__*.png
```

Optional distribution plots are available through:

```bash
python scripts/pipeline5_sop/plot_theta_rg_density_correlation.py --include-supplements
```

## What Is Not Default

These are intentionally excluded from the default P5 SOP:

```text
abs transform
standalone selected_component_perievent_mean_by_band figure
paired top-window sheets for every method-k
P2/P5 top30 window QC
two-subprocess consensus-state trajectory
```

Use them only as sensitivity/QC/showcase outputs after the minimal SOP has
identified the method-k settings worth inspecting.

## Scripts

```text
run_p5_sop.py
    One-command orchestrator.

analyze_band_event_response.py
    Computes strict P2-band event labels for each P5 dimred component and
    writes the four default P5 gate figures.

plot_event_diagnostic_plate.py
    Plots selected theta/RG component peri-event means plus representative top
    windows. Default runner calls it only for pair-pass method-k settings.

analyze_theta_rg_density_correlation.py
    Writes theta/RG density-correlation CSVs. By default it does not write the
    old exploratory summary_figures.

plot_theta_rg_density_correlation.py
    Writes the SOP theta/RG density-correlation heatmap. Supplement
    distribution plots require --include-supplements.
```

## Minimality Check

The formal runner only calls:

```text
P5 band-event strict label analysis
P5 pair-pass diagnostic plates
optional theta/RG density-correlation tables and SOP figures
```

It does not rerun P5 reductions, P5 density backfill, P8, P10, or any BOLD
analysis. Those must already exist.
