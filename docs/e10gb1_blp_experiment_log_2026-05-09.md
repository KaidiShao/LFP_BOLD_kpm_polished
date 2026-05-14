# E10gb1 BLP Experiment Log

Date: `2026-05-09`

This note records the main `E10gb1` BLP experiments carried out during the
current ResKoopNet / projected-`vlambda` investigation. The main goal was to
identify a solver/training setup that:

- works well on `E10gb1` BLP data,
- stays reasonably close to the official ResKoopNet training flavor,
- and still preserves the key original-paper validation behavior on
  `pendulum` and `turbulence`.

## Data And Validation Context

- BLP dataset: `e10gb1_low50_high250_g2_abs_single.mat`
- Official pendulum data now confirmed to match:
  `tmp/figure2d_lab/original_data/data_pendulum_90.mat`
- Official turbulence input now confirmed to match:
  `tmp/turbulence_lab/original_data/pressure_data.mat`

Related reference notes:

- `docs/reskoopnet_reference_data_manifest_2026-05-08.md`
- `docs/reskoopnet_solver_versions_2026-05-08.md`

## 1. Initial Seed Sweep

The first pass was a seed sweep on the earlier `projected_vlambda` setup.

Summary:

| Seed | Best epoch | Best val metric |
| --- | ---: | ---: |
| `1234` | `8` | `0.0415676` |
| `2026` | `2` | `0.0614613` |
| `3407` | `11` | `0.0479781` |
| `7777` | `7` | `0.0280172` |
| `9999` | `2` | `0.0738381` |

Takeaway:

- The setup was clearly seed-sensitive.
- `seed 7777` was the best run from this early sweep.
- At that point, the configuration still did not look robust enough for
  large-scale production runs.

Artifacts:

- `E:\autodl_results\e10gb1\mlp\seed_sweep_summary\e10gb1_abs_seed_sweep_outer_metrics.png`
- `E:\autodl_results\e10gb1\mlp\seed_sweep_summary\e10gb1_abs_seed_sweep_summary.csv`

## 2. Learning-Rate Sweep

The next step compared `lr=1e-5` and `lr=3e-5` at fixed
`reg=0.3`, `inner_epochs=2`.

Summary:

| LR | Seed | Best epoch | Best val | Final val |
| --- | --- | ---: | ---: | ---: |
| `1e-5` | `1234` | `6` | `0.0225012` | `0.0225012` |
| `1e-5` | `2026` | `6` | `0.0272500` | `0.0272500` |
| `1e-5` | `3407` | `6` | `0.0853635` | `0.0853635` |
| `3e-5` | `1234` | `5` | `0.0312899` | `0.0906753` |
| `3e-5` | `2026` | `6` | `0.0299661` | `0.0299661` |
| `3e-5` | `3407` | `6` | `0.1055652` | `0.1055652` |

Takeaway:

- `1e-5` was clearly better and more stable than `3e-5`.
- The old `1e-4` line was effectively abandoned after this stage.
- This established the new mainline LR candidate: `1e-5`.

Artifacts:

- `E:\autodl_results\e10gb1\mlp\hparam_sweep_summary\e10gb1_hparam_sweep_outer_metrics.png`
- `E:\autodl_results\e10gb1\mlp\hparam_sweep_summary\e10gb1_hparam_sweep_loss_curve_contact_sheet.png`
- `E:\autodl_results\e10gb1\mlp\hparam_sweep_summary\e10gb1_hparam_sweep_summary.csv`

## 3. Regularization Sweep At LR = 1e-5

With `lr=1e-5` fixed, the next sweep compared `reg=0.2`, `0.3`, and `0.5`
using two seeds.

Summary:

| Reg | Seed | Best epoch | Best val | Final val |
| --- | --- | ---: | ---: | ---: |
| `0.2` | `1234` | `5` | `0.0214383` | `0.0686170` |
| `0.2` | `2026` | `6` | `0.0377136` | `0.0377136` |
| `0.3` | `1234` | `6` | `0.0220203` | `0.0220203` |
| `0.3` | `2026` | `6` | `0.0260069` | `0.0260069` |
| `0.5` | `1234` | `6` | `0.0212876` | `0.0212876` |
| `0.5` | `2026` | `5` | `0.0386253` | `0.0519706` |

Takeaway:

- `reg=0.3` looked like the best balanced choice across seeds.
- `reg=0.2` could reach a good point, but one run bounced badly later.
- `reg=0.5` was strong on one seed but less stable overall.

Artifacts:

- `tmp/hparam_sweep_20260506_summary/e10gb1_reg_sweep_loss_comparison.png`
- `tmp/hparam_sweep_20260506_summary/e10gb1_reg_sweep_loss_comparison.csv`

## 4. Long-Term Solver / Schedule Comparison

Once `lr=1e-5` and `reg=0.3` were fixed, the focus shifted to long-term
training behavior and the effect of the solver schedule.

The main conditions compared were:

- `dual-sync dynamic, inner=2`
- `torch-like dynamic, inner=2`
- `torch-like dynamic, inner=4`
- `static target engineered, inner=2`
- `torch-like dynamic, inner=1`

Final comparison summary:

| Condition | Completed outer epochs | Best outer epoch | Best outer val | Final outer val | Min inner val per-dim | Final inner val per-dim |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Dual-sync dynamic, `inner=2` | `14` | `6` | `0.0227029` | `0.1063262` | `4.4303e-05` | `1.11196e-04` |
| Torch-like dynamic, `inner=2` | `30` | `8` | `0.0150615` | `0.0472528` | `3.6646e-05` | `1.14970e-04` |
| Torch-like dynamic, `inner=4` | `15` | `5` | `0.0189667` | `0.0352194` | `4.6148e-05` | `8.5692e-05` |
| Static target engineered, `inner=2` | `30` | `1` | `0.0489045` | `0.1092352` | `9.2126e-05` | `9.5700e-05` |
| Torch-like dynamic, `inner=1` | `30` | `12` | `0.0158882` | `0.0411072` | `3.8657e-05` | `1.00017e-04` |

Interpretation:

- `torch-like, inner=2` gave the best overall outer metric:
  `best outer val = 0.0150615 @ epoch 8`.
- `torch-like, inner=1` was surprisingly competitive, but still slightly worse
  on best outer value.
- `torch-like, inner=4` was not the deepest best point, but it ended with a
  lower final inner loss than the `inner=1/2` dynamic runs.
- The engineered static-target solver reproduced the old smooth inner-loss
  behavior, but did not give competitive dynamic outer metrics.

Main plot:

- `tmp/live_plots/e10gb1_condition_comparison_final_20260509.png`
- `tmp/live_plots/e10gb1_condition_comparison_final_20260509.csv`

The old local-TF-like feeling was most visible in the static-target run:

- inner loss descended smoothly,
- but the dynamic outer metric did not continue improving.

That supports the interpretation that the older local behavior was strongly
connected to a more static spectral target.

## 5. Static-Target Engineered Solver

A new engineered solver was added in order to emulate the older local TF flavor
while keeping the current batching / export infrastructure:

- `python_scripts/autodl/solver_resdmd_batch3_static_target.py`

Behavior:

- training target `V / Lambda` is frozen after the initial synchronization,
- inner training remains stable and smooth,
- but outer evaluation still uses the dynamically recomputed spectral state.

Purpose:

- isolate whether the BLP instability was coming from the inner optimizer,
- or from repeatedly changing the spectral target in the outer loop.

Result:

- inner optimization itself was healthy,
- the larger issue was the moving-target outer schedule.

## 6. Pendulum And Turbulence Reproducibility Checks

The solver choice was also tested against the original-paper validation tasks.

### Turbulence

For the torch-like branch, the most important `mode 1` remained intact.

Key result:

- `focus_mode_abs_corr = 0.9965`
- `focus_mode_corr_pass = true`
- `focus_mode_rank_by_residue_ours_1based = 1`

Artifact:

- `tmp/active_validation/outputs/turbulence_active_torchlike_20260508/comparison_vs_official_mode1check.json`

Interpretation:

- torch-like scheduling is acceptable for the key turbulence mode-1 criterion.

### Pendulum

The pendulum story was subtler:

- simple matrix-level comparisons looked different between branches,
- but when the correct official-drive pendulum data was used,
  the torch-like branch still produced a very similar pseudospectrum shape.

Artifacts:

- Torch-like official-data pseudospectrum:
  `tmp/active_validation/outputs/pendulum_active_torchlike_officialdrive_20260508/pseudospectrum/pseudospectrum_notebook_style.png`
- Mainline baseline pseudospectrum:
  `tmp/active_validation/outputs/pendulum_active_mainline/pseudospectrum/pseudospectrum_notebook_style.png`

Interpretation:

- the torch-like branch does not preserve every internal quantity exactly,
- but it still preserves the pendulum pseudospectrum structure well enough to
  count as a successful reproduction in the visual sense that matters here.

## 7. Current Recommendation

For the current `E10gb1` BLP mainline, the best working recommendation is:

- solver: polished TF solver in torch-like mode
- schedule: `spectral_sync_mode=pre_only`
- `inner_epochs = 2`
- `lr = 1e-5`
- `reg = 0.3`
- `seed = 1234`
- export results from the **best checkpoint**, not from the final epoch

Why:

- it gives the best outer metric among the tested conditions,
- it is the closest practical match to the official paper / Torch schedule,
- and it still preserves the key original-paper validation behavior on
  `turbulence`, while also reproducing the pendulum pseudospectrum using the
  official data.

## 8. Exported E10gb1 Best-Checkpoint Eigenfunction Outputs

The best `torch-like, inner=2` checkpoint was restored and exported with full
eigenfunction chunks.

Current location:

- `E:\autodl_results_new\e10gb1\mlp\outputs\mlp_obs_blp_vlambda_torchlike_longterm_full_20260507_e10gb1_seed1234_projected_vlambda_abs`

Notes:

- this export was generated from the best checkpoint,
- not from the final epoch,
- and currently contains the summary file plus full `*_outputs_*.mat` chunk
  files.

