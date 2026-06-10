# P4 Standardization And Mixed-GPU Solver Archive

Date: 2026-05-21

This note archives the current P4 standardization experiments and the solver
acceleration changes made around the BLP/BOLD ResKoopNet runs. The main goal is
to keep the experimental branches distinguishable from existing mainline and
legacy runs.

## Short Summary

- BOLD `standardize_data=true` runs were completed for `roi_mean` and
  `global_svd100` across the current BOLD set using the low-learning-rate
  branch `stdT_l1e4_r1e3_b2000_i2_20260521_allbold`.
- BLP `e10fV1 abs` was used as the GPU speed probe because it is large enough
  to expose the bottleneck.
- A separate solver, `resdmd_batch_mixedgpu`, now keeps the DATA tensor on GPU
  and gathers snapshot rows by index. This does not replace the main solver.
- The mixed-GPU path is much faster: one outer epoch dropped from about
  880 seconds for the earlier full-float64 GPU probe to about 90-105 seconds.
- A long `batch=5000` mixed-GPU BLP run is active as of this archive. At
  `2026-05-21 16:33 +08:00`, it reached epoch 9 with best validation loss
  `8.418481e-04`.
- A future complex observable branch should use
  `complex_split + std_complex_pair`: the solver remains real-valued, while
  real/imag spectrogram columns share a complex-pair scale during input
  standardization.

## BOLD Standardized P4 Runs

Run tag:

```text
stdT_l1e4_r1e3_b2000_i2_20260521_allbold
```

Common settings:

```text
standardize_data = true
lr               = 1e-4
reg              = 0.001
batch_size       = 2000
inner_epochs     = 2
seed             = 1234
residual_form    = projected_vlambda
export_mode      = summary_only
```

Loss curve figures:

- `tmp/monitor_reports/bold_p4_stdT_lowlr_global_svd100_loss_curves_20260521.png`
- `tmp/monitor_reports/bold_p4_stdT_lowlr_roi_mean_loss_curves_20260521.png`

Summary CSV:

- `tmp/monitor_reports/bold_p4_stdT_lowlr_loss_curves_20260521_summary.csv`

Results archived from that CSV:

| Dataset | Observable | Status | Epochs | Best Epoch | Best Val | Final Val | Stale |
|---|---|---:|---:|---:|---:|---:|---:|
| e10aw1 | global_svd100 | stale_done | 3000 | 548 | 0.036275 | 0.126092 | 2452 |
| e10bv1 | global_svd100 | stale_done | 3000 | 787 | 0.024152 | 0.067348 | 2213 |
| e10fV1 | global_svd100 | stale_done | 6500 | 4310 | 0.061269 | 0.262790 | 2190 |
| e10gb1 | global_svd100 | stale_done | 6500 | 4166 | 0.053479 | 0.110932 | 2334 |
| e10gh1 | global_svd100 | stale_done | 5500 | 3147 | 0.040502 | 0.070734 | 2353 |
| e10gw1 | global_svd100 | stale_done | 3500 | 1056 | 0.042981 | 0.125395 | 2444 |
| f12m01 | global_svd100 | stale_done | 3000 | 944 | 0.039681 | 0.101371 | 2056 |
| k13m17 | global_svd100 | stale_done | 3500 | 1321 | 0.062116 | 0.167337 | 2179 |
| k13m23 | global_svd100 | missing_observable | | | | | |
| e10aw1 | roi_mean | stale_done | 2500 | 116 | 0.005043 | 0.007984 | 2384 |
| e10bv1 | roi_mean | stale_done | 3000 | 960 | 0.006314 | 0.020717 | 2040 |
| e10fV1 | roi_mean | stale_done | 2500 | 323 | 0.014083 | 0.024033 | 2177 |
| e10gb1 | roi_mean | stale_done | 2500 | 414 | 0.007967 | 0.020729 | 2086 |
| e10gh1 | roi_mean | stale_done | 5500 | 3367 | 0.013297 | 0.044809 | 2133 |
| e10gw1 | roi_mean | stale_done | 3500 | 1235 | 0.010324 | 0.017422 | 2265 |
| f12m01 | roi_mean | stale_done | 2500 | 470 | 0.006817 | 0.010522 | 2030 |
| k13m17 | roi_mean | stale_done | 3000 | 664 | 0.005875 | 0.014786 | 2336 |
| k13m23 | roi_mean | stale_done | 3500 | 1146 | 0.010531 | 0.017202 | 2354 |

Interpretation at archive time:

- Standardization did not automatically make all BOLD loss curves smooth.
- Many runs reached the `best unchanged for about 2000 epochs` stopping rule,
  but several validation curves still bounced or degraded after the best point.
- These losses are in standardized coordinates and should not be compared
  directly to old raw-coordinate losses.

## BLP Standardized Speed Probes

Dataset and observable:

```text
dataset_stem    = e10fV1
observable_mode = abs
seed            = 1234
standardize     = zscore_by_observable
```

Key probe results:

| Run | Batch | Numeric Policy | Outer Epoch 1 Time | Outer Compute_K | Val | Spectrum |
|---|---:|---|---:|---:|---:|---|
| full float64 GPU baseline | 2000 | train/analysis/gram/spectral float64 | 879.99s | about 247s | 0.002665 | not exported |
| mixedgpu one-epoch probe | 5000 | train float32, gram/spectral float64 | 104.35s | 11.91s | 0.003049 | stable |
| mixedgpu one-epoch probe | 10000 | train float32, gram/spectral float64 | 89.92s | 11.35s | 0.003641 | stable |

Spectrum comparison figure:

- `tmp/monitor_reports/blp_e10fV1_abs_stdT_gpu_mixedgpu_i2_spectrum_compare_20260521.png`

Spectrum statistics:

| Batch | Best Val | Final Eigvec Cond | Max Abs Lambda | Median Abs Lambda | Count Abs Lambda > 1 |
|---:|---:|---:|---:|---:|---:|
| 5000 | 0.003048662 | 2.421e5 | 0.99999838 | 0.736479 | 0 |
| 10000 | 0.003640518 | 1.012e5 | 0.99999840 | 0.736449 | 0 |

Interpretation:

- Both one-epoch probes had all eigenvalues inside the unit circle.
- `batch=10000` was faster and had a lower final eigenvector condition number,
  but `batch=5000` had better validation loss after the first outer epoch.
- The next useful test is a longer `batch=5000` run.

## Active Long BLP Mixed-GPU Run

Script:

- `tmp/run_blp_stdT_gpu_mixedgpu_i2_b5000_e100_e10fV1_abs_20260521.sh`

Experiment:

```text
blp_vlambda_stdT_gpu_mixedgpu_i2_b5000_e100_abs_20260521_e10fV1_seed1234
```

Log:

- `tmp/run_logs/blp_vlambda_stdT_gpu_mixedgpu_i2_b5000_e100_abs_20260521_e10fV1_seed1234.log`

Settings:

```text
solver_name        = resdmd_batch_mixedgpu
training_policy    = float32
analysis_dtype     = float64
gram_dtype         = float64
spectral_dtype     = float64
batch_size         = 5000
inner_epochs       = 2
epochs             = 100
lr                 = 1e-4
reg                = 0.1
train_shuffle      = true
spectral_sync_mode = pre_only
export_mode        = summary_only
seed               = 1234
```

Progress at `2026-05-21 16:33 +08:00`:

| Epoch | Train | Val | Eigvec Cond | Outer Time |
|---:|---:|---:|---:|---:|
| 1 | 3.049029e-03 | 3.048662e-03 | 1.419e4 | 103.33s |
| 2 | 1.539321e-03 | 1.538637e-03 | 7.621e4 | 88.15s |
| 3 | 1.259063e-03 | 1.259084e-03 | 4.914e4 | 86.23s |
| 4 | 1.130569e-03 | 1.131080e-03 | 4.854e4 | 85.39s |
| 5 | 1.098430e-03 | 1.099406e-03 | 5.268e4 | 88.40s |
| 6 | 9.742625e-04 | 9.750038e-04 | 6.416e4 | 87.09s |
| 7 | 8.958632e-04 | 8.966114e-04 | 5.312e4 | 87.80s |
| 8 | 9.174928e-04 | 9.179100e-04 | 8.580e4 | 92.99s |
| 9 | 8.410273e-04 | 8.418481e-04 | 5.976e4 | 88.72s |

Important caveat:

- This long run was started before the later runner-side load-dtype optimization.
  It therefore loaded the CPU-side DATA matrix as float64, then cached the base
  tensor on GPU as float32. The result is still valid, but the CPU memory usage
  is higher than future mixed-GPU runs should be.

## Solver Acceleration Changes

New solver:

- `python_scripts/autodl/solver_resdmd_batch3_mixedgpu.py`

Solver registry:

- `python_scripts/autodl/edmd_utils.py`

Runner:

- `python_scripts/autodl/run_autodl_reskoopnet_mlp.py`

Key implementation points:

- `resdmd_batch_mixedgpu` is separate from the production solver.
- It detects `IndexedRowsView` snapshot pairs.
- It caches one base DATA tensor on GPU as `float32`.
- Training datasets are built from row-index tensors; x/y batches are produced
  with `tf.gather(base, row_idx)`.
- `compute_K` and `evaluate_spectral_losses` also gather rows from the cached
  base tensor instead of materializing NumPy batches.
- Training uses `float32`, while analysis/Gram/spectral operations can still
  cast to `float64`.
- Current speedup comes mostly from eliminating repeated NumPy materialization
  and host-to-device transfers during pair batching and `compute_K`.

Future mixed-GPU runner optimization:

- `run_autodl_reskoopnet_mlp.py` now resolves load dtype specially for
  `solver_name=resdmd_batch_mixedgpu` plus `training_policy=float32`.
- Future mixed-GPU runs will load observable DATA as float32 even when
  `analysis_dtype=float64`, because the solver performs analysis casts inside
  `compute_K`.
- This should reduce CPU memory pressure, especially for `complex_split`.

## Complex Split Standardization

New runner option:

```bash
--standardize-data --standardize-mode std_complex_pair
```

Intended use:

```text
observable_mode = complex_split
solver          = real-valued solver4 or mixedgpu solver
```

Definition:

- BLP raw columns are still z-scored independently.
- Spectrogram real/imag columns are paired using the neighboring
  `*_obs_info.csv`.
- For each complex pair, real and imaginary columns keep their own means but
  share one scale:

```text
z_re_std = (z_re - mean_t(z_re)) / sqrt(std_t(z_re)^2 + std_t(z_im)^2)
z_im_std = (z_im - mean_t(z_im)) / sqrt(std_t(z_re)^2 + std_t(z_im)^2)
```

This keeps the solver real-valued while avoiding independent real/imag scaling
that would distort the relative phase geometry.

Validation:

- A small WSL import test using
  `e10fV1_low50_high250_g2_complex_split_single_obs_info.csv` found 302
  spectrogram real/imag pairs.
- The first pair mapped to zero-based columns 8 and 159, and both columns were
  assigned the same pair scale.

Prepared but not launched because the active abs mixed-GPU run is using about
21 GB of GPU memory:

- `tmp/run_blp_std_complex_pair_gpu_mixedgpu_i2_b2000_e5_e10fV1_complex_split_20260521.sh`

Proposed first complex probe:

```text
dataset_stem      = e10fV1
observable_mode   = complex_split
standardize_mode  = std_complex_pair
solver_name       = resdmd_batch_mixedgpu
batch_size        = 2000
inner_epochs      = 2
epochs            = 5
seed              = 1234
export_mode       = summary_only
```

## Reproducibility Pointers

BLP mixed-GPU one-epoch probes:

- `tmp/run_blp_stdT_gpu_mixedgpu_i2_b5000_e10fV1_abs_20260521.sh`
- `tmp/run_blp_stdT_gpu_mixedgpu_i2_b10000_e10fV1_abs_20260521.sh`

BLP mixed-GPU long run:

- `tmp/run_blp_stdT_gpu_mixedgpu_i2_b5000_e100_e10fV1_abs_20260521.sh`

BLP spectrum plotting:

- `tmp/plot_blp_gpu_mixedgpu_spectrum_20260521.py`
- `tmp/plot_blp_gpu_mixedgpu_spectrum_compare_20260521.py`

BOLD standardized low-lr plots:

- `tmp/plot_bold_p4_stdT_lowlr_loss_curves_20260521.py`

Complex-pair probe:

- `tmp/run_blp_std_complex_pair_gpu_mixedgpu_i2_b2000_e5_e10fV1_complex_split_20260521.sh`

## Open Decisions

- Whether to promote `resdmd_batch_mixedgpu` from probe to a regular solver
  option after the long `batch=5000` curve and spectrum are inspected.
- Whether the BLP standardized coordinate should become the default for abs
  and complex_split P4.
- Whether BOLD standardized P4 needs different parameters rather than just
  standardization, because the current BOLD stdT low-lr runs are not uniformly
  smooth.
- Whether `std_complex_pair` should also be added to BLP/BOLD controller
  scripts, or kept as a manual runner option until the first complex probe is
  inspected.
