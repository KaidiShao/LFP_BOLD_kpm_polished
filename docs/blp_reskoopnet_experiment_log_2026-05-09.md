# BLP ResKoopNet Experiment Log (2026-05-09)

This note records the main BLP ResKoopNet experiments discussed and run during the current iteration, with emphasis on `E10gb1`, solver variants, paper-faithfulness, and current recommendations.

## Scope

- Main task: stabilize `projected_vlambda` training for BLP data, especially `E10gb1 abs`.
- Main comparison axes:
  - seed sensitivity
  - learning rate
  - regularization
  - outer schedule (`dual-sync` vs `torch-like`)
  - `inner_epochs`
  - static-target vs dynamic spectral target
- Paper-validation side tasks:
  - `pendulum`
  - `turbulence`

## Reference Data Provenance

- Official pendulum Drive `data_pendulum_90.mat` matches local:
  - `tmp/figure2d_lab/original_data/data_pendulum_90.mat`
- Official turbulence Drive `pressure_data.mat` matches local:
  - `tmp/turbulence_lab/original_data/pressure_data.mat`
- The repo-regenerated pendulum variant is different:
  - `tmp/figure2d_lab/data_pendulum_90.mat`
- Detailed provenance is documented in:
  - `docs/reskoopnet_reference_data_manifest_2026-05-08.md`

## Solver Variants Considered

- Old local TF solver:
  - `code/algorithms/new_kpm_code_yuanchao_gpu20240702/solver_resdmd2.py`
  - Behavior: effectively trains against a relatively static spectral target.
- Current polished TF solver, dual-sync:
  - `python_scripts/autodl/solver_resdmd_batch3.py`
  - Behavior: pre-train spectral sync plus post-train spectral sync.
- Current polished TF solver, torch-like:
  - same file as above, with `spectral_sync_mode=pre_only`
  - Behavior: pre-train spectral sync only, closer to paper/upstream Torch schedule.
- Current engineered static-target solver:
  - `python_scripts/autodl/solver_resdmd_batch3_static_target.py`
  - Behavior: uses current engineered pipeline but freezes training `V/Lambda` after initial sync.
- Upstream reference solvers:
  - `tmp/ResKoopNet_upstream/solver_resdmd_tf.py`
  - `tmp/ResKoopNet_upstream/solver_resdmd_torch.py`

The higher-level solver-version summary is documented in:
- `docs/reskoopnet_solver_versions_2026-05-08.md`

## Early Seed Sweep on Old BLP Mainline

Configuration:
- residual form: `projected_vlambda`
- observable: `abs`
- old mainline regime before later torch-like tuning

Results:

| seed | best outer val | best epoch |
|---|---:|---:|
| 1234 | 0.0415676 | 8 |
| 2026 | 0.0614613 | 2 |
| 3407 | 0.0479781 | 11 |
| 7777 | 0.0280172 | 7 |
| 9999 | 0.0738381 | 2 |

Takeaway:
- seed sensitivity was clearly large
- the configuration was not yet stable enough for production

## Learning-Rate Coarse Sweep

Setup:
- compared `lr=3e-5` vs `lr=1e-5`
- fixed `reg=0.3`
- fixed `inner_epochs=2`
- short-run diagnostics on `E10gb1`

Key outcome:
- `1e-5` was clearly better than `3e-5`
- `3e-5` tended to show stronger rebound / instability
- `1e-5` became the main candidate learning rate

Representative results:
- `lr=1e-5, seed1234`: best val about `0.02250`
- `lr=1e-5, seed2026`: best val about `0.02725`
- `lr=3e-5, seed1234`: briefly good but then rebounded badly

## Regularization Sweep at `lr=1e-5`

Setup:
- fixed `lr=1e-5`
- fixed `inner_epochs=2`
- compared `reg=0.2, 0.3, 0.5`
- used seeds `1234` and `2026`

Results:

| reg | seed | best outer val | best epoch | final outer val |
|---|---:|---:|---:|---:|
| 0.2 | 1234 | 0.021438 | 5 | 0.068617 |
| 0.2 | 2026 | 0.037714 | 6 | 0.037714 |
| 0.3 | 1234 | 0.022020 | 6 | 0.022020 |
| 0.3 | 2026 | 0.026007 | 6 | 0.026007 |
| 0.5 | 1234 | 0.021288 | 6 | 0.021288 |
| 0.5 | 2026 | 0.038625 | 5 | 0.051971 |

Takeaway:
- `reg=0.3` was the most balanced and stable choice across seeds
- `reg=0.2` could hit good points but drifted more
- `reg=0.5` was too harsh for one of the seeds

## Long-Term `E10gb1` Comparisons

The most important long-term `E10gb1 abs` comparisons are:

| solver variant | key setup | best outer val | best epoch | qualitative behavior |
|---|---|---:|---:|---|
| dual-sync dynamic | `inner=2`, `seed=1234` | 0.0227029 | 6 | good early, later drifted strongly |
| torch-like dynamic | `inner=2`, `seed=1234` | 0.0150615 | 8 | best overall result so far, drift remains but is milder |
| torch-like dynamic | `inner=4`, `seed=1234` | 0.0189667 | 5 | slightly more stable final state, but worse best point |
| torch-like dynamic | `inner=1`, `seed=1234` | 0.0158882 | 12 | very competitive, but still not better than `inner=2` |
| static target | `inner=2`, `seed=1234` | 0.0489045 | 1 | smooth inner loss, poor dynamic outer metric |

Main takeaways:
- `torch-like + inner=2` currently gives the best BLP result on `E10gb1`
- `inner=1` is surprisingly competitive, but still slightly behind `inner=2`
- `inner=4` makes final drift milder, but the best solution is not as deep
- the static-target solver confirms that inner optimization itself is healthy; the main difficulty is dynamic spectral consistency

## Important Plots

Current `inner=1` full global curve:
- `tmp/live_plots/e10gb1_torchlike_ie1_full_global_loss_curve.png`

Final `inner=1` versus `inner=2` comparison:
- `tmp/live_plots/e10gb1_torchlike_inner1_vs_inner2_full_comparison.png`

Earlier solver comparison panels:
- `tmp/live_plots/e10gb1_solver_inner_loss_comparison_with_ie1.png`
- `tmp/live_plots/e10gb1_solver_outer_metric_comparison_with_ie1.png`

These plots were used to support the current choice of:
- `torch-like`
- `inner_epochs=2`
- `lr=1e-5`
- `reg=0.3`

## Static-Target Experiment

We introduced a new engineered static-target solver:
- `python_scripts/autodl/solver_resdmd_batch3_static_target.py`

Purpose:
- keep current batch/chunk/checkpoint/export engineering
- mimic the older TF behavior where the training target is effectively static

Result:
- inner loss becomes smooth and steadily decreasing
- outer dynamic metric does not improve correspondingly

Interpretation:
- the old smooth-loss feel likely came from the static spectral target
- this supports the idea that current instability is mainly driven by dynamic outer spectral updates, not by failure of inner optimization

## Paper-Faithfulness Validation

### Turbulence

Torch-like validation remained strong on the most important paper target:
- mode 1 correlation with official result stayed around `0.9965`
- mode 1 remained the top low-residue mode

Conclusion:
- turbulence mode-1 reproduction is preserved under the torch-like schedule

### Pendulum

When the official/original pendulum data were used:
- torch-like training reproduced the key annulus-like pseudospectrum structure well enough to count as a practical reproduction

Important nuance:
- raw matrix/eigenvalue comparisons can differ noticeably
- pseudospectrum geometry is the more meaningful comparison target here

## Checkpoint and Export Notes

- For completed runs, both `best` and `final` checkpoints exist.
- For interrupted runs, usually only `best` can be relied on.
- The current best `inner=2` run was exported from its best checkpoint into full eigenfunction chunks.

Current full `E10gb1` eigenfunction export from best `inner=2`:
- output folder:
  - `E:\\autodl_results_new\\e10gb1\\mlp\\outputs\\mlp_obs_blp_vlambda_torchlike_longterm_full_20260507_e10gb1_seed1234_projected_vlambda_abs`
- chunk count:
  - `1199` files matching `*_outputs_*.mat`

## Current Recommendation

For the BLP mainline at the current stage, the recommended working configuration is:

- solver: `resdmd_batch`
- schedule: `torch-like` (`spectral_sync_mode=pre_only`)
- seed: `1234`
- learning rate: `1e-5`
- regularization: `0.3`
- `inner_epochs=2`
- use the `best checkpoint` for export, not the final epoch by default

Reasoning:
- this version currently gives the best `E10gb1` result
- it is the closest practical compromise between paper-faithfulness and BLP performance
- it reproduces the critical paper validation experiments well enough

## Open Questions

- whether `batch_size > 2000` can further stabilize torch-like training
- whether `inner=1` can become competitive or superior under a different batch-size regime
- whether a better early-stopping / export policy should be formalized for BLP production runs

## Related Files Written During This Iteration

- `docs/reskoopnet_reference_data_manifest_2026-05-08.md`
- `docs/reskoopnet_solver_versions_2026-05-08.md`
- `python_scripts/autodl/solver_resdmd_batch3_static_target.py`
- `tmp/live_plots/e10gb1_torchlike_ie1_full_global_loss_curve.png`
- `tmp/live_plots/e10gb1_torchlike_inner1_vs_inner2_full_comparison.png`
