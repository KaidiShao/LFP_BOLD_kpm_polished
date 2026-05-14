## ResKoopNet Solver Versions

Date: `2026-05-08`

This note records the main solver variants discussed and tested in the current
ResKoopNet reproduction and BLP adaptation work. The goal is to keep the
different codepaths, training schedules, and observed behaviors separate.

### Summary Table

| Version | File | Core schedule | Main behavior | Current role |
| --- | --- | --- | --- | --- |
| Old local TF solver | `code/algorithms/new_kpm_code_yuanchao_gpu20240702/solver_resdmd2.py` | older TF graph; effectively more static spectral target in `projected_vlambda` training | felt stable on `E10gb1`, but is less faithful to the alternating optimization in the paper | historical local reference |
| Polished TF solver, dual-sync | `python_scripts/autodl/solver_resdmd_batch3.py` | `compute_K -> eig -> sync -> train_psi -> compute_K -> eig -> sync` every outer epoch | preserves current `pendulum` and `turbulence` reproduction behavior, but long BLP runs drift after finding a good early point | current reproduction mainline |
| Polished TF solver, torch-like | `python_scripts/autodl/solver_resdmd_batch3.py` with `spectral_sync_mode=pre_only` | `compute_K -> eig -> sync -> train_psi`, then skip post-train spectral refresh | improves early-to-mid BLP behavior on `E10gb1`, but still drifts later; does not preserve the current `pendulum` solution family | BLP experiment branch |
| Upstream TF solver | `tmp/ResKoopNet_upstream/solver_resdmd_tf.py` | official TensorFlow implementation from the repository | useful as a reference implementation, but not the current engineering mainline for BLP | reference only |
| Upstream Torch solver | `tmp/ResKoopNet_upstream/solver_resdmd_torch.py` | official PyTorch implementation; small-batch, `float64/complex128`, outer alternating updates | closest in spirit to the paper, but in the current workspace it does not automatically reproduce the local `pendulum` baseline | reference only |

### 1. Old Local TF Solver

File:
`code/algorithms/new_kpm_code_yuanchao_gpu20240702/solver_resdmd2.py`

Observed properties:

- This is the older local line that previously felt stable on `E10gb1`.
- In `projected_vlambda`, the effective training target is more static than in
  the newer alternating versions.
- Because the spectral target is not refreshed in the same aggressive way as
  the polished solver, optimization can look smoother on BLP data.

Current interpretation:

- Good as a historical behavior reference.
- Not the best choice for paper-faithful validation.

### 2. Polished TF Solver, Dual-Sync Mainline

File:
`python_scripts/autodl/solver_resdmd_batch3.py`

Default outer schedule:

1. `compute_K`
2. `eig_decomp`
3. sync spectral quantities into the training graph
4. `train_psi`
5. `compute_K`
6. `eig_decomp`
7. sync again

Observed properties:

- This version currently preserves the local `pendulum` pseudospectrum baseline.
- This version also preserves the key `turbulence` reproduction behavior,
  especially the first mode.
- On BLP, especially `E10gb1` and `E10gh1`, it can find a good early point but
  often drifts on longer runs.

Current interpretation:

- Best current choice for reproduction work that must preserve `pendulum`.
- Not ideal for long BLP runs without additional stabilization.

### 3. Polished TF Solver, Torch-Like Branch

File:
`python_scripts/autodl/solver_resdmd_batch3.py`

Relevant switches:

- `spectral_sync_mode=pre_only`
- optional `train_shuffle=True`

Torch-like outer schedule:

1. `compute_K`
2. `eig_decomp`
3. sync spectral quantities into the training graph
4. `train_psi`
5. skip the post-train spectral refresh inside that outer epoch

Observed properties:

- On BLP `E10gb1`, this branch produced a better best point than the dual-sync
  version in the tested long-run comparisons.
- It reduced the severity of late drift, but did not eliminate long-term drift.
- It did not preserve the current `pendulum` solution family under the present
  validation setup.
- It still preserved the most important `turbulence` mode-1 behavior.

Current interpretation:

- Good candidate branch for BLP-specific experimentation.
- Not suitable as the only global solver mode if `pendulum` must remain intact.

### 4. Upstream TensorFlow Solver

File:
`tmp/ResKoopNet_upstream/solver_resdmd_tf.py`

Observed properties:

- Official repository TensorFlow implementation.
- Useful for understanding the repository's original TensorFlow design.
- Not currently used as the workspace mainline for BLP production runs.

Current interpretation:

- Keep as a reference implementation, not as the default working solver.

### 5. Upstream Torch Solver

File:
`tmp/ResKoopNet_upstream/solver_resdmd_torch.py`

Observed properties:

- Official repository PyTorch implementation.
- Uses `torch.float64` for real quantities and `torch.complex128` for complex
  spectral quantities.
- Uses a smaller-batch training recipe than the current BLP-oriented TF runs.
- More closely reflects the paper-style alternating flavor than the local old
  solver.

Observed behavior in the current workspace:

- The important `turbulence` first mode remains reproducible.
- The current local `pendulum` baseline is not automatically recovered just by
  switching to the official Torch implementation.
- This means framework alone does not explain the differences; the full recipe
  and data choice matter.

Current interpretation:

- Important paper reference.
- Not recommended as a drop-in replacement for the current production workflow.

### Dataset Notes

These solver comparisons depend on using the correct official reference data.

- Official Drive `pendulum_data/data_pendulum_90.mat` matches:
  `tmp/figure2d_lab/original_data/data_pendulum_90.mat`
- Official Drive `turbulence_data/pressure_data.mat` matches:
  `tmp/turbulence_lab/original_data/pressure_data.mat`

See also:

- `docs/reskoopnet_reference_data_manifest_2026-05-08.md`
- `tmp/reference_data/reskoopnet/manifest.json`

### Current Recommendation

If the priority is:

1. Preserve `pendulum` and `turbulence` reproduction:
   use the polished TF solver in dual-sync mode.
2. Explore better BLP training dynamics:
   use the polished TF solver in torch-like mode as an experiment branch.
3. Understand paper-faithful behavior:
   compare against the upstream Torch implementation, but do not assume it can
   replace the current mainline without additional recipe alignment.
