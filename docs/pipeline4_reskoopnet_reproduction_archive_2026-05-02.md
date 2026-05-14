# Pipeline 4 ResKoopNet Reproduction Archive

Date: `2026-05-02`

This note archives the ResKoopNet reproduction work completed in this thread.
The focus was to compare the current project `solver3`, an isolated
`solver3_updated` high-precision variant, and the upstream ResKoopNet solver on
the public pendulum and turbulence experiments.

## High-Level Conclusion

- The project should treat `solver3_updated`-style `projected_vlambda` as the
  canonical pipeline 4 line.
- `projected_kv` remains useful as a comparison branch, but not as the default
  publication-style reproduction branch.
- Pendulum reproduces well with the updated high-precision solver.
- Turbulence is more sensitive. The current project `solver3` and the updated
  `projected_kv` branch recover the leading mode shape, but their low-residue
  mode ordering is not as reliable as the upstream/original-style baseline.

## Reference Data And Labs

Main isolated workspaces used during reproduction:

- `tmp/figure2d_lab`
- `tmp/turbulence_lab`
- `tmp/ResKoopNet_upstream`

Main reference data files:

- `tmp/figure2d_lab/original_data/data_pendulum_90.mat`
- `tmp/turbulence_lab/original_data/pressure_data.mat`
- `tmp/turbulence_lab/original_data/turbulence_resdmd_250basis_official.mat`

## Test Matrix

| ID | Experiment | Solver variant | Residual form | Status | Main artifact root |
|---|---|---|---|---|---|
| P1 | Pendulum | project `solver3` | `projected_vlambda` | reproduced with imperfect annulus | `tmp/figure2d_lab/outputs/solver3_current_originaldata_vlambda_seed1234` |
| P2 | Pendulum | `solver3_updated` | `projected_vlambda` | reproduced, clean annulus | `tmp/figure2d_lab/outputs/solver3_updated_originaldata_seed1234` |
| P3 | Pendulum | `solver3_updated` | `projected_kv` | reproduced, slightly thicker annulus | `tmp/figure2d_lab/outputs/solver3_updated_originaldata_kv_seed1234` |
| T1 | Turbulence | upstream solver | upstream/original branch | reproduced and matched official baseline closely | `tmp/turbulence_lab/outputs/turbulence_reskoopnet_full` |
| T2 | Turbulence | project `solver3` | `projected_kv` | partial reproduction only | `tmp/turbulence_lab/outputs/turbulence_solver3_kv_full` |
| T3 | Turbulence | project `solver3` | `projected_vlambda` | partial reproduction only, nearly identical to T2 | `tmp/turbulence_lab/outputs/turbulence_solver3_vlambda_full` |
| T4 | Turbulence | `solver3_updated` | `projected_kv` | partial reproduction, leading mode shape recovered | `tmp/turbulence_lab/outputs/turbulence_solver3_updated_kv_full` |

## Pendulum Summary

### P1. Project solver3 + original data + projected_vlambda

Artifact root:

- `tmp/figure2d_lab/outputs/solver3_current_originaldata_vlambda_seed1234`

Key files:

- `solver_resdmd_batch3_current_seed1234.mat`
- `solver_resdmd_batch3_current_seed1234_summary.json`
- `pseudospectrum_full/pseudospectrum_threshold_style.png`

Key observations:

- `best_val_metric = 0.3540107421875`
- `condG = 2.468157886199638e+20`
- `max_abs_eigenvalue = 0.9986763349849235`
- visually close to Figure 2(d), but annulus has gaps and rough edges

### P2. solver3_updated + original data + projected_vlambda

Artifact root:

- `tmp/figure2d_lab/outputs/solver3_updated_originaldata_seed1234`

Key files:

- `solver_resdmd_batch3_updated_seed1234.mat`
- `solver_resdmd_batch3_updated_seed1234_summary.json`
- `pseudospectrum_full/pseudospectrum_threshold_style.png`

Key observations:

- `best_val_metric = 0.4033256607629385`
- `condG = 7.251642355939045e+19`
- `max_abs_eigenvalue = 0.9996312177552443`
- produced the cleanest annulus-style continuous spectrum in this thread

### P3. solver3_updated + original data + projected_kv

Artifact root:

- `tmp/figure2d_lab/outputs/solver3_updated_originaldata_kv_seed1234`

Key files:

- `solver_resdmd_batch3_updated_seed1234.mat`
- `solver_resdmd_batch3_updated_seed1234_summary.json`
- `pseudospectrum_full/pseudospectrum_threshold_style.png`

Key observations:

- `best_val_metric = 0.42841701933263016`
- `condG = 1.8843663756888552e+20`
- `max_abs_eigenvalue = 0.9996313171040336`
- still produced a complete annulus
- annulus was slightly thicker and less clean than `projected_vlambda`

### Pendulum comparison notes

- `solver3_updated(projected_vlambda)` is the strongest pendulum reproduction
  line.
- `solver3_updated(projected_kv)` is acceptable as a comparison branch.
- the current project `solver3` can approach the target shape on the original
  data, but it is less clean and less numerically faithful.

## Turbulence Summary

### T1. Upstream solver baseline

Artifact root:

- `tmp/turbulence_lab/outputs/turbulence_reskoopnet_full`

Key files:

- `turbulence_resdmd_250basis.mat`
- `summary.json`
- `comparison_vs_official.json`
- `plots/modes/koopman_mode_1.png`

Key observations:

- `N_dict = 250`
- `runtime_seconds = 406.35574436187744`
- `max_abs_eigenvalue = 1.0024483222606626`
- `psi_x_rel_l2_diff vs official = 8.563128054892824e-05`
- `mode 1 correlation vs official mode 1 = 0.9965135534966361`
- this was the strongest turbulence reproduction and the reference baseline for
  later comparisons

### T2. Project solver3 + projected_kv

Artifact root:

- `tmp/turbulence_lab/outputs/turbulence_solver3_kv_full`

Key files:

- `turbulence_solver3_projected_kv_250basis.mat`
- `summary.json`
- `plots/modes/summary.json`

Key observations:

- `outer_epochs_completed = 3`
- `best_val_metric = 48932992.535564855`
- top residue modes: `[3, 4, 5, 1]`
- recovered the leading physical mode shape, but low-residue mode ordering was
  not close to the upstream baseline

### T3. Project solver3 + projected_vlambda

Artifact root:

- `tmp/turbulence_lab/outputs/turbulence_solver3_vlambda_full`

Key files:

- `turbulence_solver3_projected_vlambda_250basis.mat`
- `summary.json`
- `plots/modes/summary.json`

Key observations:

- `outer_epochs_completed = 3`
- `best_val_metric = 31943465.774058577`
- numerically almost identical to T2 in spectrum and mode content
- changing `kv -> vlambda` inside the current project solver alone did not fix
  turbulence reliability

### T4. solver3_updated + projected_kv

Artifact root:

- `tmp/turbulence_lab/outputs/turbulence_solver3_updated_kv_full`

Key files:

- `turbulence_solver3_projected_kv_250basis.mat`
- `summary.json`
- `plots/modes/summary.json`
- `plots/modes/koopman_mode_3.png`

Key observations:

- `outer_epochs_completed = 3`
- `best_val_metric = 31943457.0026029`
- top residue modes: `[3, 4, 5, 1]`
- top residues: `[1.9951765278121564e-06, 0.021851391459334932, 0.021851391459334932, 0.021887624362469525]`
- recovered the leading spatial mode shape well
- `mode 3 correlation vs official mode 1 = -0.9922727587740682`
- still did not reproduce the upstream low-residue mode ranking

### Turbulence comparison notes

- `projected_kv` can recover the dominant turbulence mode shape, but it is not
  the most reliable line for low-residue mode ranking.
- The upstream/original-style branch remains the strongest turbulence
  reproduction reference in this thread.
- For pipeline defaults, `projected_vlambda` should be preferred as the main
  line and `projected_kv` should remain optional.

## Files Used To Support The Pipeline 4 Migration

Archived pre-change pipeline 4 core files:

- `python_scripts/archive/pipeline4_reskoopnet_core_pre_solver3_updated_2026-05-02/`

Current active canonical solver target after this thread:

- `python_scripts/autodl/solver_resdmd_batch3.py`

Canonical pipeline 4 defaults intended after this thread:

- default residual form: `projected_vlambda`
- default numeric policy: `float64`
- keep `projected_kv` available as an explicit comparison branch
