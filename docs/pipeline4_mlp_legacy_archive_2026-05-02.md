# Pipeline 4 MLP Legacy Archive

Date: `2026-05-02`

This note records the training-data archive performed after switching the
pipeline 4 ResKoopNet core to the `solver3_updated`-style high-precision
`projected_vlambda` mainline.

## Scope

The archive covered existing `mlp` training artifacts under:

- `E:\autodl_results`
- `E:\autodl_results_local`

For each discovered `mlp` result root, existing children under these buckets
were moved into:

- `mlp\legacy\pre_solver3_updated_2026-05-02\outputs`
- `mlp\legacy\pre_solver3_updated_2026-05-02\logs`
- `mlp\legacy\pre_solver3_updated_2026-05-02\checkpoints`
- `mlp\legacy\pre_solver3_updated_2026-05-02\console_logs`

This keeps the old runs intact while clearing the active `outputs`, `logs`,
`checkpoints`, and `console_logs` buckets for the new mainline.

## Summary

Total moved items: `395`

Archived `mlp` roots and moved bucket counts:

- `E:\autodl_results\bold\e10gb1\mlp`
  - outputs: `3`
  - logs: `3`
- `E:\autodl_results\e10fV1\mlp`
  - outputs: `8`
  - logs: `8`
  - checkpoints: `7`
  - console_logs: `6`
- `E:\autodl_results\e10gb1\mlp`
  - outputs: `7`
  - logs: `6`
  - checkpoints: `6`
- `E:\autodl_results\e10gh1\mlp`
  - outputs: `4`
  - logs: `4`
  - checkpoints: `5`
- `E:\autodl_results\e10gw1\mlp`
  - outputs: `2`
  - logs: `3`
  - checkpoints: `1`
- `E:\autodl_results\f12m01\mlp`
  - outputs: `4`
  - logs: `1`
  - checkpoints: `2`
- `E:\autodl_results_local\bold_wsl\e10fV1\mlp`
  - outputs: `18`
  - logs: `18`
  - checkpoints: `18`
  - console_logs: `18`
- `E:\autodl_results_local\bold_wsl\e10gb1\mlp`
  - outputs: `24`
  - logs: `24`
  - checkpoints: `24`
  - console_logs: `27`
- `E:\autodl_results_local\bold_wsl\e10gh1\mlp`
  - outputs: `18`
  - logs: `18`
  - checkpoints: `18`
  - console_logs: `18`
- `E:\autodl_results_local\bold_wsl\f12m01\mlp`
  - outputs: `18`
  - logs: `18`
  - checkpoints: `18`
  - console_logs: `18`

## New Single-Run Launcher

The new narrow launcher for the current BLP mainline is:

- `python_scripts/local/run_one_blp_reskoopnet_vlambda_wsl.ps1`

It is intentionally limited to:

- residual form: `projected_vlambda`
- observable modes: `abs`, `complex_split`

Change `DatasetStem`, `AbsFilename`, or `ComplexSplitFilename`, then run it to
launch a fresh pair of BLP ResKoopNet jobs on the current pipeline 4 mainline.
