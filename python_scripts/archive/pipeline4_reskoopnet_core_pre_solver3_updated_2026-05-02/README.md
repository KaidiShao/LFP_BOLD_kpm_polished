# Pipeline 4 ResKoopNet Core Archive

Date: `2026-05-02`

This directory stores pre-migration copies of the pipeline 4 ResKoopNet core
files before the canonical MLP training/export path was updated to the
`solver3_updated`-style default behavior.

Reason for archive:

- preserve the exact pre-change pipeline 4 training/export code
- allow side-by-side comparison with the new `projected_vlambda + float64`
  canonical path
- keep rollback/reference copies outside the active `autodl/` and `local/`
  directories

Archived file groups:

- `python_scripts/autodl/solver_resdmd_batch3.py`
- `python_scripts/autodl/run_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/multi_dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `python_scripts/autodl/run_blp_observables_mlp_reskoopnet.ps1`
- `python_scripts/autodl/run_bold_observables_mlp_reskoopnet.ps1`
- `python_scripts/local/run_blp_observables_mlp_reskoopnet_local.ps1`
- `python_scripts/local/run_bold_observables_mlp_reskoopnet_local.ps1`
- `python_scripts/local/run_blp_observables_mlp_reskoopnet_wsl.ps1`
- `python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.ps1`
- `python_scripts/local/run_blp_observables_mlp_reskoopnet_wsl.sh`
- `python_scripts/local/run_bold_observables_mlp_reskoopnet_wsl.sh`

Migration intent of the active code after this archive:

- make `projected_vlambda` the default canonical residual form
- make high-precision `float64 / complex128` the default canonical numeric path
- keep `projected_kv` available as an explicit comparison branch
