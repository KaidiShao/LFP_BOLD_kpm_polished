# Local Entry

This directory is the canonical local side of the pipeline 4 training family.
Read it as two current entry branches:

- `pipeline 4.1`: BLP training/export
- `pipeline 4.2`: BOLD training/export

Use these entry points first:

- `run_blp_observables_mlp_reskoopnet_local.ps1`
- `run_blp_observables_mlp_reskoopnet_wsl.ps1`
- `run_blp_observables_mlp_reskoopnet_wsl.sh`
- `run_one_blp_reskoopnet_vlambda_wsl.ps1`
- `run_bold_observables_mlp_reskoopnet_local.ps1`
- `run_bold_observables_mlp_reskoopnet_wsl.ps1`
- `run_bold_observables_mlp_reskoopnet_wsl.sh`

`run_E10fV1_blp_reskoopnet_4conditions_local.ps1` is an archived one-off
launcher, not the canonical pipeline 4.1 entry.

`run_one_blp_reskoopnet_vlambda_wsl.ps1` is the narrow single-dataset launcher
for the current BLP mainline when you only want:

- `projected_vlambda`
- `abs`
- `complex_split`
