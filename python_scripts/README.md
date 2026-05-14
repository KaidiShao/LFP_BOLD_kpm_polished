# Python Scripts Index

This folder is easiest to read as the pipeline 4 training family plus a few
support/archival directories.

## Canonical Pipeline 4 Family

### Pipeline 4.1: BLP training/export

Remote entry:

- `autodl/run_blp_observables_mlp_reskoopnet.ps1`

Local entries:

- `local/run_blp_observables_mlp_reskoopnet_local.ps1`
- `local/run_blp_observables_mlp_reskoopnet_wsl.ps1`
- `local/run_blp_observables_mlp_reskoopnet_wsl.sh`

Current canonical observable modes:

- `abs`
- `complex_split`

### Pipeline 4.2: BOLD training/export

Remote entry:

- `autodl/run_bold_observables_mlp_reskoopnet.ps1`

Local entries:

- `local/run_bold_observables_mlp_reskoopnet_local.ps1`
- `local/run_bold_observables_mlp_reskoopnet_wsl.ps1`
- `local/run_bold_observables_mlp_reskoopnet_wsl.sh`

Current canonical observable modes:

- `eleHP`
- `HP`
- `roi_mean`
- `slow_band_power`
- `svd`
- `HP_svd100`
- `global_svd100`

### Shared controller core

Both branches delegate to the same low-level controller stack:

- `autodl/run_autodl_reskoopnet_mlp.py`
- `autodl/controller_autodl_reskoopnet_mlp.py`
- `autodl/batch_controller_autodl_reskoopnet_mlp.py`
- `autodl/dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `autodl/multi_dataset_batch_controller_autodl_reskoopnet_mlp.py`

## Archived / Other Directories

### `autodl_local_e10gb1`

Archived local experiment reproduction material. Not canonical pipeline 4.1 or 4.2.

### `cnn_reskoopnet`

Optional CNN comparison branch. Not the canonical MLP training line.

### `settings`

Environment/setup files such as conda environment YAMLs.

### `original`

Archived upstream/original reference files kept for comparison.
