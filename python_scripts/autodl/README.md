# AutoDL Remote Entry

This directory is the canonical remote/cloud side of the pipeline 4 training
family. Read it as two current entry branches:

- `pipeline 4.1`: BLP training/export
- `pipeline 4.2`: BOLD training/export

Use these launchers first:

- `run_blp_observables_mlp_reskoopnet.ps1`
- `run_bold_observables_mlp_reskoopnet.ps1`

The shared controller stack under this folder supports the current canonical
observable modes for both branches:

- BLP: `abs`, `complex_split`
- BOLD: `eleHP`, `HP`, `roi_mean`, `slow_band_power`, `svd`, `HP_svd100`, `global_svd100`

Shared controller layers:

- `controller_autodl_reskoopnet_mlp.py`
- `batch_controller_autodl_reskoopnet_mlp.py`
- `dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `multi_dataset_batch_controller_autodl_reskoopnet_mlp.py`
- `run_autodl_reskoopnet_mlp.py`

Files such as `run_E10fV1_4conditions_60epoch.ps1` or
`run_E10gW1_4conditions_60epoch.ps1` are archived dataset-specific launchers,
not the canonical pipeline 4.1 / 4.2 entry.
