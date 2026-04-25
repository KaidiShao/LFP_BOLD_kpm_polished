# AutoDL Automatic Pipeline Setup

Date: `2026-04-18`

## Scope

This note records the current working automatic AutoDL pipeline for the
ResKoopNet MLP workflow under:

- `python_scripts/autodl/`

It covers:

- local controller structure
- remote directory assumptions
- current default training/export settings
- tested automation layers
- current production launcher for `E10fV1` and `E10gH1`
- current remote/local workload split

This note is about the Python AutoDL training/export pipeline, not the upstream
MATLAB preprocessing steps.

## Current Script Layers

Main training/export entry:

- `python_scripts/autodl/run_autodl_reskoopnet_mlp.py`

Per-run controller:

- `python_scripts/autodl/controller_autodl_reskoopnet_mlp.py`

Per-observable batch controller:

- `python_scripts/autodl/batch_controller_autodl_reskoopnet_mlp.py`

Per-dataset batch controller:

- `python_scripts/autodl/dataset_batch_controller_autodl_reskoopnet_mlp.py`

Multi-dataset batch controller:

- `python_scripts/autodl/multi_dataset_batch_controller_autodl_reskoopnet_mlp.py`

Current real-run launcher for two new datasets:

- `python_scripts/autodl/run_E10fV1_E10gH1_60epoch.ps1`

Current local complex-split batch entry:

- `python_scripts/autodl_local_e10gb1/run_local_complex_split_batch.py`
- `python_scripts/autodl_local_e10gb1/run_E10gb1_E10fV1_E10gH1_complex_split_60epoch.sh`

Main solver used by the current pipeline:

- `python_scripts/autodl/solver_resdmd_batch3.py`

Support utilities:

- `python_scripts/autodl/edmd_utils.py`

## Remote Assumptions

The current pipeline assumes:

- AutoDL host is reachable by SSH
- remote code is already uploaded manually
- current code root on AutoDL is:
  - `/root/autodl-tmp/code/reskoopnet_mlp`
- remote Python is:
  - `/root/miniconda3/envs/reskoopnet/bin/python`
- writable remote data/output roots are:
  - `/root/autodl-tmp/data`
  - `/root/autodl-tmp/outputs`
  - `/root/autodl-tmp/checkpoints`
  - `/root/autodl-tmp/logs`

The pipeline intentionally avoids `/root/autodl-pub`.

## SSH Assumptions

Current Windows-side setup assumes:

- host: `connect.westb.seetacloud.com`
- port: `19241`
- user: `root`
- SSH key:
  - `C:\Users\kdsha\.ssh\id_ed25519_autodl`

The controllers can also use the SSH alias if desired, but the current
production commands use explicit host/port/key values.

## Input/Output Hierarchy

The automation is structured in three layers:

1. `dataset`
2. `observable package`
3. `run job`

Meaning:

- one `dataset` contains two observable packages:
  - `abs`
  - `complex_split`
- one observable package expands into two run jobs:
  - `projected_kv`
  - `projected_vlambda`

So one dataset normally produces four training conditions.

## Current Default Training Settings

The currently working defaults in the automation stack are:

- model family: `MLP`
- solver: `resdmd_batch`
- file type: `.h5`
- field name: `obs`
- observable modes:
  - `abs`
  - `complex_split`
- residual forms:
  - `projected_kv`
  - `projected_vlambda`
- layer sizes:
  - `100 100 100`
- `n_psi_train = 100`
- `train_ratio = 0.7`
- `reg = 0.1`
- `rounds = 1`
- `batch_size = 2000`
- `lr = 1e-4`
- `log_interval = 1`
- `lr_decay_factor = 0.8`
- `inner_epochs = 5`
- `end_condition = 1e-9`
- `chunk_size = 5000`
- selected device: `gpu`

Current production epoch count for the new real runs:

- `epochs = 60`

## Export Settings

Current export behavior is:

- export is performed from the `best checkpoint`
- `Psi_X/Psi_Y` export is disabled by default
- normal output chunks and summary MAT files are still exported

The following metadata is written into exported outputs:

- `observable_tag`
- `observable_mode`
- `data_filename`
- `dataset_stem`
- `run_label`

## Cleanup Policy

Current cleanup policy is intentionally split into two levels.

Run-level cleanup:

- enabled after successful local download and minimal verification
- deletes remote:
  - `outputs/<run_label>`
  - `checkpoints/<run_label>`
  - `logs/<run_label>`

Observable-level cleanup:

- optional
- when enabled, runs after all residual-form jobs for that observable succeed
- deletes remote input observable files under:
  - `/root/autodl-tmp/data/<remote_data_subdir>/...`
- also attempts to remove the now-empty observable input directory

Current real-run launcher for `E10fV1` and `E10gH1` enables both:

- `--delete-remote-run-after-download`
- `--delete-remote-input-after-observable`

## Current Remote/Local Split

Current recommended split is now:

- remote AutoDL runs:
  - `abs`
- local GPU runs:
  - `complex_split`

Reason:

- formal `complex_split` observables are much larger than `abs`
- uploading them to AutoDL was becoming the main bottleneck
- `abs` remains the more practical branch for cloud execution under current upload bandwidth

## Local Results Layout

By default, local downloads are organized under:

- `E:\autodl_results\<dataset_stem>\mlp`

Typical contents include:

- controller manifest JSON
- batch manifest JSON
- dataset-batch manifest JSON
- downloaded `outputs/`
- downloaded `logs/`

## Tested Automation Levels

The following levels have been tested successfully with smoke inputs:

- single run controller
- single observable batch
- single dataset batch
- multi-dataset batch

Verified smoke coverage reached:

- `2 dataset x 2 observable x 2 residual_form`

Smoke experiment names used during validation:

- `smoke_multi01`
- `smoke_multi02`

Important note:

- the `f12m01` smoke inputs used in the multi-dataset orchestration test were
  placeholder smoke files for scheduling validation only
- this does not affect the logic of the automation stack itself

## Real Dataset Launcher: E10fV1 And E10gH1

Current production launcher:

- `python_scripts/autodl/run_E10fV1_E10gH1_60epoch.ps1`

This script sequentially runs:

- dataset `e10fV1`
- dataset `e10gh1`

For each dataset it runs:

- `abs + projected_kv`
- `abs + projected_vlambda`

It is currently configured with:

- `epochs = 60`
- fresh checkpoints by default
- resume mode available through `-Resume`
- remote run cleanup enabled
- remote observable-input cleanup enabled

This remote launcher no longer runs `complex_split`.

## Local Complex-Split Launcher

Current local complex-split batch covers:

- `e10gb1`
- `e10fV1`
- `e10gh1`

For each dataset it runs:

- `complex_split + projected_kv`
- `complex_split + projected_vlambda`

Current WSL shell launcher:

- `python_scripts/autodl_local_e10gb1/run_E10gb1_E10fV1_E10gH1_complex_split_60epoch.sh`

## Current Real Input Paths

`E10fV1`:

- `E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_abs_single.mat`
- `E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_abs_single_obs_info.csv`
- `E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_complex_split_single.mat`
- `E:\DataPons_processed\E10fV1\reskoopnet_dictionary\e10fV1_low50_high250_g2_complex_split_single_obs_info.csv`

`E10gH1`:

- `E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_abs_single.mat`
- `E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_abs_single_obs_info.csv`
- `E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_complex_split_single.mat`
- `E:\DataPons_processed\E10gH1\reskoopnet_dictionary\e10gh1_low50_high250_g2_complex_split_single_obs_info.csv`

Note:

- the `E10gH1` filesystem names are currently lowercase as `e10gh1_*`
- the launcher already matches the actual on-disk filenames

## Recommended Usage

Start a fresh real run:

```powershell
powershell -ExecutionPolicy Bypass -File "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished\python_scripts\autodl\run_E10fV1_E10gH1_60epoch.ps1" -ExperimentName real_batch60_run01
```

Resume an existing experiment name:

```powershell
powershell -ExecutionPolicy Bypass -File "D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished\python_scripts\autodl\run_E10fV1_E10gH1_60epoch.ps1" -ExperimentName real_batch60_run01 -Resume
```

## Known Caveats

- remote code sync is not automatic; code must already be uploaded manually to
  `/root/autodl-tmp/code/reskoopnet_mlp`
- `Psi` export is intentionally disabled to save space
- the current automation is optimized for reliable upload/run/download/delete
  flow, not for minimizing total wall time across datasets
- if a run is interrupted before local verification, remote run products may
  remain and should be inspected before manual deletion

## Bottom Line

The current automatic pipeline is now stable enough for:

- local preprocessing
- automatic upload of observables
- remote training/export
- local download and minimal verification
- automatic deletion of remote run artifacts
- automatic deletion of remote input observables after each observable branch finishes

The current recommended production entrypoint for the two new datasets is:

- `python_scripts/autodl/run_E10fV1_E10gH1_60epoch.ps1`
