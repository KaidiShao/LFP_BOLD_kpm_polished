# Local AutoDL Copy For E10.gb1

This folder is a local WSL-ready copy of the original `python_scripts/autodl` workflow.
It keeps the original cloud-facing AutoDL files untouched and rewires the main paths for `E10.gb1`.

## Default paths

- Code folder:
  `/mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/python_scripts/autodl_local_e10gb1`
- Observable input root:
  `/mnt/e/DataPons_processed/e10gb1/reskoopnet_dictionary`
- Results root:
  `/mnt/e/autodl_results/e10gb1/mlp`

## Main entry

Open this notebook in WSL/Jupyter:

- `ResKoopNet_pipeline_e10gb1_local_mlp_obs.ipynb`

Example:

```bash
cd /mnt/d/Onedrive/ICPBR/Alberta/koopman_events/LFP_BOLD_kpm_polished/python_scripts/autodl_local_e10gb1
jupyter lab ResKoopNet_pipeline_e10gb1_local_mlp_obs.ipynb
```

## Current training defaults

- `observable_mode = "abs"`
- `EXPERIMENT_NAME = "e10gb1_260415_shuffle_plateau"`
- `extra_epochs = 40`
- `inner_epochs = 3`
- `train_shuffle = True`
- `outer_lr_patience = 8`
- `resume_mode = "fresh"`

The notebook appends `residual_form` automatically when building the effective experiment label.

## Resume behavior

The notebook now supports three explicit modes:

- `resume_mode = "fresh"` starts a new run and clears stale in-memory history
- `resume_mode = "final"` continues from the latest final checkpoint
- `resume_mode = "best"` restores the current best checkpoint for fine-tuning or export

Trainer state is saved separately under:

- `final/training_state.json`
- `best/training_state.json`

These files store loss history, outer-history, best-metric metadata, optimizer learning rate, and scheduler counters.

## Local notes

- `solver_resdmd_batch3.py` is the main default target for `solver_name = "resdmd_batch"`
- summary files are autosaved after each training call
- final and best checkpoints both keep trainer-state sidecars
- inner-loop early stopping restores the best inner weights
- training batches are shuffled by default
- outer learning-rate decay is controlled by a best-based plateau rule

## Change log

See:

- `training_flow_update_2026-04-15.md`
