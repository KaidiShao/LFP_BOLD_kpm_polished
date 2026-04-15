# Local AutoDL Training Flow Update

Date: `2026-04-15`

## Scope

This note records the local-only updates applied in `python_scripts/autodl_local_e10gb1`.

Files updated:

- `ResKoopNet_pipeline_e10gb1_local_mlp_obs.ipynb`
- `solver_resdmd_batch3.py`
- `README.md`

## Why This Local Copy Was Changed

The original goal of the local copy was to run the AutoDL workflow against the `E10.gb1` observables under WSL without modifying the original cloud-facing `python_scripts/autodl` directory.

Later, the local workflow was refined to make continuing, exporting, and comparing experiments safer and easier.

## Local Path Rewiring

The local notebook and solver are configured around:

- observable input root:
  `/mnt/e/DataPons_processed/e10gb1/reskoopnet_dictionary`
- writable results root:
  `/mnt/e/autodl_results/e10gb1/mlp`

This keeps local runs separate from the original AutoDL defaults.

## Main Functional Changes

### 1. Experiment naming was cleaned up

The notebook uses a clearer experiment label and appends `residual_form` automatically to the effective run name.

Current default:

- `EXPERIMENT_NAME = "e10gb1_260415_shuffle_plateau"`

This avoids duplicated `abs` tags inside experiment names while still keeping the observable type in the final run label.

### 2. Resume behavior was made explicit

The notebook now uses:

- `resume_mode = "fresh"`
- `resume_mode = "final"`
- `resume_mode = "best"`

This replaced the older "resume if checkpoints exist" style flow.

### 3. Trainer state is saved separately

Each run now stores:

- `final/training_state.json`
- `best/training_state.json`

These sidecars save:

- `loss_history`
- `val_loss_history`
- `outer_history`
- best-outer metrics
- optimizer learning rate
- plateau counters
- run metadata

This makes resumed runs much easier to reason about after a kernel restart.

### 4. Optimizer state now survives checkpoint restore

The local solver now places `optimizer` inside `tf.train.Checkpoint(...)`, instead of restoring only the model and Koopman-state tensors.

This was added so resumed training is closer to the original optimization path.

### 5. Inner-loop training was stabilized

The solver now uses:

- `EarlyStopping(..., restore_best_weights=True)`

This prevents the last inner-epoch weights from replacing a better inner validation point within the same outer epoch.

### 6. Training-order randomness was added

The local solver now supports shuffled training batches with:

- `train_shuffle`
- `shuffle_buffer_size`
- `shuffle_seed`
- `reshuffle_each_iteration`

The current notebook default is:

- `train_shuffle = True`

### 7. Outer learning-rate decay was changed

The old flow reacted too quickly to short-term metric changes.

The local solver now uses a best-based plateau rule with:

- `outer_lr_patience`
- `outer_lr_cooldown`
- `outer_min_lr`

Current notebook defaults:

- `outer_lr_decay_factor = 0.5`
- `outer_lr_patience = 8`
- `outer_lr_cooldown = 2`
- `outer_min_lr = 1e-5`

### 8. Default training settings were adjusted

The current local notebook defaults are:

- `extra_epochs = 40`
- `inner_epochs = 3`
- `train_batch_size = 2000`
- `initial_lr = 1e-4`

This reflects the decision to shorten inner epochs and use shuffle plus plateau-based decay for the next round of experiments.

### 9. Summary autosave was added

After each training call, the notebook now autosaves a summary MAT file containing:

- loss histories
- outer metric histories
- best outer epoch
- best validation metric
- checkpoint references
- Koopman outputs that are already available

### 10. Best-checkpoint export was preserved as a workflow

The local notebook keeps the pattern:

- restore `best` checkpoint
- export summary and `efuns`
- skip `Psi_X / Psi_Y` when those files are too large

This was useful for the `E10.gb1` run because the best outer checkpoint was preferred over the final plateau-state checkpoint.

## Cleanup Included In This Pass

The local `README.md` was updated because it still mentioned the old flag:

- `allow_resume_existing_run`

That flag is no longer the main notebook interface. The README now points to `resume_mode` instead.

## Review Notes

No blocking logic issues were identified during the local review pass.

One small cleanup was also applied:

- the writable-path guard in the local solver now reports the configured read-only root dynamically instead of using a stale hard-coded path in the error string

## Recommended Usage

For a new local experiment:

- use `resume_mode = "fresh"`

For continuing the same local run:

- use `resume_mode = "final"`

For fine-tuning or exporting from the best checkpoint:

- use `resume_mode = "best"`

Always restart the Jupyter kernel after editing solver code so the notebook imports the latest solver implementation.
