# AutoDL Training Flow Update

Date: `2026-04-15`

## Scope

This note records the training-flow changes applied to both the local AutoDL clone and the original cloud-facing AutoDL files.

Updated files:

- `python_scripts/autodl_local_e10gb1/solver_resdmd_batch3.py`
- `python_scripts/autodl_local_e10gb1/ResKoopNet_pipeline_e10gb1_local_mlp_obs.ipynb`
- `python_scripts/autodl/solver_resdmd_batch3.py`
- `python_scripts/autodl/ResKoopNet_pipeline_autodl_mlp_obs.ipynb`

## Main Changes

### 1. Resume modes were made explicit

The notebook now supports:

- `resume_mode = "fresh"`
- `resume_mode = "final"`
- `resume_mode = "best"`

This replaces the earlier implicit resume flow and makes it clear whether a run should:

- start from scratch
- continue from the latest final checkpoint
- fine-tune from the current best checkpoint

### 2. Trainer state is now saved alongside checkpoints

Each run now writes:

- `final/training_state.json`
- `best/training_state.json`

The saved state includes:

- `loss_history`
- `val_loss_history`
- `outer_history`
- best-metric metadata
- current optimizer learning rate
- plateau-scheduler counters
- run metadata such as experiment name, residual form, and round number

This was added so resume behavior is not limited to model weights alone.

### 3. Optimizer state is included in checkpoints

`tf.train.Checkpoint(...)` now stores:

- model
- optimizer
- Koopman matrix `K`
- eigenvectors
- eigenvalues
- regularization term

This makes resumed training much closer to the original optimizer trajectory than the earlier weight-only restore.

### 4. Inner-loop training was stabilized

The inner training callback now uses:

- `EarlyStopping(..., restore_best_weights=True)`

This avoids carrying the final inner-epoch weights forward when a better validation point occurred earlier in the same outer epoch.

### 5. Training-order randomness was added

The training dataset builder now supports:

- `train_shuffle`
- `shuffle_buffer_size`
- `shuffle_seed`
- `reshuffle_each_iteration`

The default notebook settings now enable shuffling for training batches.

### 6. Outer learning-rate decay now uses best-based plateau logic

The old rule was effectively "decay when the current outer metric is worse than the previous one."

The new logic tracks:

- epochs since last best outer metric
- cooldown after a decay
- minimum learning rate

New solver arguments:

- `outer_lr_patience`
- `outer_lr_cooldown`
- `outer_min_lr`

Notebook defaults currently use:

- `inner_epochs = 3`
- `train_shuffle = True`
- `outer_lr_patience = 8`
- `outer_lr_cooldown = 2`
- `outer_min_lr = 1e-5`

### 7. Summary autosave was added after each training call

After each `solver.build(...)` call, the notebook now autosaves a summary MAT file containing:

- loss histories
- outer metric histories
- best outer epoch
- best validation metric
- checkpoint paths
- Koopman modes and eigenvalues when available

### 8. Notebook history mixing is now guarded

The notebook now clears stale in-memory histories when:

- `resume_mode = "fresh"`
- the experiment context changes

This reduces the chance of mixing loss curves from different experiments in the same kernel session.

### 9. Cloud notebook naming was updated

The default cloud experiment name was changed from:

- `f12m01_base`

to:

- `f12m01_260415_shuffle_plateau`

This matches the newer training strategy and makes output directories easier to interpret.

## Review Result

No blocking logic issues were found during this review pass.

One small cleanup was made during review:

- the cloud solver still had a local-only error string mentioning `/mnt/e/DataPons_processed`
- this was changed to use the actual configured read-only root dynamically

## Validation Performed

The following lightweight checks were run:

- notebook JSON parse for the cloud notebook
- Python AST parse for the cloud solver
- Python AST parse for the local solver

These checks confirm that the modified notebook file is structurally valid and that the modified solver sources parse correctly.

## Residual Risks

The updated flow was not re-run end-to-end in this edit pass, so the remaining risk is runtime-only rather than syntax-level.

Main residual risks:

- resume behavior still depends on the actual checkpoint and trainer-state files present on disk
- export cells still contain their own summary/chunk-writing logic, separate from the new autosave helper
- if `log_interval` is changed away from `1`, the saved loss history will remain tied to the logging cadence rather than every possible inner epoch

## Recommended Usage

For a new experiment:

- use `resume_mode = "fresh"`

For continuing the same run:

- use `resume_mode = "final"`

For fine-tuning from the best known checkpoint:

- use `resume_mode = "best"`

Always restart the kernel after updating solver code so the notebook imports the latest solver implementation.
