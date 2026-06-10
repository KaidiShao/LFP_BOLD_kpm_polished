# BOLD P4 Solver4 Time-Standardize Parameter Branch

Date: 2026-05-20

This note records a small side-branch experiment for BOLD P4 ResKoopNet training.
The goal was to test whether the current solver4 can produce a smooth decreasing
loss curve when the input observable is standardized along time, and to identify
which training parameters matter most.

## Motivation

The Kelly side project showed a much smoother ResKoopNet loss curve when each
observable dimension was standardized over time. The current BOLD P4 mainline
had mostly been trained without this input standardization, with much smaller
learning rates and `inner_epochs=1`. This made it unclear whether the unstable
or slow loss curves were caused by the solver itself, the lack of normalization,
or the old training hyperparameters.

The specific question here was:

- Can the current P4 `resdmd_batch4` solver also produce a smooth loss curve with
  time-wise standardization?
- If yes, which parameters control the improvement most strongly?

## Scope

Dataset:

- `E10.gb1`

Observables:

- `roi_mean`
- `gsvd100_ds`

This experiment is a side branch only. It does not overwrite or relabel any
existing mainline or legacy runs.

## Code And Outputs

Run script:

- `tmp/run_bold_solver4_stdT_param_branch_20260520.py`

Plot script:

- `tmp/plot_bold_solver4_stdT_param_branch_20260520.py`

Run log:

- `tmp/run_logs/bold_solver4_stdT_param_branch_20260520.log`

Summary JSON:

- `tmp/run_logs/bold_solver4_stdT_param_branch_20260520_summary.json`

Loss curve figure:

- `tmp/figures/bold_solver4_stdT_param_branch_e10gb1_20260520.png`
- `tmp/figures/bold_solver4_stdT_param_branch_e10gb1_20260520.pdf`

Run labels use the prefix:

```text
mlp_obs_bold_batch4_stdTparam_<observable>_<variant>_20260520_e10gb1_projected_vlambda_<observable>
```

## Common Settings

All runs used:

```text
solver_name        = resdmd_batch4
selected_device    = gpu
training_policy    = float64
analysis_dtype     = float64
gram_dtype         = float64
spectral_dtype     = float64
residual_form      = projected_vlambda
standardize_data   = true
standardize_eps    = 1e-8
train_ratio        = 0.7
seed               = 1234
spectral_sync_mode = pre_only
export_mode        = summary_only
fresh_checkpoints  = true
```

Time-wise standardization means:

```text
X_std[:, j] = (X[:, j] - mean_t(X[:, j])) / std_t(X[:, j])
```

The mean and scale are saved in the run metadata/export payload by the current
runner. Because the coordinate system changes, standardized losses should not be
compared directly to raw-coordinate mainline losses.

## Variants

The baseline was chosen to match the Kelly-style training regime as closely as
possible while still using the current solver4:

```text
lr           = 5e-4
reg          = 0.05
batch_size   = 1024
inner_epochs = 5
n_psi_train  = 64
layer_sizes  = [128, 128]
```

The other runs changed one factor at a time:

| Variant | Changed Factor | Epochs |
|---|---:|---:|
| `kellylike` | baseline | 100 |
| `inner1` | `inner_epochs: 5 -> 1` | 40 |
| `lr1e4` | `lr: 5e-4 -> 1e-4` | 40 |
| `lr1e3` | `lr: 5e-4 -> 1e-3` | 40 |
| `reg1e2` | `reg: 0.05 -> 0.01` | 40 |
| `reg1e1` | `reg: 0.05 -> 0.1` | 40 |
| `batch2000` | `batch_size: 1024 -> 2000` | 40 |
| `oldcap` | `n_psi_train: 64 -> 100`, `layers: [128,128] -> [100,100,100]` | 40 |

## Results

### `roi_mean`

| Variant | Best Val | Best Epoch | Final Val | Final / First |
|---|---:|---:|---:|---:|
| `kellylike` | 0.037108 | 52 | 0.048465 | 0.1493 |
| `inner1` | 0.141332 | 40 | 0.141332 | 0.2927 |
| `lr1e4` | 0.074892 | 40 | 0.074892 | 0.2228 |
| `lr1e3` | 0.038886 | 37 | 0.042567 | 0.1216 |
| `reg1e2` | 0.026190 | 39 | 0.028326 | 0.1134 |
| `reg1e1` | 0.052090 | 35 | 0.086031 | 0.2300 |
| `batch2000` | 0.057294 | 39 | 0.064097 | 0.1716 |
| `oldcap` | 0.037076 | 40 | 0.037076 | 0.2012 |

Best short-run result:

```text
roi_mean reg1e2: reg=0.01, lr=5e-4, inner=5, batch=1024
best val = 0.026190 @ epoch 39
```

### `gsvd100_ds`

| Variant | Best Val | Best Epoch | Final Val | Final / First |
|---|---:|---:|---:|---:|
| `kellylike` | 1.634820 | 100 | 1.634820 | 0.0990 |
| `inner1` | 11.609271 | 40 | 11.609271 | 0.6209 |
| `lr1e4` | 10.906720 | 32 | 10.920173 | 0.5801 |
| `lr1e3` | 1.208226 | 40 | 1.208226 | 0.0786 |
| `reg1e2` | 2.783838 | 40 | 2.783838 | 0.1688 |
| `reg1e1` | 2.966356 | 40 | 2.966356 | 0.1795 |
| `batch2000` | 4.735344 | 40 | 4.735344 | 0.2758 |
| `oldcap` | 0.828053 | 40 | 0.828053 | 0.2726 |

Best short-run result:

```text
gsvd100_ds oldcap: n_psi_train=100, layers=[100,100,100], lr=5e-4, reg=0.05
best val = 0.828053 @ epoch 40
```

The best single-factor learning-rate run was:

```text
gsvd100_ds lr1e3: lr=1e-3
best val = 1.208226 @ epoch 40
```

## Interpretation

The main conclusion is that solver4 is not the obstacle. With time-wise
standardization and Kelly-like training settings, solver4 produces smooth
decreasing validation curves.

Parameter effects:

- `standardize_data=true` changes the training regime substantially. Old raw-data
  hyperparameters should not be reused directly.
- `inner_epochs=5` is important. Reducing to `inner_epochs=1` made both
  observables much worse after 40 epochs.
- `batch_size=1024` was better than `batch_size=2000` in this test.
- `roi_mean` preferred weaker ridge regularization, with `reg=0.01` giving the
  best result.
- `gsvd100_ds` benefited from either a higher learning rate (`lr=1e-3`) or the
  older larger-capacity model (`n_psi_train=100`, `[100,100,100]`).
- Several runs still had best validation loss at the final epoch, so they are
  not fully converged.

## Follow-Up Runs To Prioritize

The most useful next longer runs are:

```text
roi_mean:
  standardize_data = true
  solver_name      = resdmd_batch4
  lr               = 5e-4
  reg              = 0.01
  batch_size       = 1024
  inner_epochs     = 5
  n_psi_train      = 64
  layer_sizes      = [128, 128]

gsvd100_ds:
  standardize_data = true
  solver_name      = resdmd_batch4
  lr               = 5e-4 or 1e-3
  reg              = 0.05
  batch_size       = 1024
  inner_epochs     = 5
  n_psi_train      = 100
  layer_sizes      = [100, 100, 100]
```

These should be run as separate labels until stale, then compared against the
current BOLD P4 mainlines using best-checkpoint validation loss and downstream
pipeline7 usability.

## Raw-Coordinate Counterpart

After the standardized branch, a direct raw-coordinate counterpart was run to
test whether similar solver4 parameters can recover the same training behavior
without time-wise standardization.

Run script:

- `tmp/run_bold_solver4_raw_param_branch_20260520.py`

Overlay plot:

- `tmp/plot_bold_solver4_raw_vs_stdT_20260520.py`
- `tmp/figures/bold_solver4_raw_vs_stdT_e10gb1_20260520.png`
- `tmp/figures/bold_solver4_raw_vs_stdT_e10gb1_20260520.pdf`

Run labels use the prefix:

```text
mlp_obs_bold_batch4_rawparam_<observable>_<variant>_20260520_e10gb1_projected_vlambda_<observable>
```

The raw-coordinate comparison tested:

| Raw Variant | Matched Standardized Variant |
|---|---|
| `kellylike_raw` | `kellylike` |
| `reg1e2_raw` | `reg1e2` |
| `lr1e3_raw` | `lr1e3` |
| `oldcap_raw` | `oldcap` |

Important caveat: raw-coordinate and standardized-coordinate losses are not on
the same numerical scale. The meaningful comparison is the curve behavior:
smoothness, spike frequency, whether the best checkpoint is near the end, and
whether final validation loss remains close to the best checkpoint.

### Raw `roi_mean`

| Variant | Best Val | Best Epoch | Final Val | Final / First |
|---|---:|---:|---:|---:|
| `kellylike_raw` | 3.020997 | 95 | 3.174333 | 0.000034 |
| `reg1e2_raw` | 1.314292 | 38 | 1.332185 | 0.000016 |
| `lr1e3_raw` | 3.020929 | 32 | 3.849339 | 0.000010 |
| `oldcap_raw` | 2.049550 | 12 | 2.270328 | 0.000038 |

The raw `roi_mean` runs can reduce loss dramatically from a very large first
epoch loss, but the curves have large transient spikes and are not as clean as
the standardized runs. The best raw `roi_mean` short run was `reg1e2_raw`, but
the standardized `reg1e2` branch had a cleaner low-scale curve and reached its
best near the end.

### Raw `gsvd100_ds`

| Variant | Best Val | Best Epoch | Final Val | Final / First |
|---|---:|---:|---:|---:|
| `kellylike_raw` | 288.007278 | 96 | 311.637924 | 0.133826 |
| `reg1e2_raw` | 379.192642 | 40 | 379.192642 | 0.162137 |
| `lr1e3_raw` | 427.479683 | 7 | 2675.343675 | 0.343685 |
| `oldcap_raw` | 7.369793 | 23 | 509.541912 | 0.034243 |

The raw `gsvd100_ds` result is much less stable. In particular, `oldcap_raw`
briefly reaches a low best checkpoint but then drifts badly by the final epoch.
This is the same best/final split pattern that made earlier raw-coordinate BOLD
P4 runs hard to interpret.

### Raw-vs-Standardized Takeaway

The raw-coordinate runs are not completely unable to learn, but they do not
reproduce the clean standardized training behavior. The strongest evidence is:

- `roi_mean` raw losses start at extremely large values and show large spikes;
  standardized curves are much smoother.
- `gsvd100_ds` raw runs have severe best/final separation, especially
  `oldcap_raw`, while standardized `oldcap` continues to improve through epoch
  40.
- The same high learning-rate settings that worked well after standardization
  are unreliable in raw coordinates.

This supports the working conclusion that time-wise standardization is a core
part of the stable BOLD P4 training regime, not just a cosmetic rescaling.
