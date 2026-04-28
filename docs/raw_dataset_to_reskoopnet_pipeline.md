# From Config To ResKoopNet

This note records the main MATLAB call sequence for taking one raw BLP dataset
from a `cfg_*.m` file to saved ResKoopNet observables.

Scope:

- start from a dataset config such as `cfg_F12m01()` or `cfg_E10gb1()`
- load raw selected BLP sessions
- compute region-mean spectrograms
- build ResKoopNet observables
- stop at saved dictionary files and observable metadata CSVs

This note does not cover downstream Python training.

## Required functions

Main functions:

- `load_blp_dataset`
- `compute_blp_region_spectrograms`
- `compute_blp_region_spectrograms_streamed`
- `build_reskoopnet_dicts`

Main config inputs:

- `configs/cfg_F12m01.m`
- `configs/cfg_E10gb1.m`

Helpful wrapper scripts:

- `scripts/script_build_reskoopnet_dicts.m`
- `scripts/script_preprocess_e10gb1_to_observables_streamed.m`

## Prerequisites

Before running the pipeline, make sure MATLAB can see:

- this repo
- EEGLAB utilities that provide `finputcheck`
- `timefreqMB`

Typical setup:

```matlab
repo_root = 'D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished';
addpath(genpath(repo_root));

eeglab_root = 'D:\Onedrive\Toolbox\eeglab10_2_5_8b\';
addpath(genpath(eeglab_root));

assert(exist('timefreqMB', 'file') == 2, 'timefreqMB not found.');
assert(exist('finputcheck', 'file') == 2, 'finputcheck not found.');
```

## Step 1: Create Or Select A Config

Start with one dataset config:

```matlab
cfg = cfg_F12m01();
% or
cfg = cfg_E10gb1();
```

The config needs to define at least:

- `cfg.dataset_id`
- `cfg.raw_data_root`
- `cfg.data_subfolder`
- `cfg.file_stem`
- `cfg.channels.sites`
- `cfg.channels.selected_labels`
- `cfg.channels.selected_all`
- `cfg.sessions(i).session_id`
- `cfg.sessions(i).include`
- `cfg.sessions(i).selected_channels`

Important assumptions used later in the pipeline:

- raw files are named like `file_stem_0001_blp.mat`
- each file contains `blp.dat` and `blp.dx`
- `load_blp_dataset` uses `blp.dat(:, selected_channels, 1)`
- concatenated sessions must use the same selected channels

## Step 2: Load Raw BLP Data

Call:

```matlab
D = io_raw.load_blp_dataset(cfg);
```

Main output fields in `D`:

- `D.data`
- `D.session_ids`
- `D.session_lengths`
- `D.session_dx`
- `D.session_start_idx`
- `D.session_end_idx`
- `D.border_idx`
- `D.selected_channels`

Notes:

- `D.data` is concatenated as `time x selected_channel`
- if `dx` differs across sessions, the function warns and leaves `D.dx` and `D.fs` empty

## Step 3: Compute Region-Mean Spectrograms

The standard path is the streamed implementation:

```matlab
output_root = io_project.get_project_processed_root();

spec_opts = struct();
spec_opts.save_precision = 'single';
spec_opts.return_data = false;
spec_opts.force_recompute = false;

S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, spec_opts);
```

If you explicitly want the full spectrogram arrays loaded back into memory
for a small run, the old convenience entry now wraps the streamed path:

```matlab
output_root = io_project.get_project_processed_root();
S = compute_blp_region_spectrograms(D, cfg, output_root);
```

Recommended spectrogram settings can be placed in the config before this step:

```matlab
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';
```

Main saved outputs:

- region-mean abs spectrogram MAT file
- region-mean complex spectrogram MAT file

The standard variable names expected downstream are:

- `tmpall_mean_abs`
- `tmpall_mean_complex`

When to prefer the streamed version:

- many sessions
- very long recordings
- `compute_blp_region_spectrograms` would OOM at final concatenation

## Step 4: Build ResKoopNet Observables

After spectrograms exist on disk, call `build_reskoopnet_dicts`.

Common parameter template:

```matlab
params = struct();
params.low_full_max_hz = 50;
params.high_max_hz = 250;
params.high_group_size = 2;
params.chunk_size = 200000;
params.precision = 'single';
```

### 4A. Build `abs` observable

```matlab
params.spec_mode = 'abs';
dict_abs = build_reskoopnet_dicts(D, cfg, output_root, params);
```

### 4B. Build `complex_split` observable

```matlab
params.spec_mode = 'complex_split';
dict_complex = build_reskoopnet_dicts(D, cfg, output_root, params);
```

What this function does:

- keeps raw BLP channels as observables
- reads saved region-mean spectrograms from disk
- keeps all frequency bins up to `low_full_max_hz`
- averages high-frequency bins above that threshold in groups of `high_group_size`
- for `complex_split`, writes real and imaginary spectrogram parts as separate observables

Main saved outputs per mode:

- one MAT file with `obs`
- one CSV file with observable metadata

Typical output folder:

```text
<output_root>\<file_stem>\reskoopnet_dictionary\
```

Typical filenames:

- `<file_stem>_low50_high250_g2_abs_single.mat`
- `<file_stem>_low50_high250_g2_abs_single_obs_info.csv`
- `<file_stem>_low50_high250_g2_complex_split_single.mat`
- `<file_stem>_low50_high250_g2_complex_split_single_obs_info.csv`

## Minimal End-To-End Example

```matlab
repo_root = 'D:\Onedrive\ICPBR\Alberta\koopman_events\LFP_BOLD_kpm_polished';
addpath(genpath(repo_root));
addpath(genpath('D:\Onedrive\Toolbox\eeglab10_2_5_8b\'));

cfg = cfg_E10gb1();
cfg.spectrogram.pad_sec = 20;
cfg.spectrogram.pad_mode = 'mirror';

output_root = io_project.get_project_processed_root();

D = io_raw.load_blp_dataset(cfg);

spec_opts = struct();
spec_opts.save_precision = 'single';
spec_opts.return_data = false;
spec_opts.force_recompute = false;

S = compute_blp_region_spectrograms_streamed(D, cfg, output_root, spec_opts);

params = struct();
params.low_full_max_hz = 50;
params.high_max_hz = 250;
params.high_group_size = 2;
params.chunk_size = 200000;
params.precision = 'single';

params.spec_mode = 'abs';
dict_abs = build_reskoopnet_dicts(D, cfg, output_root, params);

params.spec_mode = 'complex_split';
dict_complex = build_reskoopnet_dicts(D, cfg, output_root, params);
```

## Existing Script Wrappers

If you prefer wrappers instead of calling functions manually:

- `scripts/script_build_reskoopnet_dicts.m`
  - simple example using `cfg_F12m01()`
- `scripts/script_preprocess_e10gb1_to_observables_streamed.m`
  - large-dataset example for `E10gb1`
  - runs raw load, streamed spectrograms, `abs`, and `complex_split`

## Practical Notes

- `build_reskoopnet_dicts` expects the spectrogram files to match the same
  `cfg.file_stem`, `pad_mode`, and `pad_sec`
- if you change spectrogram padding settings, regenerate the spectrogram files
  before rebuilding observables
- if session sampling periods differ slightly, preprocessing may still finish,
  but downstream event analyses should be checked more carefully
- for large datasets, `compute_blp_region_spectrograms_streamed` is the safer default
