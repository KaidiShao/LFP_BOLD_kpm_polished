# BLP Dictionary Pipeline: Code Structure, Role, And Calling Mechanism

This note records the current code organization for the first major pipeline in
this repo:

- raw BLP
- region-mean spectrograms
- ResKoopNet dictionaries

The goal is to make the code structure easy to understand without reading every
script and helper file.

## Scope

This document covers only the BLP dictionary line:

1. load raw BLP sessions
2. compute region-mean spectrograms
3. build ResKoopNet dictionary files

It stops at saved dictionary MAT files and observable metadata CSV files.

It does not cover:

- event/state processing
- Python ResKoopNet training
- downstream eigenfunction postprocessing

## Pipeline Shape

The current canonical flow is:

```text
cfg_*.m
  -> io_raw.load_blp_dataset
  -> compute_blp_region_spectrograms_streamed
  -> build_reskoopnet_dicts
  -> saved dictionary MAT + obs_info CSV
```

This pipeline is intentionally script-driven at the top level, so intermediate
variables such as `cfg`, `D`, `S`, and `dict` remain visible in the MATLAB
workspace.

## Canonical Entry Scripts

There are now three different entry scripts with different roles.

### 1. Canonical interactive single-dataset entry

Use:

- `scripts/script_preprocess_one_cfg_to_observables_streamed.m`

Role:

- primary interactive entry for one manually defined dataset config
- minimal top-to-bottom script
- keeps intermediate variables in workspace
- best default choice when switching to a new `cfg_*.m`

Typical usage:

```matlab
cfg_name = 'E10gH1';
run('scripts/script_preprocess_one_cfg_to_observables_streamed.m');
```

### 2. Batch driver for multiple datasets

Use:

- `scripts/script_preprocess_cfgs_to_observables_streamed.m`

Role:

- batch driver for multiple `cfg_*.m` definitions
- adds loop, optional logging, and lightweight run-management options
- less minimal than the single-dataset script, but better for multi-dataset runs

Typical usage:

```matlab
cfg_names = {'E10fV1', 'E10gH1'};
run('scripts/script_preprocess_cfgs_to_observables_streamed.m');
```

### 3. Single-dataset example/reference script

Use:

- `scripts/script_preprocess_e10gb1_to_observables_streamed.m`

Role:

- readable fixed example for `E10gb1`
- useful as a reference pipeline page
- not the preferred reusable entry

## Core Function Ownership

### Raw BLP loading

File:

- `functions/io/+io_raw/load_blp_dataset.m`

Role:

- load selected sessions from raw BLP MAT files
- concatenate sessions
- attach session boundary metadata
- optionally support metadata-only or disk-backed loading modes

Main output:

- `D`

Important output fields:

- `D.data`
- `D.n_time`
- `D.selected_channels`
- `D.session_ids`
- `D.session_lengths`
- `D.session_dx`
- `D.session_start_idx`
- `D.session_end_idx`
- `D.border_idx`

### Spectrogram producer

Files:

- `functions/preprocessing/compute_blp_region_spectrograms_streamed.m`
- `functions/preprocessing/compute_blp_region_spectrograms.m`

Role split:

- `compute_blp_region_spectrograms_streamed` is the canonical implementation
- `compute_blp_region_spectrograms` is a thin in-memory wrapper around the
  streamed implementation

Canonical behavior:

- compute spectrograms session by session
- average channels within each selected region
- save two MAT files:
  - abs spectrogram
  - complex spectrogram

Main output:

- `S`

Important output fields:

- `S.abs_file`
- `S.complex_file`
- `S.tmpall_mean_abs_size`
- `S.tmpall_mean_complex_size`

### Shared spectrogram file contract

File:

- `functions/preprocessing/resolve_regionmean_spectrogram_files.m`

Role:

- single source of truth for the saved region-mean spectrogram file contract
- defines:
  - `pad_sec`
  - `pad_mode`
  - `pad_tag`
  - canonical save directory
  - abs/complex filenames
  - abs/complex variable names

This helper is shared by:

- `compute_blp_region_spectrograms_streamed`
- `build_reskoopnet_dicts`
- `script_rerun_e10gh1_cfg_to_observables_and_events.m`

This avoids repeating spectrogram path construction in multiple places.

### Dictionary builder

File:

- `functions/preprocessing/build_reskoopnet_dicts.m`

Role:

- consume raw BLP plus saved region-mean spectrogram files
- keep raw BLP channels as observables
- append spectrogram-based observables
- save one dictionary MAT file and one observable-info CSV file

Supported main modes:

- `abs`
- `complex_split`

Important design point:

- feature layout is now defined once through one feature plan
- that same plan drives both:
  - column construction in `obs`
  - metadata rows in `obs_info`

Main output:

- `dict`

Important output fields:

- `dict.save_file`
- `dict.info_csv_file`
- `dict.n_time`
- `dict.n_obs`

## Calling Mechanism

### Interactive single-dataset calling mechanism

The canonical single-dataset script works like this:

1. user sets `cfg_name`
2. script resolves `cfg_fun = ['cfg_' cfg_name]`
3. script evaluates the config
4. script applies spectrogram defaults
5. script calls `io_raw.load_blp_dataset`
6. script calls `compute_blp_region_spectrograms_streamed`
7. script loops over requested dictionary modes
8. script calls `build_reskoopnet_dicts` once per mode

In short:

```matlab
cfg_name -> cfg -> D -> S -> dict_outputs
```

### Batch calling mechanism

The batch driver adds one outer loop:

```matlab
cfg_names{1}
cfg_names{2}
...
    -> cfg
    -> D
    -> S
    -> dict per mode
```

The core stage calls do not change. The batch script only adds:

- looping over datasets
- optional logging
- optional raw-data caching controls

## Required Config Contract

For a new dataset to work in this pipeline, its `cfg_*.m` must match the
current BLP loader contract.

At minimum, it should define:

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

The raw BLP files are assumed to follow the current naming and storage
contract:

- file name like `<file_stem>_0001_blp.mat`
- MAT file contains `blp.dat`
- MAT file contains `blp.dx`

## Main User-Facing Variables

When the pipeline is run from a script, the main visible workspace variables
are:

- `cfg`
- `output_root`
- `D`
- `S`
- `params`
- `dict`
- `dict_outputs`

This is intentional. The top level remains a script because it is useful to
inspect these values interactively.

## Output Layout

### Spectrogram outputs

Saved under:

- `<output_root>\<file_stem>\spectrograms\`

Canonical files:

- `<file_stem><pad_tag>_regionmean_spectrograms_abs.mat`
- `<file_stem><pad_tag>_regionmean_spectrograms_complex.mat`

### Dictionary outputs

Saved under:

- `<output_root>\<file_stem>\reskoopnet_dictionary\`

Typical files:

- `<save_tag>.mat`
- `<save_tag>_obs_info.csv`

## Current Code Health

Current status of this pipeline:

- the main flow is clear
- the canonical entry points are now separated by role
- the spectrogram filename contract is now defined in one place
- the stage functions are already fairly direct

Remaining complexity is mostly real artifact complexity rather than structural
noise.

## Practical Recommendation

When working on this pipeline:

- use `script_preprocess_one_cfg_to_observables_streamed.m` for one new dataset
- use `script_preprocess_cfgs_to_observables_streamed.m` only when you really
  want a batch run
- use `script_preprocess_e10gb1_to_observables_streamed.m` as a readable
  reference example
- keep stage logic in functions
- keep pipeline wiring in scripts
