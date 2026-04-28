# BOLD IO Refactor Archive

Date: 2026-04-27

## Summary

This note archives the cleanup and contract-tightening process for the BOLD
IO layer.

The main outcome is:

- `functions/io/load_bold_dataset.m` is now treated as the formal BOLD raw-data loader.
- The loader now follows the same high-level philosophy as `load_blp_dataset.m`:
  load selected raw sessions from the source directory, then concatenate them
  in time.
- The old script-local BOLD folder loader has been removed from the batch
  observable script.
- BOLD ROI labels are no longer guessed from many possible struct fields.
  The formal label field is now `roiTs{1,r}.name`.
- Exact semantic ROI names such as `eleHP` and `HP` are now defined in each
  dataset config under `cfg.bold.role_map`, then validated against the raw
  `roiTs` files during loading.

## Why This Refactor Was Needed

Before this cleanup, the BOLD IO path had two competing entry styles:

1. `functions/io/load_bold_dataset.m`
   - expected a single pre-concatenated MAT file
   - supported broader input shapes such as:
     - concatenated `roiTs`
     - top-level `.dat`
     - `roiTs.observables`
   - accepted several metadata fallbacks

2. `scripts/script_batch_build_bold_observables.m`
   - used its own local helper `load_bold_roits_folder_dataset`
   - loaded one raw `roits` MAT file per session
   - concatenated sessions inside the script
   - guessed ROI names from many candidate fields

This had several problems:

- the real production loader for batch BOLD work was not in `functions/io`
- the formal loader and the actual batch loader used different assumptions
- ROI labels were not tied to one canonical raw-data field
- semantic roles such as `eleHP` and `HP` were partially managed in scripts,
  not in configs

## User-Level Decisions That Drove The Refactor

The refactor was based on three explicit decisions:

1. BOLD should always be loaded from the raw session files.
   - No formal support for pre-concatenated BOLD MAT files in the main loader.
   - The loader should read each session from raw `roits` files and build the
     concatenated dataset itself.

2. The raw data matrix must come from `.dat`.
   - The formal BOLD loader should read `roiTs{1,r}.dat`.
   - Other input forms such as `observables` should not be part of the main
     contract.

3. ROI labels should come from one fixed field in the raw structure.
   - The formal ROI label field is `roiTs{1,r}.name`.
   - The code should not scan a large fallback list such as
     `label / roi_name / region_name / atlas_label / ...`.

Later, one more decision was made:

4. Exact semantic ROI names such as `eleHP` / `HP` should be hardcoded in
   the dataset config because the current datasets are considered consistent.
   - These exact names should live in `cfg.bold.role_map`.
   - The loader should still validate them against the current raw
     `roiTs{1,r}.name` values.

## Final BOLD Loader Contract

The formal contract for `functions/io/load_bold_dataset.m` is now:

- expected source layout:
  - `fullfile(cfg.raw_data_root, cfg.bold.data_subfolder)`
- default BOLD subfolder:
  - `cfg.bold.data_subfolder = 'roits'`
- per-session filename:
  - `<file_stem>_%04d_roits.mat`
- expected MAT variable:
  - `cfg.bold.input_var`
  - default: `roiTs`
- expected raw BOLD container:
  - `roiTs` must be a nonempty cell array of ROI structs
- required per-ROI fields:
  - `.dat`
  - `.dx`
  - `.name`

The formal loader now:

- collects included session IDs from `cfg.sessions`
- loads one `roits` file per session
- reads selected ROI data from `roiTs{1,r}.dat`
- reads per-session TR from `roiTs{1,r}.dx`
- reads ROI names from `roiTs{1,r}.name`
- checks that selected ROI structure is consistent across sessions
- returns one standard `D` struct for downstream preprocessing and observable
  construction

## Config Convention Introduced In This Refactor

For the currently active BOLD datasets, exact semantic ROI names are now
defined in the dataset configs:

```matlab
cfg.bold = struct();
cfg.bold.data_subfolder = 'roits';
cfg.bold.input_var = 'roiTs';
cfg.bold.role_map = struct();
cfg.bold.role_map.elehp = 'eleHP';
cfg.bold.role_map.hp = 'HP';
```

The current configs updated in this pass are:

- `configs/cfg_E10gb1.m`
- `configs/cfg_E10gH1.m`
- `configs/cfg_E10fV1.m`
- `configs/cfg_F12m01.m`

Important distinction:

- `cfg.bold.role_map.*` stores the exact semantic ROI name expected for that
  dataset.
- `cfg.bold.selected_region_names` stores the exact region names selected for
  a particular loading call.

The batch observable script now uses the role map to choose the current exact
region names for branches such as:

- `eleHP`
- `HP`
- `HP_svd100`
- `slow_band_power`

## Refactor Steps

### Step 1. Tighten `load_bold_dataset`

`functions/io/load_bold_dataset.m` was rewritten from a broad multi-format
loader into a raw-session loader.

Major removals:

- support for a single pre-concatenated `cfg.bold.input_file`
- support for `roiTs.observables`
- support for top-level struct `.dat` as a formal alternate path
- support for multiple session-metadata fallback schemes
- support for broad ROI-label guessing from many field names

Major additions:

- strict session-file discovery from `cfg.raw_data_root`
- strict use of `roiTs{1,r}.dat`
- strict use of `roiTs{1,r}.dx`
- strict use of `roiTs{1,r}.name`
- cross-session checks for:
  - variable layout
  - selected ROI structure
  - full ROI-name list

### Step 2. Remove the script-local BOLD folder loader

`scripts/script_batch_build_bold_observables.m` previously contained its own
local `load_bold_roits_folder_dataset` and helper functions.

That local loader path was removed.

The batch script now calls:

- `io_raw.load_bold_dataset(cfg_mode)`

This means the formal IO logic now lives in one place instead of being split
between `functions/io` and `scripts`.

### Step 3. Tighten ROI label handling in visualization / postprocessing

Several BOLD plotting scripts had their own `local_roi_label` helper with
fallback behavior.

Those helpers were tightened so that they now require `R.name` rather than
falling back silently to generated labels such as `roi%02d`.

### Step 4. Move semantic ROI naming into configs

There was a short intermediate step where the batch script tried to infer
exact `eleHP` / `HP` names from the current raw `roiTs` name list.

That approach was later replaced by a simpler config-driven rule:

- exact semantic ROI names are defined in `cfg.bold.role_map`
- the loader validates them against the current raw `roiTs`

This keeps the exact names explicit and stable while preserving raw-data
validation.

## Current Files Most Directly Affected

Core loader:

- `functions/io/load_bold_dataset.m`

Batch BOLD observable entry:

- `scripts/script_batch_build_bold_observables.m`

Configs updated with BOLD role map:

- `configs/cfg_E10gb1.m`
- `configs/cfg_E10gH1.m`
- `configs/cfg_E10fV1.m`
- `configs/cfg_F12m01.m`

Templates aligned to the new contract:

- `scripts/templates/script_load_bold_data.m`
- `scripts/templates/script_build_bold_observables.m`
- `scripts/templates/script_preprocess_bold_sessions.m`

Visualization scripts aligned to strict `R.name` behavior:

- `scripts/plot_bold_reskoopnet_mode_activation_map_reference_style.m`
- `scripts/script_generate_missing_bold_activation_maps.m`
- `scripts/script_run_one_bold_reskoopnet_post_xcorr_activation_maps.m`

## What Was Not Done

This pass did not try to solve everything about BOLD configuration.

Specifically not done:

- no attempt to fill `cfg.bold` for every dataset config in the repo
- no claim that every dataset has verified BOLD exact ROI names already
- no broad rework of BOLD plotting logic beyond label-field tightening
- no change to the main BOLD observable construction functions
  (`build_bold_observables`, `preprocess_bold_sessions`)

## Remaining Practical Risk

The main remaining practical risk is not in the code structure itself, but in
whether the hardcoded exact semantic ROI names in `cfg.bold.role_map` match the
real `roiTs{1,r}.name` values for each dataset.

The code is now designed to fail loudly if those names do not match.

That is intentional.

## Validation Status

This refactor was completed at the code level, but it was not fully validated
by running the raw BOLD pipeline end to end in the current agent environment.

Known environment limitations during this pass:

- the current shell environment could not access the expected `E:` raw-data drive
- MATLAB startup also failed in this environment

So the current status is:

- contract cleanup: done
- code-path consolidation: done
- end-to-end raw-data verification in this environment: not done

## Intended Reading Of This Note

This note is meant to answer:

- why BOLD IO was changed
- what assumptions were removed
- what the formal raw-data contract is now
- where exact semantic ROI naming now lives
- which files to inspect first if BOLD loading breaks later
