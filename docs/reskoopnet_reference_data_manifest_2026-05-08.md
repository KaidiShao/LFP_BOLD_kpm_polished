## ResKoopNet Reference Data Manifest

Date: `2026-05-08`

This note records the canonical data naming used for the local ResKoopNet
reproduction work. The goal is to separate:

- the `pendulum` data that exactly matches the official Google Drive
  `data_pendulum_90.mat`
- the distinct `pendulum` data that matches the current repository/notebook
  generation logic
- the `turbulence` pressure-field input and official reference MAT exports

### Pendulum

There are two distinct `pendulum` input datasets in this workspace.

1. Official Google Drive pendulum `90` data

- canonical alias:
  `tmp/reference_data/reskoopnet/pendulum/pendulum_input_legacy_original_90ic_rollout.mat`
- official-drive alias:
  `tmp/reference_data/reskoopnet/pendulum/pendulum_input_official_drive_90.mat`
- source file:
  `tmp/figure2d_lab/original_data/data_pendulum_90.mat`
- used by:
  the local pendulum reproduction line that produced the annulus-style
  pseudospectrum used as the current project baseline
- origin summary:
  exact byte-for-byte match to the official Google Drive
  `pendulum_data/data_pendulum_90.mat`
  (`SHA256=0f4ecf6e468a957f11163a7f644707a2d861ab2c6c592569f933a8a76c188dba`)
  and therefore the canonical official `90`-IC pendulum input for this project

2. Repo-style/regenerated pendulum data

- canonical alias:
  `tmp/reference_data/reskoopnet/pendulum/pendulum_input_repo_style_90ic_fine1000.mat`
- explicit regenerated alias:
  `tmp/reference_data/reskoopnet/pendulum/pendulum_input_repo_regenerated_90ic_fine1000.mat`
- source file:
  `tmp/figure2d_lab/data_pendulum_90.mat`
- duplicate copies with identical array contents:
  `results/pendulum_solver3_full/data/data_pendulum_90.mat`
  `results/pendulum_solver4_full/data/data_pendulum_90.mat`
- used by:
  repository/notebook-style pendulum experiments
- origin summary:
  does not match the official Drive `data_pendulum_90.mat`; instead it matches
  the current `pendulum.m` generation logic after setting the initial condition
  grid to `90` points:
  `grid_angle_count=10`, `grid_velocity_count=10`,
  `integration_intervals=1000`, `delta_t=0.5`

### Turbulence

The turbulence reference set is simpler and is treated as canonical directly.

1. Pressure-field input

- canonical alias:
  `tmp/reference_data/reskoopnet/turbulence/turbulence_input_pressure_data.mat`
- official-drive alias:
  `tmp/reference_data/reskoopnet/turbulence/turbulence_input_pressure_data_official_drive.mat`
- source file:
  `tmp/turbulence_lab/original_data/pressure_data.mat`
- used by:
  all local turbulence reproductions and active-validation runs
- origin summary:
  exact byte-for-byte match to the official Google Drive
  `turbulence_data/pressure_data.mat`
  (`SHA256=d131ebe6b82e86952cbb93b65e7b4cd34cd98c993b750b4d2d5b91278d34991f`)

2. Official turbulence ResKoopNet exports

- canonical aliases:
  `tmp/reference_data/reskoopnet/turbulence/turbulence_official_result_150basis.mat`
  `tmp/reference_data/reskoopnet/turbulence/turbulence_official_result_250basis.mat`
  `tmp/reference_data/reskoopnet/turbulence/turbulence_official_result_300basis.mat`
- source files:
  `tmp/turbulence_lab/original_data/turbulence_resdmd_150basis_official.mat`
  `tmp/turbulence_lab/original_data/turbulence_resdmd_250basis_official.mat`
  `tmp/turbulence_lab/original_data/turbulence_resdmd_300basis_official.mat`
- main comparison baseline:
  the `250basis` official export

### Recommended Usage

- For local `pendulum` pseudospectrum reproduction comparisons:
  use `pendulum_input_official_drive_90.mat`
- For upstream/repository-style `pendulum` checks:
  use `pendulum_input_repo_regenerated_90ic_fine1000.mat`
- For turbulence:
  use `turbulence_input_pressure_data_official_drive.mat` plus
  `turbulence_official_result_250basis.mat`
